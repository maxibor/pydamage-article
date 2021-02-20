####################################################################################################
# Determine the presence of ancient DNA damage using pyDamage of randomly drawn contigs with respect
# of the contig length, GC content of the contig, coverage along the contig, read length
# distribution, and levels of ancient DNA damage.
#
# Requires the output files of simulation_shortreaddata.Snakefile
#
# Alex Huebner, 15/02/21
####################################################################################################


from glob import glob
import os

import pandas as pd

#### EXPERIMENT VARIABLES ##########################################################################
GENOMES, = glob_wildcards("genomes/{genome}.fna.gz")
# Genome lengths
GENOMELENGTHS = pd.read_csv("results/genomelength.txt", sep="\t") \
    .rename({'sample': 'genome'}, axis=1) \
    .set_index(['genome'])
# Chromosome lengths
CHROMLENGTHS = pd.read_csv("results/chromlength.txt", sep="\t", index_col=['genome'])
# Read length profiles
READLENGTHPROFILES = ['short', 'medium', 'long']
# Number of simulated reads per read length profile
NREADS = {'short': 3e7,
          'medium': 2.1e7,
          'long': 1.5e7}
# Read length distribution used for simulation in each profile
READLENGTHDISTS = pd.concat([pd.read_csv(fn, sep="\t",
                                         header=None, names=['length', 'frac']) \
                                 .assign(rl=os.path.basename(fn).split("_")[0])
                             for fn in glob("results/*_sizefreq.size.gz")]) \
    .set_index(['rl'])
# Coverage bins
COVBINS = [(1, 2), (2, 3), (3, 5), (5, 10), (10, 20), (20, 50), (50, 100), (100, 200), (200, 500)]
# Contig length bins
CLBINS = [(500, 1000), (1000, 2000), (2000, 5000),
          (5000, 10000), (10000, 20000), (20000, 50000),
          (50000, 100000), (100000, 200000), (200000, 500000)]
NUMBER_CLBINS = {'Adentalis': 8,
                 'Msmithii': 9,
                 'Tforsythia': 9}
# Damage levels
DAMAGE = {level: amount
          for level, amount in zip(list(range(0, 10)), [0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1, 0.15, 0.2])}
####################################################################################################

#### Auxilliary functions ##########################################################################

def draw_contig(chrlengths, genome, b, seed=[0, 0, 0]):
    """ Draw a contig of a contig length bin and return coordinates in BED format.

        chrlengths : Pandas DataFrame with all individual contig lengths per reference genome
        genome     : name of the reference genome
        b          : tuple with the boundaries of the contig lengths
        seed       : tuple with the three seeds for each of random number draws
    """
    # Draw contig length
    np.random.seed(seed=seed[0])
    if b[1] <= chrlengths.loc[genome, 'length'].max():
        cl = np.random.randint(b[0], b[1], 1)[0]
    else:
        cl = np.random.randint(b[0], chrlengths.loc[genome, 'length'].max(), 1)[0]
    # Subset to contigs that are long enough
    chromosomes = chrlengths.loc[chrlengths['length'] >= cl].loc[[genome]]
    np.random.seed(seed=seed[1])
    chridx = np.random.choice(chromosomes.shape[0], 1, replace=False)[0]
    np.random.seed(seed=seed[2])
    start = np.random.randint(0, chromosomes.iloc[chridx]['length'] - cl, 1)[0]
    return f"{chromosomes.iloc[chridx]['chr']}:{start}-{start + cl}"


def subsampling_fraction(covbin, rlp, n, genomelength, seed=0):
    """ Determines the fraction required for subsampling to a specific coverage.

        Parameters:
        covbin       : tuple with interval for final coverage along reference genome
        rlp          : undelrying read length profile used for simulation
        n            : number of simulated reads
        genomelength : length of reference genome
        seed         : seed used for random state initialisation
    """
    np.random.seed(seed)
    coverage = np.random.uniform(covbin[0], covbin[1], size=1)[0]
    totalbases = (rlp['length'] * rlp['frac']).sum() * n
    estbases = coverage * genomelength
    return round(seed + (estbases / totalbases), 6)

####################################################################################################


rule all:
    input: 
        expand(f"{workflow.basedir}/../05-results/pydamage_results/{genome}-{rl}.pydamage.tsv.gz", genome=['Adentalis', 'Msmithii', 'Tforsythia'], rl=READLENGTHPROFILES)

rule uncompress_fasta:
    output:
        "contigs/{genome}.fna"
    message: "Decompress FastA file of {wildcards.genome}"
    params:
        fasta = "genomes/{genome}.fna.gz"
    shell:
        "zcat {params.fasta} > {output}"

rule write_region_to_bed:
    output:
        temp("pydamage/contigs/{genome}/{genome}-{clbin}-{c}.bed")
    message: "Write region to BED file for contig {wildcards.c} of contig length bin {wildcards.clbin} for genome {wildcards.genome}"
    params:
        region = lambda wildcards: draw_contig(CHROMLENGTHS, wildcards.genome, CLBINS[int(wildcards.clbin)], seed=[int(wildcards.c)] * 3)
    shell:
        """
        echo "{params.region}" | sed "s/[:\-]/\t/g" > {output}
        """

rule extract_fa:
    input:
        fasta = "contigs/{genome}.fna",
        bed = "contigs/{genome}/{genome}-{clbin}-{c}.bed"
    output:
        "contigs/{genome}/{genome}-{clbin}-{c}.fasta"
    message: "Extract FastA sequence of contig {wildcards.c} of contig length bin {wildcards.clbin} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fasta} -bed {input.bed} > {output}
        """

rule gccontent:
    input:
        "contigs/{genome}/{genome}-{clbin}-{c}.fasta"
    output:
        temp("contigs/{genome}/{genome}-{clbin}-{c}.gccontent")
    message: "Determine GC content of contig {wildcards.c} of contig length bin {wildcards.clbin} for genome {wildcards.genome}"
    shell:
        """
        bioawk -c fastx '{{
            print $name "\t" gc($seq) "\t" length($seq)
        }}' {input} > {output}
        """

rule build_bwa_index:
    input:
        "contigs/{genome}/{genome}-{clbin}-{c}.fasta"
    output:
        amb = "contigs/{genome}/{genome}-{clbin}-{c}.amb",
        ann = "contigs/{genome}/{genome}-{clbin}-{c}.ann",
        bwt = "contigs/{genome}/{genome}-{clbin}-{c}.bwt",
        pac = "contigs/{genome}/{genome}-{clbin}-{c}.pac",
        sa = "contigs/{genome}/{genome}-{clbin}-{c}.sa"
    message: "Build index for BWA of contig {wildcards.c} of contig length bin {wildcards.clbin} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        prefix = "contigs/{genome}/{genome}-{clbin}-{c}"
    shell:
        """
        bwa index -p {params.prefix} {input}
        """

rule extract_and_subsample:
    output:
        temp("pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.fq.gz")
    message: "Subsample repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} and extract contig {wildcards.c} of contig length bin {wildcards.clbin}"
    conda: "simulation.yaml"
    params:
        bam = "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sorted.bam",
        region = lambda wildcards: draw_contig(CHROMLENGTHS, wildcards.genome, CLBINS[int(wildcards.clbin)], seed=[int(wildcards.c)] * 3),
        frac = lambda wildcards: subsampling_fraction(COVBINS[int(wildcards.covbin)], READLENGTHDISTS.loc[wildcards.rl], NREADS[wildcards.rl], GENOMELENGTHS.at[wildcards.genome, 'genomelength'], seed=int(wildcards.c))
    shell:
        """
        samtools view -bh -s {params.frac} {params.bam} {params.region} | \
        samtools fastq - | \
        bioawk -c fastx '{{
            gsub(/!/, "F", $qual);
            print "@" $name " " $comment;
            print $seq;
            print "+";
            print $qual
        }}' - | gzip > {output}
        """

rule bwa_aln:
    input:
        fq = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.fq.gz",
        amb = "contigs/{genome}/{genome}-{clbin}-{c}.amb",
        ann = "contigs/{genome}/{genome}-{clbin}-{c}.ann",
        bwt = "contigs/{genome}/{genome}-{clbin}-{c}.bwt",
        pac = "contigs/{genome}/{genome}-{clbin}-{c}.pac",
        sa = "contigs/{genome}/{genome}-{clbin}-{c}.sa"
    output:
        temp("pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sai")
    message: "Align data of repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} and extract contig {wildcards.c} of contig length bin {wildcards.clbin} using BWA aln"
    conda: "simulation.yaml"
    params:
        prefix = "contigs/{genome}/{genome}-{clbin}-{c}"
    threads: 4
    shell:
        """
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -f {output} \
            {params.prefix} \
            {input}
        """

rule bwa_samse:
    input:
        fq = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.fq.gz",
        amb = "contigs/{genome}/{genome}-{clbin}-{c}.amb",
        ann = "contigs/{genome}/{genome}-{clbin}-{c}.ann",
        bwt = "contigs/{genome}/{genome}-{clbin}-{c}.bwt",
        pac = "contigs/{genome}/{genome}-{clbin}-{c}.pac",
        sa = "contigs/{genome}/{genome}-{clbin}-{c}.sa",
        sai = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sai"
    output:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam"
    message: "Generate alignment file of repeat {wildcards.rep} of genome {wildcards.genome} and read length profile {wildcards.rl} with damage level {wildcards.damage}"
    conda: "simulation.yaml"
    params:
        prefix = "contigs/{genome}/{genome}-{clbin}-{c}"
    shell:
        """
        bwa samse \
            {params.prefix} \
            {input.sai} \
            {input.fq} | \
        samtools view -bh - |
        samtools sort -o {output} -
        """

rule bam_index:
    input:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam"
    output:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam.bai"
    message: "Index BAM file of the repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} and extract contig {wildcards.c} of contig length bin {wildcards.clbin}"
    shell:
        """
        samtools index {input}
        """

rule pydamage:
    input:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam.bai"
    output:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}/pydamage_results.csv"
    message: "Run pyDamage on the repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} and extract contig {wildcards.c} of contig length bin {wildcards.clbin}"
    conda: "simulation.yaml"
    params:
        bam = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam",
        dir = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}",
        region = lambda wildcards: draw_contig(CHROMLENGTHS, wildcards.genome, CLBINS[int(wildcards.clbin)], seed=[int(wildcards.c)] * 3)
    threads: 2
    shell:
        """
        if [[ $(samtools view -c -F 4 {params.bam}) -gt 0 ]]; then
            pydamage -w 35 \
                    -p {threads} \
                    -vv \
                    --force \
                    -o {params.outputprefix} {input.bam}
        fi
        if [[ ! -f {output} ]]; then
        echo -e
        "reference,pred_accuracy,null_model_p0,null_model_p0_stdev,damage_model_p,damage_model_p_stdev,damage_model_pmin,damage_model_pmin_stdev,damage_model_pmax,damage_model_pmax_stdev,pvalue,qvalue,RMSE,nb_reads_aligned,coverage,reflen,CtoT-0,CtoT-1,CtoT-2,CtoT-3,CtoT-4,GtoA-0,GtoA-1,GtoA-2,GtoA-3,GtoA-4" > {output}
        echo -e "{params.region},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,$(samtools view -c {params.bam}),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> {output}
        fi
        """

rule depth:
    input:
        "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam.bai"
    output:
        temp("pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.depth")
    message: "Determine depth of the repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} and extract contig {wildcards.c} of contig length bin {wildcards.clbin}"
    params:
        bam = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.sorted.bam",
    shell:
        """
        samtools depth -a {params.bam} | bioawk -t '{{d += $3}}END{{if (NR > 0){{print d / NR}} else {{print 0}}}}' > {output}
        """

rule summary:
    input:
        pydamage = expand("pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}/pydamage_results.csv", genome=config['genome'], rl=config['rl'], damage=config['damage'], rep=config['rep'], covbin=config['covbin'], clbin=config['clbin'], c=list(range(100))),
        depth = expand("pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}-{c}.depth", genome=config['genome'], rl=config['rl'], damage=config['damage'], rep=config['rep'], covbin=config['covbin'], clbin=config['clbin'], c=list(range(100)))
        gc = expand("contigs/{genome}/{genome}-{clbin}-{c}.gccontent", genome=config['genome'], covbin=config['covbin'], clbin=config['clbin'], c=list(range(100))),
    output:
        "pydamage_results/{genome}-{rl}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.pydamage.tsv"
    message: "Summarise results of the repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} subsampled to coverage bin {wildcards.covbin} looking at contigs of contig length bin {wildcards.clbin}"
    run:
        res = []
        for contig in glob(f"pydamage/{wildcards.genome}-{wildcards.rl}.{wildcards.damage}-{wildcards.rep}.{wildcards.covbin}-{wildcards.clbin}/{wildcards.genome}-{wildcards.rl}.{wildcards.damage}-{wildcards.rep}.{wildcards.covbin}-{wildcards.clbin}-*.depth"):
            mean_depth = float(next(open(contig, "rt")).rstrip())
            pydamage_res = pd.read_csv(contig.replace(".depth", "/pydamage_results.csv"), sep=",") \
                    .assign(genome=wildcards.genome) \
                    .assign(readlength=wildcards.rl) \
                    .assign(damage=DAMAGE[int(wildcards.damage)]) \
                    .assign(repeat=wildcards.rep) \
                    .assign(simuCov="-".join([str(i) for i in COVBINS[int(wildcards.covbin)]])) \
                    .assign(actualCov=mean_depth) \
                    .assign(simuContigLength="-".join([str(i) for i in CLBINS[int(wildcards.clbin)]]))
            for trans in ['GtoA-0', 'GtoA-1', 'GtoA-2', 'GtoA-3', 'GtoA-4']:
                if trans not in pydamage_res.columns:
                    pydamage_res[trans] = 0
            res.append(pydamage_res[['genome', 'readlength', 'damage', 'repeat',
                                     'simuCov', 'simuContigLength',
                                     'reference', 'null_model_p0', 'null_model_p0_stdev',
                                     'damage_model_p', 'damage_model_p_stdev',
                                     'damage_model_pmin', 'damage_model_pmin_stdev',
                                     'damage_model_pmax', 'damage_model_pmax_stdev',
                                     'pvalue', 'qvalue', 'RMSE', 'nb_reads_aligned',
                                     'coverage', 'actualCov',
                                     'CtoT-0', 'CtoT-1', 'CtoT-2', 'CtoT-3', 'CtoT-4',
                                     'GtoA-0', 'GtoA-1', 'GtoA-2', 'GtoA-3', 'GtoA-4']])
        pd.concat(res).to_csv(output[0], sep="\t", index=False, na_rep = "NA")
            
rule tarball_pdfs:
    input:
        "pydamage_results/{genome}-{rl}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.pydamage.tsv"
    output:
        "pydamage_results/{genome}-{rl}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.plots.tar.bz2"
    message: "Generate tarball from plots of the repeat {wildcards.rep} from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} subsampled to coverage bin {wildcards.covbin} looking at contigs of contig length bin {wildcards.clbin}"
    params:
        prefix = "pydamage/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}"
    shell:
        """
        find {params.prefix} -name "*.png" | tar -cvjf {output} --files-from -
        """


#rule write_config:
    #output:
        #temp("04-analysis/pydamage/configs/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.json")
    #message: "Write config file for {wildcards.genome}-{wildcards.rl}.{wildcards.damage}-{wildcards.rep}.{wildcards.covbin}-{wildcards.clbin}."
    #run:
        #config_dict = {'prefix': "04-analysis/pydamage",
                       #'genome': wildcards.genome,
                       #'rl': wildcards.rl,
                       #'damage': wildcards.damage,
                       #'rep': wildcards.rep,
                       #'covbin': wildcards.covbin,
                       #'clbin': wildcards.clbin,
                       #}
        #with open(output[0], "wt") as outfile:
            #json.dump(config_dict, outfile)

#rule run_pydamage_pipeline:
    #input:
        #"04-analysis/pydamage/configs/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.json"
    #output:
        #"04-analysis/pydamage/pydamage_results/{genome}-{rl}/{genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.pydamage.tsv"
    #message: "Run Pydamage pipeline on simulated data with parameters {wildcards.genome}-{wildcards.rl}.{wildcards.damage}-{wildcards.rep}.{wildcards.covbin}-{wildcards.clbin}."
    #params:
        #Snakefile = f"{workflow.basedir}/PYD_pydamage_simudata_pipeline.Snakefile"
    #threads: 24
    #shell:
        #"""
        #snakemake -s {params.Snakefile} \
                #--configfile {input} \
                #--cores {threads} \
                #--nolock \
                #-k --latency-wait 15
        #"""

#rule summarise_reports:
    #input:
        #lambda wildcards: [f"04-analysis/pydamage/pydamage_results/{wildcards.genome}-{rl}/{wildcards.genome}-{rl}.{damage}-{rep}.{covbin}-{clbin}.pydamage.tsv" for rl in READLENGTHPROFILES for damage in DAMAGE.keys() for rep in [0] for covbin in range(len(COVBINS)) for clbin in NUMBER_CLBINS[wildcards.genome]]
    #output:
        #f"{workflow.basedir}/../05-results/pydamage_results/{genome}-{rl}.pydamage.tsv.gz"
    #message: "Summarise reports of genome {wildcards.genome} and read length profile {wildcards.rl}"
    #run:
        #reports = [pd.read_csv(f"04-analysis/pydamage/pydamage_results/{wildcards.genome}-{wildcards.rl}/{wildcards.genome}-{wildcards.rl}.{damage}-0.{covbin}-{clbin}.pydamage.tsv",
                               #sep="\t")
                   #for damage in DAMAGE.keys()
                   #for covbin in range(len(COVBINS))
                   #for clbin in range(NUMBER_CLBINS[wildcards.genome])]
        #pd.concat(reports).to_csv(output[0], sep="\t", index=False, compression="gzip", na_rep="NA", float_format="%.5f")
