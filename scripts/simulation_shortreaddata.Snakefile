####################################################################################################
# Simulation of short-read sequenicng data with ancient DNA damage for three microbial taxa 
#
# Alex Huebner, 15/02/21
####################################################################################################

from glob import glob
import os
import re

from scipy.stats import lognorm
import numpy as np
import pandas as pd
import pyfastx

#### Microbial species #############################################################################
URLS = {'Msmithii': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1/GCF_000016525.1_ASM1652v1_genomic.fna.gz',
        'Tforsythia': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz',
        'Adentalis': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/429/225/GCF_000429225.1_ASM42922v1/GCF_000429225.1_ASM42922v1_genomic.fna.gz',
       }
####################################################################################################

wildcard_constraints:
    genome = "[A-Za-z0-9_]+",
    rl = "[a-z]+",
    damage = "[1-9]"

rule all:
    input: 
        "reference_genomes.done",
        "simulate.done"

#### Prepare reference genomes #####################################################################

rule reference_genomes:
    input:
        expand("genomes/{genome}.fna.gz.fai", genome=URLS.keys()),
        "results/gccontent.txt",
        "results/genomelength.txt",
        "results/chromlength.txt"
    output:
        touch("reference_genomes.done")

rule download_genomes:
    output:
        "genomes/{genome}.fna.gz"
    message: "Download genome of species {wildcards.genome}"
    params: 
        url = lambda wildcards: URLS[wildcards.genome]
    shell:
        """
        curl {params.url} | zcat - | bgzip > {output}
        """

rule index_fastas:
    input:
        "genomes/{genome}.fna.gz"
    output:
        "genomes/{genome}.fna.gz.fai"
    message: "Index downloaded FastA file of {wildcards.genome}"
    conda: "simulation.yaml"
    shell:
        "samtools faidx {input}"

rule determine_gc_content:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/gccontent.txt'
    message: "Calculate GC content of the downloaded genomes"
    conda: "simulation.yaml"
    shell:
        """
        echo -e 'sample\tGC' > {output}
        for g in {input}; do
            gc=$(bioawk -c fastx '{{
                    gcs = gcs + gc($seq)
                 }} END {{
                    print gcs / NR
                 }}' ${{g}})
            echo -e "$(basename ${{g}} .fna.gz)\t${{gc}}"
        done >> {output}
        """

rule determine_genome_length:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/genomelength.txt'
    message: "Calculate genome lengths of reference genomes"
    params:
        dir = 'genomes'
    run:
        with open(output[0], "wt") as outfile:
            outfile.write("sample\tgenomelength\n")
            for ref in list(URLS.keys()):
                length = 0
                for name, seq in pyfastx.Fasta(f"{params.dir}/{ref}.fna.gz", build_index=False):
                    length += len(seq)
                outfile.write(f"{ref}\t{length}\n")

rule determine_chrom_length:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/chromlength.txt'
    message: "Calculate the individual chromosome lengths of reference genomes"
    params:
        dir = 'genomes'
    run:
        chrom_lengths = pd.concat([pd.DataFrame([(name, len(seq))
                                                 for name, seq in pyfastx.Fasta(f"{params.dir}/{genome}.fna.gz",
                                                 build_index=False)], columns=['chr', 'length'])
                                     .assign(genome=genome)
                                   for genome in list(URLS.keys())])
        chrom_lengths[['genome', 'chr', 'length']] \
            .to_csv(output[0], sep="\t")

####################################################################################################

#### Simulate short-read sequencing data using gargammel ###########################################

# Parameters for simulations with gargammel
## read length profiles
READLENGTHS = {'short': (0.25, 75),  # sigma, scale for lognorm.rvs
               'medium': (0.35, 100),
               'long': (0.45, 150)}
## number of reads per read length profile 
NREADS = {'short': 30000000,
          'medium': 21000000,
          'long': 15000000}
## damage levels
DAMAGE = {level: amount
          for level, amount in zip(list(range(1, 10)), [0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1, 0.15, 0.2])}
# Briggs et al. aDNA model parameters for each profile
BRIGGSPARAMS = pd.DataFrame.from_dict({'readlength': ['short', 'medium', 'long'],
                                       'l': [0.5, 0.5, 0.2],
                                       'intercept': [0.00285351934051996, 0.006, 0.000806560431244065],
                                       's': [0.510462904248573, 0.46, 0.793827313143463]}).set_index('readlength')
# Repeats: 5
REPEATS = list(range(5))

## Auxilliary functions for gargammel
def infer_s(rl, ctfreq):
    """Infer the parameter s for the amount of C-to-T substitutions for the
       specified readlength profile."""
    p = BRIGGSPARAMS.loc[rl]
    return round((ctfreq - p.intercept) / p.s, 3)

rule simulate:
    input:
        expand("pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sorted.bam.bai", genome=URLS.keys(), rl=READLENGTHS.keys(), damage=DAMAGE.keys())
    output:
        touch("simulate.done")

rule generate_readlength_profile:
    output:
        "results/rldistribution_{rl}_sizefreq.size.gz"
    message: "Generate fragment size distribution for read length profile {wildcards.rl}"
    run:
        sigma, scale = READLENGTHS[wildcards.rl_profile]
        lengths_df = pd.DataFrame.from_dict({'length': lognorm.rvs(sigma, loc=-10,
                                                                   scale=scale,
                                                                   size=100000,
                                                                   random_state=1)})
        lengths_df['length'] = np.floor(lengths_df['length']).astype(int)
        lengths_df = lengths_df.query('length >= 35 & length < 300')
        frequency = lengths_df.groupby(['length'])[['length']].count() \
            .rename({'length': 'n'}, axis=1)
        frequency['freq'] = frequency['n'] / frequency['n'].sum()
        frequency[['freq']].to_csv(output[0], sep="\t", header=False,
                                   compression="gzip", float_format="%.6f")

rule fragSim:
    input:
        rl = expand("results/rldistribution_{rl}_sizefreq.size.gz", rl=READLENGTHS.keys()),
        refgenome = "reference_genomes.done"
    output:
        temp("gargammel/{genome}-{rl}.bam")
    message: "Simulate DNA molecules with read length distribution {wildcards.rl} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        rl_profile = "results/rldistribution_{rl}_sizefreq.size.gz",
        nreads = lambda wildcards: NREADS[wildcards.rl],
        genome = "genomes/{genome}.fna.gz"
    shell:
        """
        fragSim \
            -n {params.nreads} \
            -b {output} \
            -u \
            -m 35 \
            -M 300 \
            -f {params.rl_profile} \
            {params.genome} 
        """

rule deamSim:
    input:
        "gargammel/{genome}-{rl}.bam"
    output:
        "gargammel/{genome}-{rl}.{damage}.bam"
    message: "Simulate reads of read length profile {wildcards.rl} with damage level {wildcards.damage} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        s = lambda wildcards: infer_s(wildcards.rl, DAMAGE[int(wildcards.damage)]),
        l = lambda wildcards: BRIGGSPARAMS.loc[wildcards.rl, 'l']
    shell:
        """
        deamSim \
            -b {output} \
            -damage 0.03,{params.l},0.01,{params.s} \
            {input}
        """

rule bam2fq:
    input:
        "gargammel/{genome}-{rl}.{damage}.bam"
    output:
        temp("pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.fq.gz")
    message: "Convert BAM file of genome {wildcards.genome} and read length profile {wildcards.rl} with damage level {wildcards.damage} to FastQ"
    conda: "simulation.yaml"
    shell:
        """
        samtools fastq {input} > {output}
        """

rule bwa_aln:
    input:
        "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.fq.gz"
    output:
        temp("pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sai")
    message: "Align reads of genome {wildcards.genome} and read length profile {wildcards.rl} with damage level {wildcards.damage} to reference genome"
    conda: "simulation.yaml"
    params:
        bwaidx = "genomes/{genome}.fna.gz"
    threads: 8
    shell:
        """
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -f {output} \
            {params.bwaidx} \
            {input}
        """

rule bwa_samse:
    input:
        fq = "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.fq.gz",
        sai = "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sai"
    output:
        "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sorted.bam"
    message: "Generate alignment file of genome {wildcards.genome} and read length profile {wildcards.rl} with damage level {wildcards.damage}"
    conda: "simulation.yaml"
    params:
        bwaidx = "genomes/{genome}.fna.gz"
    shell:
        """
        bwa samse \
            {params.bwaidx} \
            {input.sai} \
            {input.fq} | \
        samtools view -bh - |
        samtools sort -o {output} -
        """

rule index_bam:
    input:
        "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sorted.bam"
    output:
        "pydamage/aligned_reads/{genome}/{genome}-{rl}.{damage}.sorted.bam.bai"
    message: "Create index for alignment file of genome {wildcards.genome} and read length profile {wildcards.rl} with damage level {wildcards.damage}"
    conda: "simulation.yaml"
    shell:
        """
        samtools index {input}
        """

####################################################################################################
