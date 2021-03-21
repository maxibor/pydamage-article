# Kraken version 2.1.1

all: kraken2_damaged kraken2_non_damaged
	@echo "Done"
    
kraken2_damaged:
	kraken2 --db /projects1/users/borry/02_db/kraken2/kraken2_plus_pfp_25_01_2021 \
	--threads 8 --output ../../results/kraken2/ZSM028_kraken2_damaged.out \
	--report-minimizer-data --report ../../results/kraken2/ZSM028_kraken2_damaged.kreport \
	../../data/ZSM028.pydamage_filtered.fa
	cut -f1-3,6-8 ../../results/kraken2/ZSM028_kraken2_damaged.kreport > ../../results/kraken2/ZSM028_kraken2_damaged.old_kreport

kraken2_non_damaged:
	kraken2 --db /projects1/users/borry/02_db/kraken2/kraken2_plus_pfp_25_01_2021 \
	--threads 8 --output ../../results/kraken2/ZSM028_kraken2_non_damaged.out \
	--report-minimizer-data --report ../../results/kraken2/ZSM028_kraken2_non_damaged.kreport \
	../../data/ZSM028.pydamage_not_filtered.fa
	cut -f1-3,6-8 ../../results/kraken2/ZSM028_kraken2_non_damaged.kreport > ../../results/kraken2/ZSM028_kraken2_non_damaged.old_kreport
    
clean:
	rm -rf ../../results/kraken2/*.kreport ../../results/kraken2/*.old_kreport ../../results/kraken2/*.out
