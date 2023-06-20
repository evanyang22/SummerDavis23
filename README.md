# SummerDavis23

## miRNA alignment using miRDeep2
Functions of mirDeep2 (https://github.com/rajewsky-lab/mirdeep2)
1. Mapper.pl : processes reads and maps them to reference genome, is a pre-requisite for miRDeep2.pl and quantifier.pl
2. miRDeep2.pl :  detects novel microRNAs
3. Quantifier.pl : maps and quantifiers microRNAs based on microRNA precursors

#### Steps 

## totalRNA alignment using STAR


---------------------




Index builds using bowtie-build: $ bowtie-build reference_sequence.fasta index_name

Generate clean .fa files using perl: perl -plane 's/\s+.+$//' < genome.fa > new_genome.fa

https://www.biostars.org/p/147955/ - replace unknown characters with N

Use VIM and this command: :%s/R/N/g 
