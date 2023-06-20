# SummerDavis23

## miRNA alignment using miRDeep2

## totalRNA alignment using STAR


---------------------




Index builds using bowtie-build: $ bowtie-build reference_sequence.fasta index_name

Generate clean .fa files using perl: perl -plane 's/\s+.+$//' < genome.fa > new_genome.fa

https://www.biostars.org/p/147955/ - replace unknown characters with N

Use VIM and this command: :%s/R/N/g 
