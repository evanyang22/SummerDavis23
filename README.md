# SummerDavis23

## miRNA alignment using miRDeep2
Functions of mirDeep2 (https://github.com/rajewsky-lab/mirdeep2)
1. Mapper.pl : processes reads and maps them to reference genome, is a pre-requisite for miRDeep2.pl and quantifier.pl
2. miRDeep2.pl :  detects novel microRNAs
3. Quantifier.pl : maps and quantifiers microRNAs based on microRNA precursors

#### Steps 

## totalRNA alignment using STAR

#### How to Use STAR General Instructions

Install trimgalore, fastqc, star,cutadapt
1. Download sheep reference genome from Ensembl database: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000298735.2/ Note: This can be any reference genome, also download the annotation files in GTF format
2. Activate cutadapt env : “conda activate cutadaptenv”
3. Run a quality control : “fastqc input.fastq”
4. Trim adapters and low quality bases from the input fastq file : “trim_galore --quality 20 --fastqc --illumina input.fastq”
5. Generate reference genome index file (Use GCF ones)

 Note : Files with “.fna” may need to be renamed to “.fa” Note : Better to use full directories rather than abbreviated ones Command : “STAR --runThreadN 48 -- runMode genomeGenerate --genomeDir genomeDirectory --genomeFastaFiles directory to fasta file(s) --sjdbGTFfile directory to genome GTF file" Note : Only needs to be done once

6. Use STAR to align genome to reference genome Command : “STAR --genomeDir directory to reference genome --readFilesIn directory to fastq files --quantMode GeneCounts”
In this case, directory to reference genome is : “/home/ehyang4/ncbi_dataset/data”

#### How to Use My Program

1. Unzip all files into .fastq format
2. Open terminal and cd to afolder with all the fastq files
3. Activate cutadapt env
4. Start Python in terminal via "python"
5. Copy and paste the #finished loop program from the jupyter notebook into terminal
6. Run it and it should begin processing all the fastq files in the folder and store them into an output folder automatically


---------------------




Index builds using bowtie-build: $ bowtie-build reference_sequence.fasta index_name

Generate clean .fa files using perl: perl -plane 's/\s+.+$//' < genome.fa > new_genome.fa

https://www.biostars.org/p/147955/ - replace unknown characters with N

Use VIM and this command: :%s/R/N/g 
