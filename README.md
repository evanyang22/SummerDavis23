# SummerDavis23

For questions: email evan.yang@emory.edu

Link to extensive documentation : https://tutormd.notion.site/Aligners-Documentation-e4fee84eba2448fa9df00fd82dd48e9a 

## miRNA alignment using miRDeep2
Functions of mirDeep2 (https://github.com/rajewsky-lab/mirdeep2)
1. Mapper.pl : processes reads and maps them to reference genome, is a pre-requisite for miRDeep2.pl and quantifier.pl
2. miRDeep2.pl :  detects novel microRNAs
3. Quantifier.pl : maps and quantifiers microRNAs based on microRNA precursors

#### Steps (6_14_Evan_miRNA.ipynb)
1. Trim files using trimmer (given in miRNA_trimgalore_script.py)
2. Store trimmed files in a folder named "Input"
3. Index builds using bowtie-build in the terminal: $ bowtie-build reference_sequence.fasta index_name

Should generate .EBWT files
 
5. In the Jupyter notebook code, change the precursors, species, matures, output_dir, genome, and outputDir to fit your needs

precursors/matures= fasta files with the precursor miRNA and mature miRNA sequences, respectively

species= 3 letter index of the species of interest, used to extract the precursor and mature miRNAs from the above fasta files

output_dir= where the results of mapper.pl are stored

genome= prefix of the .ebwt files

outputDir= output directory with the final results of each sample

5. Generate clean .fa files using perl: perl -plane 's/\s+.+$//' < genome.fa > new_genome.fa

This is primarily necessary for the miRDeep2.pl function

This is necessary for 'hsa_hairpin_clean.fa' and 'hsa_mature_clean.fa' to remove white spaces
This is also necessary for 'newHuman-Copy.fa' to remove white spaces in the reference genome

6. Open a terminal and "cd" to the correct directory, use "python" to open a python terminal, and copy and paste the first block of jupyter notebook into the terminal

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

#### How to Use My Program (6_15_Troubleshoot_tRNA.ipynb)

1. Unzip all files into .fastq format
2. Open terminal and cd to afolder with all the fastq files
3. Activate cutadapt env
4. Start Python in terminal via "python"
5. Copy and paste the first chunk of code from the jupyter notebook into terminal
6. Run it and it should begin processing all the fastq files in the folder and store them into an output folder automatically


## Acknowledgements

miRNA_trimgalore_script.py and TotalRNA_trimgalore_script.py created by Felipe

---------------------






Generate clean .fa files using perl: perl -plane 's/\s+.+$//' < genome.fa > new_genome.fa

https://www.biostars.org/p/147955/ - replace unknown characters with N

Use VIM and this command: :%s/R/N/g 
