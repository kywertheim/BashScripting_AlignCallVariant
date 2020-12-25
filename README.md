# BashScripting_AlignCallVariant
Context: This is the solution to the third exam of the Coursera course entitled 'Command Line Tools for Genomic Data Science'. The sequencing reads produced in a WGS study of an Arabidopsis thaliana strain are in the file 'wu_0_A_wgs.fastq'. The Bash script catalogues their genetic variations from the reference genome in the file 'wu_0.v7.fas'.

About:
1. Working within a Unix shell, the Bash script 'analysis.sh' indexes the reference genome, aligns the reads to the indexed genome, compiles a list of candidate sites of variation, and calls the variants.
2. Some of the commands are redundant. For example, 'bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_align.sam' and 'bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.full.sam' generate equivalent results.

Files:
1. 'gencommand_proj3_data.tar.gz' is the archive file. Due to its size, it is not available in the repository. However, the Bash script includes a command that downloads it.
2. 'analysis.sh' should be placed in the same directory as 'gencommand_proj3_data.tar.gz'.

Software:
1. samtools v.1.2.
2. bowtie v.2.2.5.
3. bcftools v.1.2.
