#Download the archive.
wget https://d396qusza40orc.cloudfront.net/gencommand/gencommand_proj3_data.tar.gz

#Extract the files from the archive.
tar -xf gencommand_proj3_data.tar.gz

#1. How many sequences are there in the genome file?
grep ">" wu_0.v7.fas | wc -l
grep -c "^>" wu_0.v7.fas

#2. What is the name of the first sequence in the genome file?
grep ">" wu_0.v7.fas
grep "^>" wu_0.v7.fas | head -3 | tail -1

#3. What is the name of the last sequence in the genome file?
grep ">" wu_0.v7.fas
grep "^>" wu_0.v7.fas | tail -1

#Learn more about bowtie2-build.
bowtie2-build &> bowtie2-build.log

#Generate the bowtie2 index for the genome file.
bowtie2-build wu_0.v7.fas ./wu_0

#4. How many index files does this operation create?
ls | grep "bt2" | wc -l

#5. What is the 3-character extension for the index files?
ls

#6. How many reads are in the FASTQ file?
wc -l wu_0_A_wgs.fastq

#Run bowtie2 to produce end-to-end read alignments.
bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_align.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.full.sam

#Run bowtie2 to produce potential partial read alignments.
bowtie2 -p 4 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_lalign.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.local.sam --local

#7. How many end-to-end read alignments does this operation create?
bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_align.sam
cat out.full.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l

#8. How many partial read alignments does this operation create?
bowtie2 -p 4 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_lalign.sam
cat out.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l

#9. How many reads does the operation map with the full-match setting?
bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_align.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.full.sam

#10. How many reads does the operation map with the local setting?
bowtie2 -p 4 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_lalign.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.local.sam --local

#11. How many reads does the operation map to multiple alignments with the full-match setting?
bowtie2 -p 4 -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_align.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.full.sam

#12. How many reads does the operation map to multiple alignments with the local setting?
bowtie2 -p 4 --local -x wu_0 wu_0_A_wgs.fastq -S wu_0_A_lalign.sam
bowtie2 -x wu_0 -U wu_0_A_wgs.fastq -S out.local.sam --local

#13. How many alignments with insertions and/or deletions does the operation produce with the full-match setting?
less wu_0_A_align.sam | cut -f6 | grep "I" | grep "D" | wc -l
less wu_0_A_align.sam | cut -f6 | grep "I" | wc -l
less wu_0_A_align.sam | cut -f6 | grep "D" | wc -l
cat out.full.sam | cut -f6 | grep -c "[I,D]"

#14. How many alignments with insertions and/or deletions does the operation produce with the local setting?
less wu_0_A_lalign.sam | cut -f6 | grep "I" | grep "D" | wc -l
less wu_0_A_lalign.sam | cut -f6 | grep "I" | wc -l
less wu_0_A_lalign.sam | cut -f6 | grep "D" | wc -l
cat out.local.sam | cut -f6 | grep -c "[I,D]"

#Convert the sam file to a bam file, sort it, and then index it.
samtools view -b -T wu_0.v7.fas wu_0_A_align.sam > wu_0_A_align.bam
samtools sort wu_0_A_align.bam wu_0_A_align.sorted
samtools index wu_0_A_align.sorted.bam

#Identify the candidate sites of variation.
samtools mpileup -f wu_0.v7.fas wu_0_A_align.sorted.bam > wu_0_A.mpileup
samtools mpileup -uv -f wu_0.v7.fas wu_0_A_align.sorted.bam > wu_0_A.vcf

#15. How many entries does the operation produce?
cat wu_0_A.vcf | cut -f1 | grep "Chr3" | wc -l
cat wu_0_A.vcf | grep -v "^#" | cut -f1 | grep -c "^Chr3"

#16. How many entries with 'A' as the corresponding genome letter does the operation produce?
cat wu_0_A.vcf | cut -f4 | awk '$1 == "A"' | wc -l
cat wu_0_A.vcf | grep -v "^#" | cut -f4 | grep -P "^A$" | wc -l

#17. How many entries with exactly 20 supporting reads (read depth) does the operation produce?
cat wu_0_A.vcf | cut -f8 | grep "DP=20" | wc -l
cat wu_0_A.vcf | grep -v "^#" | grep -c "DP=20"

#18. How many entries representing indels does the operation produce?
cat wu_0_A.vcf | cut -f8 | grep "INDEL" | wc -l
cat wu_0_A.vcf | grep -v "^#" | grep -c "INDEL"

#19. How many entries does the operation produce at position 175672 on Chr1?
cat wu_0_A.vcf | cut -f1,2,3,4,5,6,7,8 | grep "Chr1" | awk '$2 == "175672"' | wc -l
cat wu_0_A.vcf | grep -v "^#" | cut -f1,2 | grep "Chr1" | cut -f2 | grep -P "^175672$" | wc -l

#Identify the candidate sites of variation and record them in a bcf file.
samtools mpileup -g -f wu_0.v7.fas wu_0_A_align.sorted.bam > wu_0_A.bcf

#Call the variants with the multi-allelic caller, presenting only the variant sites in an uncompressed VCF file.
bcftools call -m -v -O v wu_0_A.bcf > wu_0_A.final.vcf

#20. How many variants does this operation call on Chr3?
cat wu_0_A.final.vcf | cut -f1 | grep "Chr3" | wc -l
cat wu_0_A.final.vcf | grep -v "^#" | cut -f1 | sort | uniq -c | grep "Chr3"

#21. How many variants representing A-to-T SNPs does this operation call?
cat wu_0_A.final.vcf | cut -f4,5 | awk '$1=="A" && $2=="T"' | wc -l
cat wu_0_A.final.vcf | grep -v "^#" | cut -f4,5 | grep -P "^A\tT$" | wc -l

#22. How many variants representing indels does this operation call?
cat wu_0_A.final.vcf | cut -f8 | grep "INDEL" | wc -l
cat wu_0_A.final.vcf | grep -v "^#" | grep -c "INDEL"

#23. How many variants with exactly 20 supporting reads (read depth) does the operation call?
cat wu_0_A.final.vcf | cut -f8 | grep "DP=20" | wc -l
cat wu_0_A.final.vcf | grep -v "^#" | grep -c "DP=20"

#24. What kind of variant does the operation call at position 11937923 on Chr3?
cat wu_0_A.final.vcf | grep "Chr3" | awk '$2 == "11937923"'
cat wu_0_A.final.vcf | grep -v "^#" | cut -f1-5 | grep "Chr3" | grep "11937923"