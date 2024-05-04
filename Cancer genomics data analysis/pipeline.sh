#!/bin/bash

# unzipping data files
gunzip hg19.fa.gz tu.r1.fq.gz tu.r2.fq.gz wt.r1.fq.gz wt.r2.fq.gz

# indexing reference genome
./bwa index -p  hg19 hg19.fa

# aligning reads to the reference

./bwa mem -t 12 reference/hg19 raw_data/tu.r1.fq raw_data/tu.r2.fq >results/bwa/tu_aligned.sam
./bwa mem -t 12 reference/hg19 raw_data/wt.r1.fq raw_data/wt.r2.fq >results/bwa/wt_aligned.sam

# conversion form SAM to BAM and sorting

samtools view -S -b tu_aligned.sam > tu_aligned.bam
samtools view -S -b wt_aligned.sam > wt_aligned.bam

# sorting of the BAM files

samtools sort -o tu_aligned.sorted.bam tu_aligned.bam
samtools sort -o wt_aligned.sorted.bam wt_aligned.bam

# indexing of the sorted BAM files

samtools index -o tu_aligned.sorted.bam
samtools index -o wt_aligned.sorted.bam

# cutting out just the chromosome of interest

samtools view -b tu_aligned.sorted.bam chrX:20000000-40000000 > tu_region_output.bam
samtools view -b wt_aligned.sorted.bam chrX:20000000-40000000 > wt_region_output.bam

# getting the read depth from the BAM files

samtools depth tu_region_output.bam > tu_depth.txt
samtools depth wt_region_output.bam > wt_depth.txt
