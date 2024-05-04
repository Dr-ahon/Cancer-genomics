

# Cancer genomics data 

### Data preparation

downloading datasets

```bash
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
wget "https://gear-genomics.embl.de/data/.exercise/tu.r1.fq.gz"
wget "https://gear-genomics.embl.de/data/.exercise/tu.r2.fq.gz"
wget "https://gear-genomics.embl.de/data/.exercise/wt.r1.fq.gz"
wget "https://gear-genomics.embl.de/data/.exercise/wt.r2.fq.gz"
```

unzipping data files
``
```bash
gunzip hg19.fa.gz tu.r1.fq.gz tu.r2.fq.gz wt.r1.fq.gz wt.r2.fq.gz
```

### Analysis

indexing reference genome

```bash
 ./bwa index -p  hg19 hg19.fa
```

aligning reads to the reference

```bash
./bwa mem -t 12 reference/hg19 raw_data/tu.r1.fq raw_data/tu.r2.fq >results/bwa/tu_aligned.sam
./bwa mem -t 12 reference/hg19 raw_data/wt.r1.fq raw_data/wt.r2.fq >results/bwa/wt_aligned.sam
```

conversion form SAM to BAM and sorting

```bash
samtools view -S -b tu_aligned.sam > tu_aligned.bam
samtools view -S -b wt_aligned.sam > wt_aligned.bam
```

sorting of the BAM files

```bash
samtools sort -o tu_aligned.sorted.bam tu_aligned.bam
samtools sort -o wt_aligned.sorted.bam wt_aligned.bam
```

indexing of the sorted BAM files

```bash
samtools index -o tu_aligned.sorted.bam
samtools index -o wt_aligned.sorted.bam
```

cutting out just the chromosome of interest

```bash
samtools view -b tu_aligned.sorted.bam chrX:20000000-40000000 > tu_region_output.bam
samtools view -b wt_aligned.sorted.bam chrX:20000000-40000000 > wt_region_output.bam
```

getting the read depths from the BAM files

```bash
samtools depth tu_region_output.bam > tu_depth.txt
samtools depth wt_region_output.bam > wt_depth.txt
```

###  Plotting

running the script.py script, that divides the reads into windows, counts them and plots them
