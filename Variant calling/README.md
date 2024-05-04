# Variant calling

### Data preparation

unzip data files

```bash
gunzip chr7.fa.gz
gunzip R1.fastq.gz
gunzip R2.fastq.gz
```

install bwa

```bash
git clone https://github.com/lh3/bwa.git
cd bwa
make
```

first we need to index our reference data 
- -p to set the prefix of output files 

```bash
 ./bwa index -p chr7 chr7.fa
```

we have reads 150bp long, we will therefore use BWA-MEM to align reads 
- -t for set number of threads it will be running on

```bash
./bwa mem -t 12 reference_data/chr7 raw_data/R1.fastq raw_data/R2.fastq >results/bwa/R.sam
```

We convert the resulting SAM file to BAM format, to reduce it in size and enable indexing and sort it.

```bash
samtools view -S -b R.sam > R.bam
samtools sort -o R.sorted.bam R.bam
```

### Variant calling

Now we do the variant calling
-  -Ou set output type to uncompressed BCF 
- ls
- ls-v output variants sites only
-  -Ob set output type to compressed BCF
- -m choose multiallelic caller (should be better than -c consensus caller)

```bash
bcftools mpileup -Ou -f ../../reference_data/chr7.fa R.sorted.bam | bcftools call -mv -Ov -o calls.vcf
```

Now run the result through Variant Effect Predictor (https://www.ensembl.org/Tools/VEP)

### Dscussion

For the causative variant I'd nominate the variant on the location 7:2915243-2915243. It causes an appearance of a new stop codon, the most severe consequence of them all, which can easily cause some pathogenic behavior.
