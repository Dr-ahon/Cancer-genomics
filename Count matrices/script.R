library(DESeq2)
library(pheatmap)

# Load count matrix
x = read.table("sample.counts", row.names=1, header=T, sep=",")
s = read.table("sample.info", header=T, row.names=1, colClasses=c("character", "factor"))

# Create DESeq2 object
dds = DESeqDataSetFromMatrix(countData = x, colData = s, design = ~ condition)

# Run a differential expression analysis (Tumour vs. Normal) using a log-fold change threshold of 1
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds <- DESeq(dds)
res <- results(dds, name="condition_Tumour_vs_Normal", lfcThreshold=1)

# Generate an MA-plot
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results
plotMA(res, ylim=c(-5,5))

# Plot the normalized counts for the GJB2 gene
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts
plotCounts(dds, gene="GJB2", intgroup="condition", normalized = TRUE)


# Generate a PCA plot of the samples using the transformed count data
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extracting-transformed-values
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))


# Visualize the differential gene expression results as a heatmap
# Take the top 20 genes according to the adjusted p-value
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
top_genes <- head(order(res$padj), 20)
pheatmap(assay(vsd)[top_genes,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Export the significant results (padj < 0.01) to a CSV file
# Tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exporting-results-to-csv-files
signifGenes <- subset(res, padj < 0.01)
write.csv(as.data.frame(signifGenes), file = "count_matrices_sign_result.csv")
