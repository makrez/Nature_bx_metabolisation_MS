###########################################################################################
# Extension for IBU RNAseq pipeline for microbial datasets
# Simone Oberhaensli IBU, 2021-07-22
#
# This script performs the following tasks:
# reads in experimental design, count table and feature table
# perform DESeq2 analysis
# generate output tables and QC plots
#
# on IBU cluster run like this:
# Rscript --vanilla RNAseq_microbial_DESeq2.R expdesign.txt finalCounts.txt features.csv
##########################################################################################


## Preamble --------------------------------------------------------------------
# test if input files are there
args = commandArgs(trailingOnly=TRUE)
if (length(args)<=2) {
  stop("please indicate input files: experimental_design.txt, *finalCounts.txt, features.csv", 
       call.=FALSE)
} 


# prep: load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


# function: read experimental design and create coldata
read_exp <- function(exp_design){
  t.exp <- read.table(exp_design, header=FALSE)
  colnames(t.exp) <- c("sample", "condition")
  return(t.exp)
}

# function: based on experimental desing generate coldata
create_coldata <- function(counts, exp){
  samples <- colnames(counts)
  condition <- rep("dummy_content", length(samples))
  cdat <- data.frame(condition, row.names=samples)
  cdat$condition <- exp$condition[match(rownames(cdat), exp$sample)]
  return(cdat)
}

# function: generate table with normalized counts and features
crt_normcounts <- function(features){
  counts_norm <- counts(dds, normalized=T)
  d.counts_norm <- as.data.frame(counts_norm)
  d.counts_norm$gene <- rownames(d.counts_norm)
  # add gene symbol and function to it
  d.counts_norm_ant <-  merge(d.counts_norm, features, by = "gene", 
                          all.x = TRUE, all.y = TRUE) 
  write.csv(d.counts_norm_ant, "normalized_counts_fun.csv", row.names = FALSE)
}


# generate contrast output files with annotations info
crt_DEG_files <- function(contrast, features){
  # simples analysis, no LFCshrinkeage
  res_ctr <- results(dds, name=contrast, alpha=0.05) 
  summary(res_ctr)
  d.res_ctr  <- data.frame(res_ctr)
  d.res_ctr$gene <- rownames(d.res_ctr)
  # add gene symbol and function to it
  d.res_ctr_ant <-  merge(d.res_ctr, features, by = "gene", 
                          all.x = TRUE, all.y = TRUE) 
  d.res_ctr_ant.srt <- d.res_ctr_ant[order(d.res_ctr_ant$padj, decreasing = T),]
  write.csv(d.res_ctr_ant.srt, paste0(contrast, "_fun.csv"), row.names = FALSE)
}


## Analysis ---------------------------------------------------------------------

sink('analysis-output.txt')

## data preparation
f.expdesign <- args[1] # exp. design
f.counts <- args[2] # counts table
f.features <- args[3]

# nice entry in log
cat("\nexperimental design: ")
cat(f.expdesign)
cat("\ncount table: ")
cat(f.counts)
cat("\ntable with features or annotations: ")
cat(f.features)

# read features file and the experimental design 
d.feat <- read.csv(f.features, header = F)
colnames(d.feat) <- c("gene","symbol","product")
d.exp <- read_exp(f.expdesign)

# read counts table and create coldata 
d.counts <- read.table(f.counts, header =T, row.names = 'Geneid')
# remove trailing underscore from sample names
colnames(d.counts) <- sub("_$","",colnames(d.counts))
cat("\nColdata looks like this, please check:\n")
(coldata <- create_coldata(d.counts, d.exp))

# rownames coldata must have the same order als colnames of count table
cat("\ncolnames count table correspond to rownames in coldata: ")
cat(all(rownames(coldata) == colnames(d.counts)))

## DESeq object generation
# very simple one factor design
dds <- DESeqDataSetFromMatrix(countData = d.counts,
                              colData = coldata,
                              design = ~ condition)


dds <- DESeq(dds)
saveRDS(dds, "dds.Rds")

crt_normcounts(d.feat)

cat("\nAvailable contrasts:\n")
resultsNames(dds)

# generate result files for all contrasts
for (i in resultsNames(dds)[-1]){
  cat("\nRunning contrast: ")
  print(i)
  crt_DEG_files(i,d.feat)
}

sink()

## Additional QC plots ----------------------------------------------------------

# raw count transformations
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
# this gives log2(n + 1)
#ntd <- normTransform(dds)

# pca plot
pdf("PCA.pdf")
plotPCA(vsd, intgroup="condition")
dev.off()

# Heatmap of the sample-to-sample distances
pdf("heatmap_sampledist.pdf")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


## packages and versions used
sink('analysis-output.txt', append=TRUE)
cat("\nSession info:\n")
sessionInfo()
sink()