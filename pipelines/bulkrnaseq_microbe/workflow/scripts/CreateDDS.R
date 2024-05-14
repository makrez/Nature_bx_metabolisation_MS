## The script expects these inputs:

## 1) path to a table of counts containing EnsemblID as first column and counts for each sample in the remaining columns
## 2) path to a table specifying the experimental group for each sample. 
## 3) A table with additional info on genes from biomaRt

## The following outputs are produced:
## DESEqDataSet with the results from DESeq() call
## DESeqTransform object with results from rlog()
## pdf with PCA plot of 500 most variable genes
## pdf with heatmap of 50 most highly expressed genes

args <- commandArgs(TRUE)

CreateDDS <- function(count_table, design, gene_table, delim, factorial.design, name.dds,
name.rld, name.heatmap, name.pca){
  library("DESeq2")
  library("RColorBrewer")
  library("gplots")
  library("ggplot2")
  
  if(factorial.design == "twofactorial"){
    unifactorial <- FALSE
  } else {
    unifactorial <- TRUE
  }
  
  data <- read.table(count_table, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  
  if(unifactorial){
    expDesign <- read.table(design, header=FALSE, sep="\t", colClasses = c("character", "character"))
    colnames(expDesign) <- c("Sample", "Group")
    expDesign$Group <- as.factor(expDesign$Group)
  } else {
    expDesign <- read.table(design, header=FALSE, sep="\t", colClasses = c("character", "character", "character"))
    colnames(expDesign) <- c("Sample", "Group", "Group2")
    expDesign$Group <- as.factor(expDesign$Group)
    expDesign$Group2 <- as.factor(expDesign$Group2)
  }
  
  # double-check that the sample ID ends with delim, of not add. 
  for (i in 1:length(expDesign$Sample)){
    if(!endsWith(expDesign$Sample[i], delim)) { expDesign$Sample[i] <- paste(expDesign$Sample[i], delim, sep="") }
  }

  # Make sure the order of samples is the same in data and expDesign
  # Note: If the sample names do not match exactly between expDesign and data, the reordered expDesign will 
  # contain NA in this row, and the DESeqDataSetFromMatrix step will fail
  expDesign <- expDesign[match(colnames(data), as.character(expDesign$Sample)),]
  
  if(unifactorial){
    diff <- DESeqDataSetFromMatrix(data, expDesign, ~Group)
  } else {
    diff <- DESeqDataSetFromMatrix(data, expDesign, ~Group2 + Group)
  }
  
#  ## Add additional info on genes as extra columns in the DESeqDa --------
  # gene_id (=ensemblID)
  # gene_symbol(s). Depending on the species, there can be multiple, e.g. that for the species
  # of interest + human. All are retrieved. 
  # entrezgene     
  # chromosome_name
  # heatmap.id (=gene_symbol if exists, gene_id otherwise)
  
  annotation <- read.table(gene_table, header = TRUE, sep = "\t", quote = "")
  # reorder so that the rows are ligned up correctly
  annotation <- annotation[match(rownames(diff), annotation$gene_id),]
  # Add gene info to the DESeqDataSet
  mcols(diff) <- cbind(mcols(diff), annotation)

  
  # DE analysis -------------------------------------------------------------
  
  dds <- DESeq(diff, betaPrior = TRUE)
  saveRDS(dds, name.dds)
  
  rld <- rlog(dds, blind=TRUE) # apply a regularized log transformation, ignoring information about experimental groups
  saveRDS(rld, name.rld)  
  
# PCA -------------------------------------------------------------------
  if(unifactorial){
    p <- plotPCA(rld, intgroup=c("Group"))
    pdf(name.pca)
    print(p, ntop=500)
    dev.off()
  } else {
    rv <- rowVars(assay(rld))
    ntop<-500   # number of genes to use
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    results <- prcomp(t(assay(rld)[select, ]))
    results<-as.data.frame(results$x)
    
    p <- ggplot(results, aes(x=PC1, y=PC2)) +
      geom_point(aes(colour = rld$Group, shape = rld$Group2), size = 2) +
      scale_colour_discrete(name='Group') +
      scale_shape_discrete(name='Confounder')
    pdf(name.pca)
    print(p)
    dev.off()
  }
  

# Heatmap of most highly expressed genes ----------------------------------
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  pdf(name.heatmap)
  heatmap.2(assay(rld)[select,], col = hmcol, trace="none", margin=c(10, 6),
            labCol=colnames(dds), cexRow = 0.4)
  dev.off()

}


CreateDDS(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])




