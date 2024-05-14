# Sample ids CANNOT start with number

args<-commandArgs(TRUE)


###########################################################################################
# Deseq2_Contrasts_Rlog()
###########################################################################################
# function that performs pairwise comparison through contrasts and writes the results to file
# DESeqDataSet is an object produced by DESeq() of DESeq2
# factorName, level1 and level2 specify the contrast to be run
# TableGeneNames is a table that has at least one column called "gene_id" which contains the same gene ids as the counts table, typically ensembl ids, and a column called "gene_symbol" (plus optional columns with other info for each gene)
# rlog_in is an object produced by rlog() of DESeq2

# Also, creates a table that specifies the experimental group for each sample needed by shiny app named expGroup_filename
# and a separate table with all rlog values which will be used by createTopTables() and will be named rlog_filename

###########################################################################################

Deseq2_Contrasts_Rlog<-function(DESeqDataSet, factorName, level1, level2, TableGeneNames, rlog_in, rlog_filename, expGroup_filename, cutoff){
  library(DESeq2)
  temp<-results(DESeqDataSet, contrast=c(factorName, level1, level2), alpha = cutoff)
  factorVector<-DESeqDataSet@colData@listData[[factorName]]
  f1_columns <- which(factorVector==level1)
  f2_columns <- which(factorVector==level2)

  countsT<-cbind(apply(counts(DESeqDataSet, normalized=TRUE)[, f1_columns, drop = FALSE], 1, mean),
                 apply(counts(DESeqDataSet, normalized=TRUE)[, f2_columns, drop = FALSE], 1, mean))
  countsT<-as.data.frame(countsT)
  countsT$V1<-as.numeric(format(countsT$V1, digits=3, scientific=FALSE))
  countsT$V2<-as.numeric(format(countsT$V2, digits=3, scientific=FALSE))
  countsT$V3<-as.numeric(format((apply(countsT, 1, mean)+1), digits=3, scientific=FALSE))
  countsT<-cbind(countsT, counts(DESeqDataSet)[, f1_columns],
                 counts(DESeqDataSet)[, f2_columns],
                 counts(DESeqDataSet, normalized=TRUE)[, f1_columns],
                 counts(DESeqDataSet, normalized=TRUE)[, f2_columns])
  # add rlog transformed counts for the groups of interest
  countsT<-cbind(countsT, assay(rlog_in)[, f1_columns], assay(rlog_in)[, f2_columns])

  colnames(countsT)<-c(paste('normMean', level1, sep='.'), paste('normMean', level2, sep='.'),
                       'baseMean1',
                       colnames(DESeqDataSet)[f1_columns],
                       colnames(DESeqDataSet)[f2_columns],
                       paste(colnames(DESeqDataSet)[f1_columns], 'norm', sep='.'),
                       paste(colnames(DESeqDataSet)[f2_columns], 'norm', sep='.'),
                       paste(colnames(assay(rlog_in))[f1_columns], 'rlog', sep='.'),
                       paste(colnames(assay(rlog_in))[f2_columns], 'rlog', sep='.'))
  resultsTable<-merge(countsT, as.data.frame(temp), by="row.names", all.x=TRUE, all.y=TRUE)
  colnames(resultsTable)[1]<-'gene_id'
  resultsTable$gene_id<-as.character(resultsTable$gene_id) # strange formatting because column was produced from rownames, it seems..


  temp<-merge(TableGeneNames, resultsTable, by.x='gene_id', by.y='gene_id', all.x=FALSE, all.y=TRUE)


  write.table(temp[order(as.character(temp$gene_id)),], paste("results/DESeq", "Contrasts", paste(level1, level2, "DEResults", "original", "rlog", "txt",sep="."), sep="/"), quote=F, row.names=F, sep="\t")
#  write.table(temp[order(temp$padj),], paste(level1, level2, 'DEResults', "rlog", "txt",sep="."), quote=F, row.names=F, sep="\t")

	if(!file.exists(expGroup_filename)){
	expGroups<-data.frame(rownames(colData(DESeqDataSet)), colData(DESeqDataSet)[factorName])
    colnames(expGroups)<-c("id", "group")
	write.table(expGroups, expGroup_filename, sep="\t", quote=FALSE, row.names=FALSE)
	}

# Write a table containing rlog values for all samples, followed by a column with heatmap.id
	if(!file.exists(rlog_filename)){
		temp_rlog <- TableGeneNames[, c("gene_id", "heatmap.id")]
		temp_rlog <- merge(assay(rlog_in), temp_rlog, by.x = "row.names", by.y = "gene_id", all.x = TRUE, all.y = FALSE)
		rownames(temp_rlog) <- temp_rlog$Row.names
		temp_rlog <- temp_rlog[, -1]
		write.table(temp_rlog, rlog_filename, sep="\t", quote=FALSE, row.names=TRUE)
	}

}


################################################################################################################

library(DESeq2)
dds <- readRDS(args[1])
rld <- readRDS(args[5])


# Table with gene info
GeneTable <- rowData(dds)[, 1:5]
colnames(GeneTable)[2] <- "gene_symbol"

Deseq2_Contrasts_Rlog(dds, args[2], args[3], args[4], GeneTable, rld, args[6], args[7], args[8])
