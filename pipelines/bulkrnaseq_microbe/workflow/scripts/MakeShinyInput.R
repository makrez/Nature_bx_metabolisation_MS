args <- commandArgs(TRUE)

###########################################################################################
# MakeShinyInput()
###########################################################################################

# workingDir must end with "/"

MakeShinyInput <- function(rlog, topgo_in, topgo_out, de_tables, genes_in_terms, output){
  if(topgo_out != "NULL"){
    TopGoTables(topgo_in, topgo_out)
  }
  AddGoToTopTables(de_tables, rlog, genes_in_terms, output)
}

###########################################################################################
# TopGoTables()
###########################################################################################
# takes as input *topGoResults.txt files produced by RunTopGO()

TopGoTables<-function(results_files, output){
  # create *topGoResults.rds
  topGoResults<-lapply(results_files, read.table, header=TRUE, sep="\t", quote="")
  names(topGoResults)<-c(sub("\\.topGoResults\\.txt", "", basename(results_files)))
  saveRDS(topGoResults, file = output)
}
  


#######################################################################################
# readTopTables()
#######################################################################################
# In results tables from DESeq2, keep only genes that have a gene_id and a padj value
# If a gene_id occurs multiple times, keep only the first occurrence
# Finally, order by increasing pvalue

readTopTable <- function(inputFile) {
  topTable<-read.table(inputFile, sep="\t", header=TRUE, 
                       stringsAsFactors=FALSE, colClasses=c(entrezgene="character"), quote="")
  topTable<-topTable[!is.na(topTable$gene_id) & !is.na(topTable$padj),]
  topTable<-topTable[!duplicated(topTable$gene_id),]
  # when a gene_id occurs multiple times, only the first occurrence is kept
  topTable<-topTable[order(topTable$pvalue),]
  topTable
}
#######################################################################################


#######################################################################################
# createTopTables()
#######################################################################################
# Creates a list of dataframes where the first N elements are the DE results tables
# for each of N pairwise comparisons
# Element N+1 contains the rlog values for all samples

# rlog_path is path to a table containing the rlog values for all samples (produced by Deseq2_Contrasts_Rlog())

createTopTables<-function(DESeqTables, rlog){
  topTables<-lapply(DESeqTables, readTopTable)
  names(topTables) <- c(sub("\\.DEResults\\.original\\.rlog\\.txt", "", basename(DESeqTables)))
  
  # as the last element in the list, add a dataframe with the rlog values for all samples
  rlog<-read.table(rlog, header=TRUE, sep='\t', quote="")
  topTables$rlog<-rlog
  
  invisible(topTables)

} 


###########################################################################################
# AddGoToTopTables()
###########################################################################################
# Produces topTables using createTopTables()
# To each table of DE results, it will add nbr.nodes*3 (for each ontology) columns of 0/1's which indicate the genes in a particular GO term

# Note: The function assumes that the order of tables is the same in topTables and genesInGoTerms. This is the case because list of input files is expanded from the same contrast vector

# inputDirectory is a directory containing one or several *DEResults.original.rlog.txt files
# rlog_path is path to a table containing the rlog values for all samples (produced by Deseq2_Contrasts_Rlog())
# TableGeneNames is a table containing gene_id and gene_symbol as first two columns



AddGoToTopTables<-function(de_tables, rlog, genes_in_terms, output){

  topTables<-createTopTables(de_tables, rlog)

  if(!is.na(genes_in_terms) && genes_in_terms != "NULL"){
    genesInGoTerms<-lapply(genes_in_terms, readRDS)
    names(genesInGoTerms)<-c(sub("\\.GenesInGoTerms\\.rds", "", basename(genes_in_terms)))
  
    onts<-c("MF", "BP", "CC")
    
    
    for (i in 1:length(genesInGoTerms)) {
      # number of comparisons
      for (j in 1:3) {
        # number of ontologies
        outputMatrix <-
          matrix(
            data = 0,
            ncol = length(genesInGoTerms[[i]][[j]]),
            nrow = nrow(topTables[[i]])
          )
        if (ncol(outputMatrix) > 0) { # account for possibility that some onts may have no significant terms
          for (k in 1:length(genesInGoTerms[[i]][[j]])) {
            # number of genes per ontology
            outputMatrix[, k][is.element(topTables[[i]]$gene_id, genesInGoTerms[[i]][[j]][[k]])] <-
              1
          }
          dimnames(outputMatrix) <-
            list(topTables[[i]]$gene_id, paste(onts[j], names(genesInGoTerms[[i]][[j]]), sep =
                                                 "_"))
          outputMatrix <- as.data.frame(outputMatrix)
          topTables[[i]] <-
            merge(
              topTables[[i]],
              outputMatrix,
              by.x = "gene_id",
              by.y = "row.names",
              all.x = TRUE,
              all.y = TRUE
            )
        } else {
          next
        }
      }
    }
  }
    
  saveRDS(topTables, file=output)
}
###########################################################################################

topgo_in <- c(strsplit(args[2], ",")[[1]])
de_tables <- c(strsplit(args[4], ",")[[1]])
genes_in_goterms <- c(strsplit(args[5], ",")[[1]])

MakeShinyInput(args[1], topgo_in, args[3], de_tables, genes_in_goterms, args[6])



