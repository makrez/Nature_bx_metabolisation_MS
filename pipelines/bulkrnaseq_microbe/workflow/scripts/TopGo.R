args <- commandArgs(TRUE)


###########################################################################################
# RunTopGO()
###########################################################################################
# Function to run TopGo and produce two output files per analysis:
# *GenesInGoTerms.rds: Contains all gene ids for the genes in the top nbr.nodes GO terms
# *topGoResults.txt: Contains the TopGo results for the top nbr.nodes GO terms from each ontology

# Runs TopGo using Weight01-Fisher, Weight01-KS and Classic-Fisher. Genes will be ranked based on Weight01.Fisher.If neither of the 3 P-values is < 0.05, the term will not be output

# The GO graph will be visualised only for a maximum of 10 nodes (or nbr.nodes if nbr.nodes<10). The nodes will be
# ranked based on Weight01.Fisher but no P-value cutoff is applied

# Arguments:
# DESeq2Table: Output of Deseq2_Contrasts_Rlog(). Must contain the columns pvalue, padj, gene_id (with ensembl ids). Expected file name is *.DEResults.original.rlog.txt
# fdr.threshold: P-value threshold below which a gene will be considered as significantly DE
# database: contains the mapping between gene ids and go terms. Must be a Bioconductor annotation package of the format org.XX.XX.db
# nbr.nodes: Number of top GO terms to be output. This is perhaps better than using a P-value cut-off because the ranks may be more meaningful than the P-values (see TopGo manual)

write_empty_plot <- function(name, DESeq2Table){
  output<-sub(".DEResults.original.rlog.txt", "", DESeq2Table)
  output_name<-paste(output, name, 'pdf', sep=".")
  pdf(output_name)
  plot.new()
  dev.off()
}

RunTopGO<-function(DESeq2Table, fdr.threshold, database, nbr.nodes, table_out, genes_out){

  fdr.threshold <- as.numeric(fdr.threshold)
  nbr.nodes <- as.integer(nbr.nodes)

  library(topGO)
  library(plyr)
  library(dplyr)
  
  DEResults<-read.table(DESeq2Table, header=TRUE, sep='\t', quote="")
  
  # Remove genes without a p-value or adjusted p-value
  DEResults<-DEResults[!is.na(DEResults$pvalue) & !is.na(DEResults$padj),]
  # Remove genes with baseMean1==1. These have 0 reads in all samples in this comparison
  # They can have a P-value of 1 if the gene has reads in some other samples (not in that comparison)
  # and sometimes they can have an adjusted P-value (this seems to happen e.g. when there is a big gap
  # in P-values between some top genes and the rest. It then doesn't matter how many genes are included in FDR correction
  DEResults<-DEResults[DEResults$baseMean1>1,]
  
  # create a named vector containing the adjusted P-value as entries and the ensembl ids as names
  geneList<-setNames(DEResults$padj, DEResults$gene_id)


  if(sum(geneList < fdr.threshold) == 0){
    warning('no differentially expressed genes')
    write('no differentially expressed genes', file = table_out)
    saveRDS(NULL, file=genes_out)
    write_empty_plot('empty', DESeq2Table)
    return('no differentially expressed genes')
  }
  
  # function to indicate DE genes
  topDiffGenes<-function(allScore, threshold=fdr.threshold){
    return(allScore<threshold)
  }
  
  onts<-c("MF", "BP", "CC")
  tab<-as.list(onts)
  names(tab)<-onts
  
  ann.genes<-list(MF=list(), BP=list(), CC=list())
  names(ann.genes)<-onts
  
  for(i in 1:3){
    write(onts[i], stdout())
    tgd<-tryCatch(new("topGOdata", 
             ontology=onts[i],
             allGenes=geneList, 
             geneSel=topDiffGenes,
             nodeSize=5,
             annot=annFUN.org,
             mapping=database,
             ID="ensembl"), error = function(e) onts[i])
    
    
    resultTopGo.weight01.Fisher<-tryCatch(runTest(tgd, algorithm="weight01", statistic="Fisher"), error =  function(e) onts[i])
    resultTopGo.weight01.KS<-tryCatch(runTest(tgd, algorithm="weight01", statistic="ks"), error =  function(e) onts[i])
    resultTopGo.classic.Fisher<-tryCatch(runTest(tgd, algorithm="classic", statistic="Fisher"), error =  function(e) onts[i])

    tab[[i]]<-tryCatch(GenTable(tgd, 
                       weight01.Fisher=resultTopGo.weight01.Fisher, 
                       weight01.KS=resultTopGo.weight01.KS,
                       classic.Fisher=resultTopGo.classic.Fisher, 
                       orderBy="weight01.Fisher", 
                       ranksOf="classic.Fisher", 
                       topNodes=nbr.nodes), error = function(e) onts[i])

    if(is.null(tab[[i]]) | length(tab[[i]]) == 1){
      write(paste('that did not work for', onts[i]), stdout())
      tab[[i]] <- onts[i]
      write_empty_plot(onts[i], DESeq2Table)
      # output<-sub(".DEResults.original.rlog.txt", "", DESeq2Table)
      # output_name<-paste(output, onts[i], 'pdf', sep=".")
      # pdf(output_name)
      # # plot(1:10, type = 'n')
      # plot.new()
      # dev.off()
      next
    }    
    # write(str(tab), stdout())
    
    # write(nrow(tab[[i]]), stdout())               
    tab[[i]]$ontology<-onts[i]

    tab[[i]]$weight01.Fisher <- as.numeric(tab[[i]]$weight01.Fisher)
    tab[[i]]$weight01.KS <- as.numeric(tab[[i]]$weight01.KS)
    tab[[i]]$classic.Fisher <- as.numeric(tab[[i]]$classic.Fisher)
	
	# Remove a GO term if none of the 3 P-values is < 0.05
    tab[[i]] <- dplyr::filter(tab[[i]], weight01.Fisher < 0.05 | weight01.KS < 0.05 | classic.Fisher < 0.05)

    if(nrow(tab[[i]]) == 0){
            write(paste('that did not work for', onts[i]), stdout())
      tab[[i]] <- onts[i]
      write_empty_plot(onts[i], DESeq2Table)
      # output<-sub(".DEResults.original.rlog.txt", "", DESeq2Table)
      # output_name<-paste(output, onts[i], 'pdf', sep=".")
      # pdf(output_name)
      # # plot(1:10, type = 'n')
      # plot.new()
      # dev.off()
      next
    }
	

	
	# Create an R data object that contains the gene IDs for all genes in the GO terms
	for (j in 1:nbr.nodes){
      ann.genes[[i]][[j]]<-unname(unlist(genesInTerm(tgd, tab[[i]]$GO.ID[j])))
    }

    names(ann.genes[[i]])<-tab[[i]]$GO.ID
	
    # visualise the results:  
    
	output<-sub(".DEResults.original.rlog.txt", "", DESeq2Table)
    output_name<-paste(output, onts[i], sep=".")
    
	nbr.nodes <- as.integer(nbr.nodes)
    nodesToVisualise<-min(10, nbr.nodes)
    printGraph(tgd, resultTopGo.weight01.Fisher, firstSigNodes=nodesToVisualise, fn.prefix=output_name, useInfo="all", pdfSW=TRUE)
    
    
  }

  saveRDS(ann.genes, file=genes_out)
  
  topGOResults <- tryCatch(rbind.fill(tab), error = function(e) NULL)
  write.table(topGOResults, table_out, sep='\t', quote=FALSE, row.names=FALSE)
  
}
###########################################################################################


RunTopGO(args[1], args[2], args[3], args[4], args[5], args[6])