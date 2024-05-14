args <- commandArgs(TRUE)


###########################################################################################
# runGSEA()
###########################################################################################
# Function to run GSEA and produce two output files per analysis:
# *gseKEGG.txt: Contains the GSEA results for KEGG pathways the top nbr.nodes GO terms from each ontology
# *gseKEGG.pdf: Graphical representation of top nbr.nodes significant KEGG pathways

# The GO graph will be visualised only for a maximum of 20 nodes (or nbr.nodes if nbr.nodes<20). The nodes will be
# ranked based on padj 

# Arguments:
# DESeq2Table: Output of Deseq2_Contrasts_Rlog(). Must contain the columns stat, entrezgene, gene_id (with ensembl ids), padj and log2FoldChange. 
#              Expected file name is *.DEResults.original.rlog.txt
# dbTranslationTable: Table with DB names: 1. column Bioconductor annotation package name (org.XX.XX.db); 2. column kegg organism code (e.g. mmu), 
#                                          3. column msigdb species name
# fdr.threshold: adjusted p-value threshold below which a KEGG pathway/GO term will be considered as significant
# database: contains the mapping between gene ids and go terms. Must be a Bioconductor annotation package of the format org.XX.XX.db
# kegg_organism: contains the kegg organism code (e.g. mmu, hsa)
# nbr.nodes: Number of top KEGG pathways to plot.
# gseKEGG_table: GSEA result table for KEGG pathways
# gseKEGG_plot:Graphical representation of top nbr.nodes significant KEGG pathwayss

# Dependencies: clusterProfiler, org.XX.CC.db, topGO, RColorBrewer, pathview, msigdbr



# make GO graph plot
# makeGOPlots <- function(DEResults, fdr.threshold, database, ontology, GOresults, outputPath){
#   # create a named vector containing the adjusted P-value as entries and the ensembl ids as names
#   geneList <- setNames(DEResults$padj, DEResults$gene_id)
#   
#   # function to indicate DE genes
#   topDiffGenes <- function(allScore, threshold=fdr.threshold){
#     return(allScore<threshold)
#   }
#   
#   # crate topGO object
#   tgd <- new("topGOdata", 
#              ontology=ontology,
#              allGenes=geneList, 
#              geneSel=topDiffGenes,
#              nodeSize=5,
#              annot=annFUN.org,
#              mapping=database,
#              ID="ensembl")
#   
#   # GO graph plot
#   GOpvalues <- GOresults$p.adjust
#   names(GOpvalues) <- GOresults$ID
#   GOpvalues <- GOpvalues[names(GOpvalues)%in%usedGO(tgd)]
#   
#   pdf(paste0(outputPath, "_GOgraphPlot.pdf"))
#   if(nrow(GOresults) > 20){
#     showSigOfNodes(tgd, GOpvalues, firstSigNodes = 20, useInfo = 'all')
#   } else {
#     showSigOfNodes(tgd, GOpvalues, firstSigNodes = nrow(GOresults), useInfo = 'all')
#   }
#   dev.off()
# }


runGSEA <- function(DESeq2Table, dbTranslationTable, fdr.threshold, database, kegg_organism, nbr.nodes, gseKEGG_table, gseKEGG_plot){

  library(clusterProfiler)
  library(ggplot2)
  library(RColorBrewer)
  library(gplots)
  library(scales)
  library(pathview)
  library(msigdbr)
  
  options(bitmapType='cairo')
  
  fdr.threshold <- as.numeric(fdr.threshold)
  nbr.nodes <- as.integer(nbr.nodes)
  
  # get all database names/ organisms
  dbTranslation <- read.table(dbTranslationTable, header=TRUE, sep='\t', quote="")
  if(database == ""){
    if(kegg_organism != "" && kegg_organism %in% dbTranslation$kegg_organism){
      database <- dbTranslation[dbTranslation$kegg_organism == kegg_organism, "topgo_db"]
    }
  }
  if(kegg_organism == ""){
    if(database != "" && database %in% dbTranslation$topgo_db){
      kegg_organism <- dbTranslation[dbTranslation$topgo_db == database, "kegg_organism"]
    }
  }
  msigdb_species <- ""
  if(kegg_organism != "" && kegg_organism %in% dbTranslation$kegg_organism){
    msigdb_species <- dbTranslation[dbTranslation$kegg_organism == kegg_organism, "msigdb_species"]
  }
  if(msigdb_species == "" && database != "" && database %in% dbTranslation$topgo_db) {
    msigdb_species <- dbTranslation[dbTranslation$topgo_db == database, "msigdb_species"]
  }
  

  DEResults <- read.table(DESeq2Table, header=TRUE, sep='\t', quote="")
  
  #remove duplicates, sort according stat and omit NA values
  DEResults <- DEResults[!is.na(DEResults$stat),]
  #remove duplicates
  DEResults <- DEResults[!duplicated(DEResults$entrezgene),]
  DEResults <- DEResults[order(DEResults$stat, decreasing = TRUE),]
  
  #remove genes without entrez gene id
  DEResults <- DEResults[!is.na(DEResults$entrezgene),]
  
  # Create a vector of the gene universe
  gene_list <- DEResults$stat
  
  # Name vector with ENTREZ ids
  names(gene_list) <- DEResults$entrezgene
  
  # store working directory
  wd <- getwd()
  
  # GSE with KEGG database ------------------------
  print("GSEA with KEGG -------------")
  tryCatch(
    {
      kegg <- gseKEGG(geneList     = gene_list,
                      organism     = kegg_organism,
                      keyType       = "ncbi-geneid",
                      pvalueCutoff = fdr.threshold,
                      pAdjustMethod = "BH")
      
      resultsKegg <- kegg@result
      resultsKegg <- resultsKegg[resultsKegg$pvalue < fdr.threshold,]
      resultsKegg <- resultsKegg[, c("ID", "Description", "setSize", "enrichmentScore", "pvalue", "p.adjust")]
      
      write.table(resultsKegg, gseKEGG_table, sep="\t", quote=FALSE, row.names = FALSE)
      
      pdf(gseKEGG_plot, width = 10, height = 10)
        ids <- resultsKegg$Description
        if(length(ids) > nbr.nodes){
          ids <- ids[1:nbr.nodes]
        }
        print(dotplot(kegg, showCategory=ids, title=paste0("Enriched KEGG Pathways (top ", nbr.nodes, ")") , split=".sign") + 
                facet_grid(.~.sign))
      dev.off()
        
    },
    error=function(e){
      #write empty table and plot
      write.table(data.frame(ID=character(), Description=character(), setSize=numeric(), enrichmentScore=numeric(), 
                             pvalue=numeric(), p.adjust=numeric()), gseKEGG_table, sep="\t", quote=FALSE, row.names = FALSE)
      
      pdf(gseKEGG_plot, width = 10, height = 10)
        plot.new()
      dev.off()
      
      return()
    }
  )
  
  tryCatch(
    {
      # Pathview plot --------
      if(nrow(resultsKegg) > 0){
        log2FoldChange <- DEResults$log2FoldChange
        names(log2FoldChange) <- DEResults$entrezgene
        
        dir.create(paste0(dirname(gseKEGG_plot), "/KEGGpathways"), showWarnings = FALSE)
        setwd(dirname(gseKEGG_plot))
        
        for(i in 1:nrow(resultsKegg)){
          pv.out <- pathview(gene.data=log2FoldChange, pathway.id=resultsKegg$ID[i],
                             species=kegg_organism, out.suffix=substr(basename(gseKEGG_plot), 1, nchar(basename(gseKEGG_plot))-4), kegg.native = TRUE,
                             kegg.dir="KEGGpathways", gene.idtype="entrez",
                             low="blue", mid="gray", high="red", same.layer=FALSE)
        }
        setwd(wd)
      }
    }
  )
  
  gseFolder <- dirname(dirname(dirname(gseKEGG_table)))
  contrastFolder <- basename(dirname(gseKEGG_table))
  fileBasename <- substr(basename(gseKEGG_table), 1, nchar(basename(gseKEGG_table))-8)
  
  
  # GSEA with MSigDB Hallmarks Collection -----------------------------------
  if(!is.na(msigdb_species) && msigdb_species != ""){
    print("GSEA with MSigDB Hallmarks Collection -------------")
    dir.create(paste0(gseFolder, "/MSigDB"), showWarnings = FALSE)
    dir.create(paste0(gseFolder, "/MSigDB/", contrastFolder), showWarnings = FALSE)
    gseMsigdb_table <- paste0(gseFolder, "/MSigDB/", contrastFolder, "/", fileBasename, "MSigDB_hallmark.txt")
    gseMsigdb_plot <- paste0(gseFolder, "/MSigDB/", contrastFolder, "/", fileBasename, "MSigDB_hallmark.pdf")
    
    tryCatch(
      {
        m_df <- msigdbr(species = msigdb_species, category = "H")[, c("gs_name", "entrez_gene")]
        
        msigdb <- GSEA(geneList = gene_list, 
                       TERM2GENE = m_df,
                       pvalueCutoff = fdr.threshold,
                       pAdjustMethod = "BH")
        
        resultsMsigdb <- msigdb@result
        resultsMsigdb <- resultsMsigdb[resultsMsigdb$pvalue < fdr.threshold,]
        resultsMsigdb <- resultsMsigdb[, c("ID", "setSize", "enrichmentScore", "pvalue", "p.adjust")]
        
        write.table(resultsMsigdb, gseMsigdb_table, sep="\t", quote=FALSE, row.names = FALSE)
        
        if(nrow(resultsMsigdb) > 0){
          ids <- resultsMsigdb$ID
          if(length(ids) > nbr.nodes){
            ids <- ids[1:nbr.nodes]
          }
          pdf(gseMsigdb_plot, width = 10, height = 10)
          print(dotplot(msigdb, showCategory=ids, title=paste0("Enriched MSigDB hallmark gene sets (top ", nbr.nodes, ")") , split=".sign") + 
                  facet_grid(.~.sign))
          dev.off()
        }
      }
    )
  }
  
  
  
  # GSE with GO database ---------------------------
  # if(!is.na(database) && database != ""){
  #   library(database, character.only = TRUE)
  #   library(topGO)
  #   print("GSEA with GO database -------------")
  
  #   dir.create(paste0(gseFolder, "/GO"), showWarnings = FALSE)
  #   dir.create(paste0(gseFolder, "/GO/", contrastFolder), showWarnings = FALSE)
  #   gseGO_table <- paste0(gseFolder, "/GO/", contrastFolder, "/", fileBasename, "GO.txt")
  #   gseGo_plot <- paste0(gseFolder, "/GO/", contrastFolder, "/", fileBasename, "GO_")
  #   
  #   write.table(data.frame(ID=character(), Description=character(), setSize=numeric(), enrichmentScore=numeric(), 
  #                          pvalue=numeric(), p.adjust=numeric(), ontology=character()), 
  #               gseGO_table, sep="\t", quote=FALSE, row.names = FALSE)
  #   
  #   onts<-c("BP", "MF", "CC")
  #   for(ont in onts){
  #     tryCatch(
  #       {
  #         print(paste0(ont, "..."))
  #         go <- gseGO(geneList = gene_list,
  #                     ont = ont,
  #                     OrgDb = get(database),
  #                     keyType = "ENTREZID",
  #                     pvalueCutoff = fdr.threshold,
  #                     pAdjustMethod = "BH")
  #         print("...done")
  #         
  #         resultsGo.ont <- go@result
  #         resultsGo.ont <- resultsGo.ont[resultsGo.ont$pvalue < fdr.threshold,]
  #         resultsGo.ont <- resultsGo.ont[, c("ID", "Description", "setSize", "enrichmentScore", "pvalue", "p.adjust")]
  #         
  #         if(nrow(resultsGo.ont) > 0){
  #           resultsGo.ont$ontology <- ont
  #           write.table(resultsGo.ont, gseGO_table, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  #           
  #           pdf(paste0(gseGo_plot, ont, ".pdf"), width = 10, height = 10)
  #             ids <- resultsGo.ont$Description
  #             if(length(ids) > nbr.nodes){
  #               ids <- ids[1:nbr.nodes]
  #             }
  #             print(dotplot(go, showCategory=ids, title=paste0("Enriched GO terms (", ont,", top ", nbr.nodes, ")") , split=".sign") + 
  #                     facet_grid(.~.sign))
  #           dev.off()
  #           
  #           # write GO graph plot
  #           makeGOPlots(DEResults, fdr.threshold, database, ont, resultsGo.ont, paste0(gseGo_plot, ont))
  #         }
  #       }
  #     )
  #   }
  # }
}


###########################################################################################


runGSEA(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])

# runGSEA(snakemake@input[["deResults"]], snakemake@input[["dbTranslationTable"]], 
#         snakemake@params[["fdr_threshold"]], snakemake@params[["topgo_db"]], snakemake@params[["kegg_organism"]], snakemake@params[["nodes"]], 
#         snakemake@output[["gseKEGG_table"]], snakemake@output[["gseKEGG_plot"]])

