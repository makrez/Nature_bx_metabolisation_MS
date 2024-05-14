# The script expects two argument:
## 1) ensembl dataset to be used, e.g. hsapiens_gene_ensembl. The script does not check if the provided value is valid. Biomart will give an informative error if it isn't
## 2) Table of counts from featureCounts. The first column must contain GeneId. This will define the genes for which info is downloaded from Biomart


args<-commandArgs(TRUE)


#########################################################################################################
# FUNCTIONS
#########################################################################################################



createGeneTable <- function(ensembl_dataset, count_table) {
  ## Assumes that count_table contains a column with ensemblIDs with header "Geneid" (as produced by featureCounts)

  library(tidyverse)
  library(biomaRt)
  # take ensembl IDs for all genes in the experiment from the count table
  # get biotype from biomart
  # Note that this will use the LATEST genome version, e.g. hg20


  data <- read.table(count_table, header = TRUE, check.names = FALSE)

  # Check if the counts table contains a column called "Geneid"
  if (!("Geneid" %in% colnames(data))) {
    # checks if there is a column containing only entries that start with "ENS"
    ens.ids <-
      apply(data, 2, function(x) {
        length(x) == length(grep("ENS", x))
      })
    # if yes, this column will be given the header "Geneid"
    if (sum(ens.ids) == 1) {
      colnames(data)[ens.ids] <- "Geneid"
    } else {
      stop("The count table does not contain a column with ensembl IDs")
    }
  }


  attributes <-
    c("ensembl_gene_id",
      "gene_biotype")


  # The script will try all four available Ensembl mirrors
  mirrors <- c(UK = "www.ensembl.org",
               `US East Coast` = "useast.ensembl.org",
               `US West Coast` = "uswest.ensembl.org",
               Asia = "asia.ensembl.org")

  m <- 0
  while(!exists("geneinfo") & m < length(mirrors)) {
    m <- m+1
    tryCatch({
      print(paste("Trying ", names(mirrors)[m], " mirror...", sep=""))
      ensembl <-
        useMart("ensembl", dataset = ensembl_dataset, host = mirrors[m]) # Note that this will use the LATEST genome version, e.g. hg20
      # some attributes may not exist depending on the species
      missingColumns <-
        attributes[!is.element(attributes, listAttributes(ensembl)$name)]
      attributes <-
        attributes[is.element(attributes, listAttributes(ensembl)$name)]
      geneinfo <-
        getBM(
          attributes,
          filters = 'ensembl_gene_id',
          values = data$Geneid,
          mart = ensembl
        )
      print(paste("Ensembl ", names(mirrors)[m], " mirror was used", sep=""))
    }, error = function(e) {
      print(paste(names(mirrors)[m], " mirror not responding", sep=""))
    })
  }


  tryCatch({
    ls(geneinfo)
  }, error = function(e) {
    stop(
      "It was not possible to set up a table with gene information because none of the Ensembl mirrors could be accessed. Please try again later."
    )
  })



  # if there were missing attributes, add column containing specific value
  for (i in missingColumns) {
    switch(
      i,
      "ensembl_gene_id" = {
        stop("Table with gene info cannot be created because ensembl gene id does not exist")
      },
      "gene_biotype" = {
        geneinfo$gene_biotype <- "artifact"
        print("Gene biotype is not known. All genes will be put in category 'other' for biotype plot")
      }
    )

  }

  invisible(geneinfo)
}


simplifyBiotype<-function(df){

  library(tidyverse)

  # set up a named vector translationTable based on
  # http://www.ensembl.org/Help/Faq?id=468 and
  # https://www.gencodegenes.org/gencode_biotypes.html
  # --> C:\Users\kelle\Documents\1-Insel\21-AnalysisImportantInfo\Snakemake/RNAseq_Biotypes.xlsx
  translationTable<-data.frame(
    gene_biotype = c("protein_coding", "lncRNA", "unprocessed_pseudogene", "miRNA", "snRNA", "misc_RNA", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "snoRNA", "rRNA_pseudogene", "TEC", "transcribed_unitary_pseudogene", "transcribed_processed_pseudogene", "scaRNA", "unitary_pseudogene", "polymorphic_pseudogene", "pseudogene", "rRNA", "IG_V_pseudogene", "scRNA", "IG_C_gene", "IG_J_gene", "IG_V_gene", "sRNA", "ribozyme", "translated_processed_pseudogene", "vaultRNA", "TR_J_gene", "TR_C_gene", "TR_V_gene", "TR_V_pseudogene", "translated_unprocessed_pseudogene", "TR_D_gene", "IG_C_pseudogene", "TR_J_pseudogene", "IG_D_gene", "IG_J_pseudogene", "IG_pseudogene", "Mt_tRNA", "Mt_rRNA"),
    reduced_biotype = c("protein_coding", "lncRNA", "pseudogene", "short_ncRNA", "short_ncRNA", "short_ncRNA", "pseudogene", "pseudogene", "short_ncRNA", "pseudogene", "protein_coding", "pseudogene", "pseudogene", "short_ncRNA", "pseudogene", "pseudogene", "pseudogene", "rRNA", "pseudogene", "short_ncRNA", "protein_coding", "protein_coding", "protein_coding", "short_ncRNA", "short_ncRNA", "pseudogene", "short_ncRNA", "protein_coding", "protein_coding", "protein_coding", "pseudogene", "pseudogene", "protein_coding", "pseudogene", "pseudogene", "protein_coding", "pseudogene", "pseudogene", "mitochondrial", "mitochondrial"),
    stringsAsFactors = FALSE
  )

  df <- left_join(df, translationTable) %>%
    dplyr::select(ensembl_gene_id, reduced_biotype) %>%
    rename(biotypes = reduced_biotype)


  # if a gene_biotype does not exist in the lookup table, reduced_biotype will be NA
  invisible(df)
}



  
#########################################################################################################################################################################################################



out <- createGeneTable(args[1], args[2])
out <- simplifyBiotype(out)
outfile <- args[3]

write.table(out, outfile, row.names = FALSE, quote = FALSE, sep="\t")
