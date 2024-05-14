## The script expects two arguments:

## 1) path to a table containing EnsemblID of all genes of interest as the first column
## 2)  ensembl dataset to be used, e.g. hsapiens_gene_ensembl. The script does
## not check if the provided value is valid. Biomart will give an informative error if it isn't

## The following output is produced:
## A table with EnsemblIDs, entrez IDs, gene symbols, chromosome, and heatmap.id for each gene of interest

  # gene_id (=ensemblID)
  # gene_symbol(s). Depending on the species, there can be multiple, e.g. that for the species
  # of interest + human. All are retrieved.
  # entrezgene
  # chromosome_name
  # heatmap.id (=gene_symbol if exists, gene_id otherwise)


args<-commandArgs(TRUE)

BiomartTableDESeq <- function(count_table, ensembl_dataset, out_table){

  library(biomaRt)

  data<-read.table(count_table, header=TRUE, sep="\t", row.names=1, check.names=FALSE)


  ## For human, gene symbols following HUGO Gene Nomenclature will be used --> hgnc_symbol
	## For mouse, MGI symbols
	## For all other species, "external_gene_name"
	if (ensembl_dataset=="hsapiens_gene_ensembl") {
		attributes<-c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', 'chromosome_name')
	} else {
		if (ensembl_dataset=="mmusculus_gene_ensembl"){
		attributes<-c('ensembl_gene_id', 'mgi_symbol', 'entrezgene_id', 'chromosome_name')
		} else {
		attributes<-c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'chromosome_name')
		}
	}

  #  ## Add additional info on genes as extra columns in the DESeqDataSet --------

  # The script will try all four available Ensembl mirrors
  mirrors <- c(UK = "www.ensembl.org",
               `US East Coast` = "useast.ensembl.org",
               `US West Coast` = "uswest.ensembl.org",
               Asia = "asia.ensembl.org")

  m <- 0

  while (!exists("genesBiomart") & m < length(mirrors)) {
    m <- m + 1
    tryCatch({
      print(paste("Trying ", names(mirrors)[m], " mirror...", sep = ""))
      ensembl <-
        useMart("ensembl", dataset = ensembl_dataset, host = mirrors[m]) # Note that this will use the LATEST genome version, e.g. hg20
      genesBiomart <-
        getBM(
          attributes,
          filters = 'ensembl_gene_id',
          values = rownames(data),
          mart = ensembl
        )
      print(paste("Ensembl ", names(mirrors)[m], " mirror was used", sep = ""))
    }, error = function(e) {
      print(paste(names(mirrors)[m], " mirror not responding", sep=""))
    })
  }

  tryCatch({
    ls(genesBiomart)
  }, error = function(e) {
    stop(
      "It was not possible to set up a table with gene information because none of the Ensembl mirrors could be accessed. Please try again later."
    )
  })

# Replace header of second column so that it is always the same across species.
  colnames(genesBiomart)[2]<-"symbol"
  colnames(genesBiomart)[3]<-"entrezgene" # Change to original header for downstream compatibility

  # replace empty cells (which occur e.g. in symbol column) by NA
  genesBiomart[genesBiomart=='']<-NA
  # remove ensembl IDs that cannot be matched unambiguously:
  genesBiomart<-genesBiomart[!(duplicated(genesBiomart$ensembl_gene_id) |
                                 duplicated(genesBiomart$ensembl_gene_id, fromLast=TRUE)),]

  #set up the annotation table to match the rows of the count table exactly
  annotation<-data.frame(gene_id=rownames(data))
  annotation<-merge(annotation, genesBiomart, by.x="gene_id", by.y="ensembl_gene_id",
                    all.x=TRUE, all.y=FALSE)
  # Add label that will be displayed in the heatmap:
  annotation$heatmap.id<-ifelse(is.na(annotation$symbol),
                                as.character(annotation$gene_id),
                                as.character(annotation$symbol))

  # there are some duplicated entries in the new column --> in these cases, replace the gene symbol with the ensembl id
  ambiguous<-rownames(as.data.frame(which(table(annotation$heatmap.id)>1)))
    for (i in 1:nrow(annotation)){
      if(is.element(annotation$heatmap.id[i], ambiguous)){
        annotation$heatmap.id[i]<-as.character(annotation$gene_id[i])
      }
    }


  write.table(annotation, out_table, quote=FALSE, row.names=FALSE, sep="\t")


}


BiomartTableDESeq(args[1], args[2], args[3])
