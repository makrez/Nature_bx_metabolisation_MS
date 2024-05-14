## The script expects three arguments:

## 1) Table of counts from featureCounts. The first column must contain GeneId, all others counts for the samples
## 2) Table with geneinfo from biomart_table_biotypeplot.rule
## 3) Sample ID

## Only genes with counts > 0 will be displayed
## Levels of biotype are manually ordered

args<-commandArgs(TRUE)



#########################################################################################################
# FUNCTIONS
#########################################################################################################


# Convert the table of counts to long format and appends gene info
makeGgplotTable <- function(df, geneinfo){
  library(tidyverse)
  dataCols <- colnames(df)[-1]
  long <- gather(df, Sample, Counts, dataCols, factor_key=FALSE)
  long <- merge(long, geneinfo, by.x="Geneid", by.y="ensembl_gene_id", all.x=TRUE, all.y=FALSE)
  # Order levels of biotypes:
  long$biotypes <- factor(long$biotypes,
              levels = c("protein_coding", "lncRNA", "short_ncRNA", "rRNA", "mitochondrial", "pseudogene"))

  invisible(long)
}

makePlot <- function(df_long, sampleID){
  library(tidyverse)
  # retain only 1 sample
  temp_sample <- filter(df_long, Sample == sampleID)

  # to add nbr of observations:
  n_fun <- function(x){
    return(data.frame(y = (median(x) + (quantile(x, 0.75) - median(x)) / 2), label = paste0("n = ",length(x))))
  }

  # plotting only counts > 0
  p <- ggplot(temp_sample, aes(x=biotypes, y=Counts)) + geom_boxplot() + scale_y_log10()
  p <- p + ggtitle(sampleID) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + stat_summary(fun.data = n_fun, geom="text")
}




#########################################################################################################


data<-read.table(args[1], header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)
library(tidyverse)

# Check if the counts table contains a column called "Geneid"
if (!("Geneid" %in% colnames(data))){
    # checks if there is a column containing only entries that start with "ENS"
    ens.ids<-apply(data, 2, function(x) {length(x)==length(grep("ENS", x))})
    # if yes, this column will be given the header "Geneid"
    if(sum(ens.ids)==1){
      colnames(data)[ens.ids]<-"Geneid"
    } else {
      stop("The count table does not contain a column with ensembl IDs")
    }
  }

geneinfo <- read.table(args[2], header=TRUE, sep="\t", stringsAsFactors = FALSE)

# If the step fails, an empty plot is returned
tryCatch({
  long <- makeGgplotTable(data, geneinfo)
  p <- makePlot(long, args[3])
  ggsave(args[4], p, device="pdf")

}, error = function(e) {
 ggsave(args[4], plot.new(), device="pdf")
}
)
