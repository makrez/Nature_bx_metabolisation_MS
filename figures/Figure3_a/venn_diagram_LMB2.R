library(tidyverse)
library(ggvenn)

ortho <- read.csv("input/orthogroup_significant_gene_list.txt", header=T)
rna <- read.csv("input/condition_MBOA_vs_control_fun.csv",
                header = T)
kmer <- read.csv("input/kmer_clusters_pvalue.csv", header = T)

kmer <- kmer %>% filter(
  str_detect(gene, "LMB2"),
  kmer_fisher_qvalue < 0.05)
  
ortho <- ortho %>% 
  filter(
    str_detect(genes, "LMB2"),
    ampo_positive == TRUE
  )

rna <- rna %>% 
  filter(
    padj < 0.05
  )

list_venn <- list(
  "Orthogroup" = ortho %>% pull(genes),
  "RNA-Seq" = rna %>% pull(gene),
  "kmer" = kmer %>% pull(gene)
)


png("output/venn_diagram.png", width = 5, height = 5, res = 300, units = "in")
ggvenn(list_venn)
dev.off()

svg("output/venn_diagram.svg", width = 5, height = 5) #, res = 300, units = "in")
ggvenn(list_venn)
dev.off()

# Make a gene list
output_table <- kmer %>% full_join(rna) %>% full_join(ortho, by = c("gene" = "genes")) %>% as_tibble()

write.csv(output_table, file = "output/gene_summary.csv", row.names = F, quote = F)
