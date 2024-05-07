#Operon  plot
library(ape)
library(tidyverse)
library(gggenes)

gff <- read.gff("input/LMB2-1.2.gff") %>% 
  as_tibble()

genes <- read.csv("input/gene_summary.csv", header=T)

gff_filt <- gff %>% filter(type == "gene") %>% 
  separate(attributes, sep = ';', into = c("id", "gene","gbkey")) %>% 
  mutate(gene = str_remove(gene, "Name="))


gff_filt <- gff_filt %>% left_join(genes)


operon <- gff_filt %>% 
  #filter(kmer_fisher_qvalue < 0.05) %>% 
  mutate(gene_number = as.numeric(str_remove(gene, "LMB2-1.2_"))) %>% 
  filter(gene_number > 2074,
         gene_number < 2095) %>%
  mutate(orientation = case_when(
    strand == "-" ~ FALSE,
    strand == "+" ~ TRUE
  )) %>% 
  mutate(
    category = case_when(
      kmer_fisher_qvalue <= 0.05 & ampo_positive == TRUE & padj <= 0.05 ~ "sign. all methods",
      kmer_fisher_qvalue <= 0.05 & ampo_positive == TRUE & padj > 0.05 ~ "kmer & orthogroup",
      kmer_fisher_qvalue <= 0.05 & 
        (ampo_positive == FALSE | is.na(ampo_positive))  & padj < 0.05 ~ "kmer & RNA-Seq",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & 
        ampo_positive == TRUE & padj <= 0.05 ~ "orthogroup & RNA-Seq",
      kmer_fisher_qvalue <= 0.05 & 
        (ampo_positive == FALSE | is.na(ampo_positive)) & 
        padj > 0.05 ~ "kmer",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & (ampo_positive == FALSE | is.na(ampo_positive)) & padj <= 0.05 ~ "RNA-Seq",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & ampo_positive == TRUE & padj > 0.05 ~ "orthogroup",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & (ampo_positive == FALSE | is.na(ampo_positive)) & 
      (padj > 0.05 | is.na(padj)) ~ "not significant",
    )
  )

operon <- operon %>% 
  mutate(
    gene_identifier = letters[1:nrow(operon)])

plot1 <- operon %>% 
  ggplot(aes(xmin=start, xmax=end, y=seqid, forward = orientation, 
             fill =category,
         label = gene)) +
  geom_text(aes(x = (start + end)/2, label = gene_identifier), angle = 0,
            vjust = -1,
            hjust = 0.5,
            size = 3.5) +
  geom_gene_arrow() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
  ) +
  xlab("Position in Genome") 

png("output/operon_plot_LMB2.png",  width = 6, height = 2, res = 300, units="in")
plot1
dev.off()
svg("output/operon_plot_LMB2.svg",  width = 6, height = 2)
plot1
dev.off()

write.csv(operon, "output/genes_in_operon_LMB2.csv", quote=F,
          row.names = F)

operon <- gff_filt %>% 
  mutate(gene_number = as.numeric(str_remove(gene, "LMB2-1.2_"))) %>% 
  mutate(orientation = case_when(
    strand == "-" ~ FALSE,
    strand == "+" ~ TRUE
  )) %>% 
  mutate(
    category = case_when(
      kmer_fisher_qvalue <= 0.05 & ampo_positive == TRUE & padj <= 0.05 ~ "sign. all methods",
      kmer_fisher_qvalue <= 0.05 & ampo_positive == TRUE & padj > 0.05 ~ "kmer & orthogroup",
      kmer_fisher_qvalue <= 0.05 & 
        (ampo_positive == FALSE | is.na(ampo_positive))  & padj < 0.05 ~ "kmer & RNA-Seq",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & 
        ampo_positive == TRUE & padj <= 0.05 ~ "orthogroup & RNA-Seq",
      kmer_fisher_qvalue <= 0.05 & 
        (ampo_positive == FALSE | is.na(ampo_positive)) & 
        padj > 0.05 ~ "kmer",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & (ampo_positive == FALSE | is.na(ampo_positive)) & padj <= 0.05 ~ "RNA-Seq",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & ampo_positive == TRUE & padj > 0.05 ~ "orthogroup",
      (kmer_fisher_qvalue > 0.05 | is.na(kmer_fisher_qvalue)) & (ampo_positive == FALSE | is.na(ampo_positive)) & 
        (padj > 0.05 | is.na(padj)) ~ "not significant",
    )
  ) %>% 
  filter(
    category != "not significant"
  )

plot2 <- operon %>% 
  ggplot(aes(xmin=start, xmax=end, y=seqid, forward = orientation, 
             fill =category,
             colour = category)) +
  geom_gene_arrow() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "None"
  ) +
  xlab("Position in Genome") +
  geom_vline(xintercept = 2146000, colour = "red", alpha = 0.3, size = 2.1)

png("output/operon_plot_LMB2_whole_genome.png",   width = 8, height = 1, res = 300, units="in")
plot2
dev.off()

svg("output/operon_plot_LMB2_whole_genome.svg",   width = 8, height = 1)
plot2
dev.off()

