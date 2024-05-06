library(tidyverse)

fisher_results <- read.table("input/fisher_result.tsv", header = T)
significant_orthogroups <- fisher_results %>% 
  filter(qval <= 0.05) %>% 
  arrange(qval) %>% as_tibble()

ortho_genes <-  read.table("input/Orthogroups.tsv",
                          header=T, sep="\t") %>% 
  as_tibble()


best_orthogroups <- significant_orthogroups %>% 
  filter((sensitivity == 100 & specificity == 100) | (specificity == 0 & sensitivity == 0))

pos_match <- best_orthogroups %>% 
  filter(sensitivity == 100) %>% pull(Gene)

neg_match <- best_orthogroups %>% 
  filter(sensitivity == 0) %>% pull(Gene)


pos_genes <- ortho_genes %>% filter(Orthogroup %in% pos_match) %>% 
  pivot_longer(-Orthogroup, values_to = "genes", names_to = "sample") %>% 
  separate_rows(genes, sep = ",") %>% 
  mutate(genes = str_remove(genes, 'gnl\\|extdb\\|')) %>% 
  mutate(genes = str_remove(genes, " ")) %>% 
  mutate(genes = str_remove(genes, '\"')) %>% 
  mutate(gene_number = as.numeric(str_remove(genes, ".*_"))) %>% 
  filter(!is.na(gene_number)) %>% 
  arrange(sample, gene_number)


write.csv(pos_genes, "output/positive_genes.csv", quote = F, row.names=F)


neg_genes <- ortho_genes %>% filter(Orthogroup %in% neg_match) %>% 
  pivot_longer(-Orthogroup, values_to = "genes", names_to = "sample") %>% 
  separate_rows(genes, sep = ",") %>% 
  mutate(genes = str_remove(genes, 'gnl\\|extdb\\|')) %>% 
  mutate(genes = str_remove(genes, " ")) %>% 
  mutate(genes = str_remove(genes, '\"')) %>% 
  mutate(gene_number = as.numeric(str_remove(genes, ".*_"))) %>% 
  filter(!is.na(gene_number)) %>% 
  arrange(sample, gene_number)

write.csv(neg_genes, "output/negative_genes.csv", quote = F, row.names=F)

