library(gggenomes)
library(ape)
library(tidyverse)
library(gggenes)
library(ggtree)
library(patchwork)
library(here)
library(aplot)
library(cowplot)

t <- readRDS("output/tree.RDS")
gff_subs <- readRDS("output/gff_subs.RDS")
gff_subs_seqs <- readRDS("output/gff_subs_seqs.RDS") 
gene_links <- readRDS("output/gene_links.RDS")
ampo_orthogroups <- read.csv("output/positive_genes.csv", header=T) %>% 
  as_tibble()

bar_dat <- read.csv("input/Groups_microbacteria_BX_metabolization_4groups.csv", header=T) %>% 
  select(-X)

samples_plot <- t$tip.label

b <- bar_dat %>% filter(Strain %in% samples_plot)

b_wide <- b %>% pivot_wider(id_cols = Strain, 
                  names_from = compound, values_from = presence) %>% 
  rename(MBOA_HMPAA = MBOA_MHPA )

tree <- ggtree(t, branch.length = "none") +
  geom_tiplab() +
  xlim(c(0,14)) +
  scale_y_continuous(expand=c(0.01,0.7,0.01,0.7))


orthogroups <- read.table("input/Orthogroups.tsv", 
                          sep = '\t',
                          header=T) %>% as_tibble()
orthogroups <- orthogroups %>% mutate_if(is.character, str_remove_all, "gnl\\|extdb\\|") %>% 
  pivot_longer(-Orthogroup, names_to = "seq_id", values_to = "feat_id") %>% 
  separate_rows(feat_id, sep = ", ") %>% 
  mutate(seq_id = str_remove(seq_id, "\\..*$"))


gene_links2 <- gene_links %>% left_join(orthogroups, by=c("feat_id")) %>% 
  select(-seq_id.y) %>% 
  rename(seq_id = seq_id.x)

operon2_genes_LMB2 <- gene_links2 %>% 
  filter(seq_id == "LMB2") %>% 
  filter(!is.na(Orthogroup)) %>% 
  select(feat_id, Orthogroup, start) %>% 
  distinct() %>% tail(n=1) %>% pull(feat_id)

# This operon contains 8 genes
operon2_genes_LMB2 <- c(
  "LMB2-1.2_002084",
  "LMB2-1.2_002085",
  "LMB2-1.2_002086",
  "LMB2-1.2_002087",
  "LMB2-1.2_002088",
  "LMB2-1.2_002089",
  "LMB2-1.2_002090",
  "LMB2-1.2_002091"
)
opron1_genes_LMB2 <- c(
  "LMB2-1.2_002079",
  "LMB2-1.2_002080",
  "LMB2-1.2_002081",
  "LMB2-1.2_002082",
  "LMB2-1.2_002083"
)

bxda <- "LMB2-1.2_002078"

# Find these genes in the orthogroups for plot
# Operon 1
operon1_orthogroups <- list()
for (gene in opron1_genes_LMB2){
  row <- orthogroups %>% filter(
    str_detect(feat_id,gene)
  )
  operon1_orthogroups[[gene]] <- row
}
operon1_orthogroups <- bind_rows(operon1_orthogroups) %>% pull(Orthogroup)

# Operon 2
operon2_orthogroups <- list()
for (gene in operon2_genes_LMB2){
  row <- orthogroups %>% filter(
    str_detect(feat_id,gene)
  )
  operon2_orthogroups[[gene]] <- row
}
operon2_orthogroups <- bind_rows(operon2_orthogroups) %>% pull(Orthogroup)

bxda <- "LMB2-1.2_002078"
bxda_orthogroups <- list()
for (gene in bxda){
  row <- orthogroups %>% filter(
    str_detect(feat_id,gene)
  )
  bxda_orthogroups[[gene]] <- row
}
bxda_orthogroups <- bind_rows(bxda_orthogroups) %>% pull(Orthogroup)
#--------------------------------------------


# Glucosidase phenotype
gluc <- read.table("input/gluc_result.tsv", 
                   header=T) %>% as_tibble()

# DMG_ampo
dmp <- read.table("input/dmp_result.tsv", 
                  header=T) %>% as_tibble()

#MAPH
maph <-read.table("input/maph_result.tsv",
                  header=T) %>% as_tibble()

dim_mboa <-read.table("input/dim_result.tsv",
                  header=T) %>% as_tibble()


interesting_glc <- gluc %>% filter(qval <= 0.05) %>% 
  pull(Gene)

interesting_maph <- maph %>% filter(qval <= 0.05) %>% 
  pull(Gene)

interesting_dim_mboa <- dim_mboa %>% filter(qval <= 0.05) %>% 
  pull(Gene)

gene_links3 <- gene_links2 %>% 
  mutate(functional_group = case_when(
    Orthogroup %in% operon1_orthogroups ~ "AMPO_cluster",
    Orthogroup %in% operon2_orthogroups ~ "AMPO_cluster",
    Orthogroup %in% interesting_glc ~ "DIMBOA_glc",
    Orthogroup %in% interesting_maph ~ "HMPAA",
    Orthogroup %in% interesting_dim_mboa ~ "DIM_MBOA",
    Orthogroup %in% bxda_orthogroups ~ "bxdA"
  ))

gene_links4 <- gene_links3 %>% 
  mutate(start_new = end,
         end_new = start,
         start_new2 = end2,
         end_new2 = start2)



# Styling of plot
#-----------------
result <- list()
for (row in 1:nrow(gene_links4)){
  new_row <- gene_links4[row,]
  if (gene_links3[row,]$seq_id == "LWH13"){
    new_row <- new_row %>% 
      mutate(start = start_new,
             end = end_new)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LWO14"){
    new_row <- new_row %>% 
      mutate(start2 = start_new2,
             end2 = end_new2)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LMS4"){
    new_row <- new_row %>% 
      mutate(start = start_new,
             end = end_new)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LWO12"){
    new_row <- new_row %>% 
      mutate(start2 = start_new2,
             end2 = end_new2)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LWH12"){
    new_row <- new_row %>% 
      mutate(start = start_new,
             end = end_new)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LTA6"){
    new_row <- new_row %>% 
      mutate(start2 = start_new2,
             end2 = end_new2)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LWH7"){
    new_row <- new_row %>% 
      mutate(start = start_new,
             end = end_new)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LWO13"){
    new_row <- new_row %>% 
      mutate(start2 = start_new2,
             end2 = end_new2)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
  if (gene_links3[row,]$seq_id == "LTA6"){
    new_row <- new_row %>% 
      mutate(start = start_new,
             end = end_new)
    result[[row]] <- new_row
  } else {
    result[[row]] <- new_row
  }
}

gene_links5 <- bind_rows(result)  

p <- gggenomes(seqs = gff_subs_seqs, genes=gff_subs) +
  geom_seq() + geom_gene() 
  

p1 <- p %>% add_links(gene_links5)  + geom_link(aes(fill=factor(functional_group)), 
                                                size = 0, alpha=0.4) +
  scale_fill_brewer("Genes", palette="Dark2", na.value="lightgrey")

b_wide <- as.data.frame(b_wide)
rownames(b_wide) <- b_wide %>% pull(Strain)
b_wide <- b_wide %>% select(-Strain)
tree2 <- gheatmap(tree,b_wide, 
                  offset=2, 
                  low = "lightgrey", 
                  high = "darkred", 
                  width = 0.5,
                  colnames_angle = 90,
                  colnames_offset_y = -1.2,
                  colnames_offset_x = 0)    +
  theme(
    legend.position = c(0.1, 0.15),
    legend.key.height = unit(0.5, 'cm')
  ) +
  ylim(c(-2,17)) +
  scale_fill_manual(values = c("lightgrey", "darkred"))


svg(filename="output/synteny_plot_tree.svg", width=15, height=8) 
tree2
dev.off()

svg(filename="output/synteny_plot_cluster.svg", width=15, height=8) 
p1 %>% pick_by_tree(tree) %>% flip_seqs(5,6,7,8,9,12,14) + 
  theme(legend.position = c(0.8,0.85))
dev.off()

png(filename="output/synteny_plot_tree.png", width=15, height=8, res=300, units="in") 
tree2
dev.off()

png(filename="output/synteny_plot_cluster.png", width=15, height=8, res=300, units = "in") 
p1 %>% pick_by_tree(tree) %>% flip_seqs(5,6,7,8,9,12,14) + 
  theme(legend.position = c(0.8,0.85))
dev.off()
