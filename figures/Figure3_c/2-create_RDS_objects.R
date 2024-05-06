library(gggenomes)
library(ape)
library(tidyverse)
library(gggenes)
library(here)
library(reprex)

interesting_orthogroups <- read.csv("output/positive_genes.csv", header=T) %>% 
  as_tibble() %>% 
  mutate(sample = str_replace(sample, "\\.","\\-"))

transform_gff <- function(gff_file){
  gff_file %>% 
    separate_rows(attributes, sep = ';') %>% 
    separate(attributes, "=", into=c("attribute", "value")) %>% 
    pivot_wider(values_from = value, names_from = attribute) %>% 
    mutate(gene_number = ifelse(!is.na(locus_tag), 
                                as.numeric(str_remove(locus_tag, ".*_")), NA))
}

extract_min_max <- function(Sample, 
                            interesting_orthogroups){
  gff <- read.gff(here(paste0("input/gff_files/", Sample, ".gff")))
  min_gene <- interesting_orthogroups %>% 
    filter(sample == Sample) %>% 
    filter(gene_number == min(gene_number))
  max_gene <- interesting_orthogroups %>% 
    filter(sample == Sample) %>% 
    filter(gene_number == max(gene_number))
  return(
    list(
      min=min_gene,
      max = max_gene
    )
  )
}

floor_gene_start_end <- function(subsetted_gff){
  subsetted_gff %>% 
    mutate(length = end -start) %>% 
    mutate(end = end - min(start),
           start = start - min(start))
}

# Parse blast output from Orthofinder
# This will produce the link file
#-----------------------------------------
blast_output_raw <- read.table(
   "input/blast_all_vs_all.tsv", header=F)

names(blast_output_raw) <- c("feat_id_blast", "feat_id2_blast","pident", "length", "mismatch", "gapopen", "start","end", "start2", "end2","evalue", "bitscore")

# Sequence ID
seq_ids <- read.csv("input/SequenceIDs.csv",
                           header=F) %>% as_tibble()
names(seq_ids) <- c("feat_id_blast", "feat_id", "description")

seq_ids <- seq_ids %>% 
  mutate(feat_id = str_remove(feat_id, ".*\\|.*\\|"),
         feat_id_blast = str_remove(feat_id_blast, ':'))

blast_output <- blast_output_raw %>% left_join(
  seq_ids, by = c("feat_id_blast" = "feat_id_blast")
) %>% 
  select(-c("description"), "feat_id_blast") %>% 
  left_join(seq_ids, by  = c("feat_id2_blast" = "feat_id_blast")) %>% as_tibble() %>% 
  rename(feat_id = feat_id.x,
         feat_id2 = feat_id.y)

blast_output <-  blast_output %>% select(any_of(c(def_names("blast"), "feat_id", "feat_id2"))) %>% 
  as_tibble()

samples <- interesting_orthogroups %>% pull(sample) %>% unique()

gff_subs_list <- list()
gff_subs_seqs_list <- list()
for (i in samples){
  gff <- read.gff(paste0("input/gff_files/",i,".gff"))
  gff <- transform_gff(gff)
  min <- extract_min_max(i, interesting_orthogroups)$min %>% 
    pull(gene_number)
  max <- extract_min_max(i, interesting_orthogroups)$max %>% 
    pull(gene_number)
  gff_subs <- gff %>% filter(gene_number >= min - 50  & gene_number <= max+50 ) %>% 
    filter(type == "CDS") %>% 
    rename(seq_id = seqid,
           feat_id = locus_tag) %>% 
    mutate(width = end -start)
  gff_subs <- floor_gene_start_end(gff_subs) %>% 
    rename(gene_name = gene)
  gff_subs_seqs <- data.frame(
    seq_id  = unique(gff_subs$seq_id),
    start = min(gff_subs$start),
    end = max(gff_subs$end)
  ) %>% 
    mutate(length = end - start)
  gff_subs_list[[i]] <- gff_subs
  gff_subs_seqs_list[[i]] <- gff_subs_seqs
}

gff_subs <- bind_rows(gff_subs_list) %>% 
  mutate(seq_id = str_remove(seq_id, "-.*$"))
gff_subs_seqs <- bind_rows(gff_subs_seqs_list)%>% 
  mutate(seq_id = str_remove(seq_id, "-.*$"))

genes_in_plot <- gff_subs %>% pull(feat_id)

gene_links <- blast_output %>% as_tibble() %>% 
  filter(feat_id %in% genes_in_plot) %>% 
  filter(feat_id2 %in% genes_in_plot) %>% 
  filter(feat_id != feat_id2) %>% 
  as_tibble()

### Rename the gene links with the appropriate seq_id (e.g. seq_id == {name}_scf1)
seqid_dict <- gff_subs %>% select(seq_id, feat_id)
gene_links <- gene_links %>% left_join(seqid_dict, by = c("feat_id" = "feat_id") ) %>% 
  distinct()  %>% 
  left_join(seqid_dict, by = c("feat_id2" = "feat_id")) %>% 
  distinct() %>% 
  rename(seq_id = seq_id.x,
         seq_id2 = seq_id.y)

columns_needed <- c("feat_id", "feat_id2", "pident", "length", "mismatch",
                    "gapopen", "start", "end", "start2", "end2", "evalue", "bitscore", "seq_id", "seq_id2")

gene_links <- gene_links %>% select(any_of(columns_needed))

gene_links <- gene_links %>% left_join(gff_subs %>% select(feat_id, start, end) %>% 
                           rename(start_feat_id1 = start,
                                  end_feat_id1 = end), by = c("feat_id" = "feat_id")) %>% 
  left_join(gff_subs %>% select(feat_id, start, end) %>% 
             rename(start_feat_id2 = start,
                    end_feat_id2 = end), by =  c("feat_id2" = "feat_id")) %>% 
  mutate(start = start_feat_id1,
         end = end_feat_id1,
         start2 = start_feat_id2,
         end2 = end_feat_id2)

# Read in tree
tree_file <- "input/SpeciesTree_rooted_at_outgroup_0.txt"
t <- read.tree(tree_file)
t$tip.label <- map_chr(t$tip.label, str_remove, "-.*$")
samples <- gff_subs %>% pull(seq_id) %>% unique()
t <- t %>% keep.tip(., samples)


saveRDS(gene_links, "output/gene_links.RDS")
saveRDS(gff_subs, "output/gff_subs.RDS")
saveRDS(gff_subs_seqs, "output/gff_subs_seqs.RDS")
saveRDS(t, "output/tree.RDS")
