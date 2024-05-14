# Bioinformatics Analysis: Maize root bacteria degrade host-specialized metabolites through the lactonase BxdA

This repository complements the bioinformatics analyses from the publication
"Maize root bacteria degrade host-specialized metabolites through the lactonase
BxdA" (Thoenen et al. 2024).

The subfolder `figures` contains the R code to produce the figure 4 of the
paper. The subfolder `pipelines` contains the code for the bullkRNA-Seq and 
the kmer pipelines.

# Kmer Analysis 

## Background

This pipeline was used to identify genes that are potentially important in the
metabolisation of MBOA to AMPO in maize root bacteria. The basic idea is to find
shared kmers between strains that are able to metabolize MBOA and AMPO, but that
are not present in strains that cannot metabolize this.

The analysis is also docuemented in the paper:



## Analysis Documentation

The strains were sequenced either with Illumina NovaSeq600 or PacBio SequelIIe and assembled using an internal genome assembly pipeline (documented in the respective repositories). The genomes were then annotated using PGAP and integrated into OpenGenomeBrowser.


### 1. Prerequisites

#### Data

For this analysis, the following data is needed per strain:

- genome fasta file (`.fna`), which contains the assembled genomes (contigs).
  This file should be stored in the folder `data`.

- gene file (`.ffn`), where the annotated nucleotide sequences of genes are
  stored in fasta. This file should be stored in `database_genes_ogb`.


#### BLAST database

Generation of BLAST database of all Microbacteria strains was performed using the script `pipelines/kmer_microbacteria/database_genes_ogb/create_database.sh`.

Before running the script, the genes files (`.ffn`) from all strains were
concatenated.

```
conda activate blast_2.10.1;
cd pipelines/kmer_microbacteria/database_genes_ogb/;
./create_database.sh;
```

### 2. Run Snakemake pipeline

Run the snakemake pipeline. This script runs until the rule `split_csv`.

### 3. Score kmers

Unfortunately, the implementation of extracting the kmer score does not work with slurm due to a misconfiguration for multiprocessing. This has not been resolved so far.

Therefore, the kmer scoring has to be started manually. In order to recreate the 
conda environment that is used for this step, the `yaml` file in
`config/conda/kmer_python_pandas.yaml` can be used.

If you use a slumr HPC, run:

```
srun --pty --mem=1200000 -c 79 --time=14:00:00 /bin/bash
conda activate kmer_python_pandas
python workflow/scripts/score_kmers_mpc.py --csv_file_path results/02_kmertable/ \ 
--phenotype_file phenotype.csv \
--output_path results/03_kmerscores/ --log_file log_manual_start.txt
```

### 4. Extract high-scoring kmers

The high-scoring kmers were extracted using the script `workflow/scripts/extract_scores.sh`.

In this case, the value 7 was considered high-scoring.

```
cd results/03_kmerscores;
../../workflow/scripts/extract_scores.sh -t $(pwd) -s 7;
```

### 5. Convert kmers to fasta file

```
cd `results/03_kmerscores
sed 's/\(.*\),\(.*\)/\2,\1/g' scores_ge_7.csv | awk '{print ">"NR"_"$0}' | sed 's/,/\n/g' | sed 's/^split.*csv://g' > scores_ge_7.fasta  
```

### 6. Run the script `search_kmers.py`

This script extracts the kmers that are matching genes.

```
conda activate biopython;
for i in $(cat genomes.txt ); do echo "Starting with sample:" $i; python workflow/scripts/search_kmers.py  database_genes_ogb/${i}.ffn  results/03_kmerscores/scores_ge_7.fasta results/03_kmerscores/out_${i}.txt; echo "Finished";
done
```

Next, get the identifiers of the genes that were found:

```
for i in $(ls results/03_kmerscores/ | grep 'out_'); 
do cat $i | cut -d, -f2 | cut -d' ' -f1 | sort | uniq | sed 's/"//g' | grep -v 'Description' > results/03_kmerscores/${i}.ids.txt; done  
```

Fetch the genes into one bigger file for use in vsearch.

```

for i in $(ls results/03_kmerscores | grep ids); do sample=$(echo $i|  sed 's/.txt.*//g' | sed 's/out_//g');
for gene in $(cat results/03_kmerscores/${i}); do /software/bin/bioawk -c fastx -v var="$gene" ' $name~var {print ">"$name" "$comment"\n"$seq}' database_genes_ogb/${sample}.ffn >> results/03_kmerscores/${sample}_genes.fasta;
done ;
done

```

After this analysis we have the following results:

- The kmers above the score of 0. We only used the kmers >= 7 
- We also have all genes per sample in which these kmers occured in results/03_kmerscores/${sample}_genes.fasta

### 7. Cluster the genese with vsearch

The next step consists of clustering the genes with an id=0.7. This step is necessary in order to fetch genes that do not have an exact kmer match but are closely related to genes with high kmerscores.

Concatenate all genes:

```
conda activate vsearch_2.17.1;
mkdir -p results/04_vsearch ;
cat results/03_kmerscores/*_genes.fasta > results/04_vsearch/all_genes.fasta;
vsearch --cluster_fast results/04_vsearch/all_genes.fasta --consout results/04_vsearch/consensus_seqs.fasta --clusters results/04_vsearch/cluster.fasta --id 0.70;
```

The centroid clusters are formatted to a file conatining which genes belong to which cluster.

```
cd results/04_vsearch ;
for i  in $(ls | grep cluster.fasta); do  genes=$(grep '>' $i | sed 's/>//g'); for gene in $(echo $genes); do echo $i","${gene} ;done; done > genes_in_clusters.csv
```
However, we still need the translation between centroid and cluster number:

```
for i in $(cat consensus_seqs.fasta | grep '>'); do genes=$(echo $i | cut -d'=' -f2 | cut -d';' -f1); for gene in $(echo $genes); do centroid=$(grep ${gene} cluster.fasta*); echo $i",    "$centroid ;done; done | sed 's/^>//g' | sed 's/:>.*//g' > translation_centroid_seq_cluster.csv
```                                                       

### 8. BLAST serach the consensus sequences to the ffn database

```
conda activate blast_2.10.1;
workflow/scripts/run_blast.sh -d database_genes_ogb/microbacteria_genes.ffn -q results/04_vsearch/consensus_seqs.fasta -o results/05_blast/consensus_blast.m8 -l results/05_blast/log_blast.txt;
```

# Bulk RNA Seq Analysis

Bioinformatics steps:

1. Run through the steps explained in 1-preparation. Save the output files in
the reference folder (`features.csv`, hisat-indices & gtf file).

2. The gtf file needed to be reformatted (see `reformat_gtf.sh`)

3. Prepare the config file

4. Run the Snakemake pipeline

5. Run the R script
