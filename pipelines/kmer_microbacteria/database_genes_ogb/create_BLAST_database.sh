#!/bin/bash

# concatenate genes
for i in $(cat ../../../metadata/genomes.txt); do 
		sample=$(echo $i | cut -d'-' -f1); 
		ln -s ../../../organisms/${sample}/genomes/${i}/${i}.ffn . ; 
done

cat *.ffn > microbacteria_genes.ffn

makeblastdb -in microbacteria_genes.ffn -dbtype 'nucl' -parse_seqids -title "microbacteria_db";

