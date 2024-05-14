awk -F ';' -f reformat_gtf.awk LMB2-1.2.gtf_original | sed 's/  /\t/g' | awk -F '\t' 'BEGIN {OFS="\t"} {$9=""; print $0}' | sed 's/locus_tag/gene_id/g' | sed 's/"//g' > LMB2-1.2.gtf

