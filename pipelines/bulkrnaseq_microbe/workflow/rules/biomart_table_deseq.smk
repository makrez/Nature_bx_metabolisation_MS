# biomart_table_deseq.rule

# Creates a table containing additional info on each gene


rule biomart_table_deseq:
	input:
		"results/counts/{}.finalCounts.txt".format(config["count_table_basename"])
	output:
		"results/DESeq/geneinfo.txt"
	params:
		ensembl_dataset=config["biomart_db"],
		R_path=config["R_path"]
	resources:
		mem_mb=5000,
		hours=2
	log:
		"logs/Deseq/biomartTable.log"
	shell:
		"""
		/bin/echo biomart_table_deseq runs on host `hostname` and has job ID "$SLURM_JOB_ID"

		# srun \

		{params.R_path}/Rscript --vanilla \
		workflow/scripts/BiomartTableDESeq.R \
		{input} \
		{params.ensembl_dataset} \
		{output} \
		2> {log}
		"""
