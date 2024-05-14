# DESeq2.rule
# Runs CreateDDS.R to produce results from DESeq(), rlog(), PCA plot and heatmap

rule DESeq2:
	input:
		counts="results/counts/{}.finalCounts.txt".format(config["count_table_basename"]),
		design=config["expDesignFile"],
		table="results/DESeq/geneinfo.txt"
	output:
		dds="results/DESeq/{}.dds.rds".format(config["count_table_basename"]),
		rld="results/DESeq/{}.rld.rds".format(config["count_table_basename"]),
		heatmap="results/DESeq/{}_heatmap-50MostHighlyExpressedGenes.pdf".format(config["count_table_basename"]),
		pca="results/DESeq/{}_PCA-500MostVariableGenes.pdf".format(config["count_table_basename"])
	params:
		R_path=config["R_path"],
		delim=delim,
		expdesign=config["expdesign"]
	resources:
		mem_mb=5000,
		hours=2
	log:
		"logs/Deseq/deseq.log"
	shell:
		"""
		/bin/echo DESeq runs on host `hostname` and has job ID "$SLURM_JOB_ID"

		# srun \
		
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/CreateDDS.R \
		{input.counts} \
		{input.design} \
		{input.table} \
		{params.delim} \
		{params.expdesign} \
		{output.dds} \
		{output.rld} \
		{output.heatmap} \
		{output.pca} \
		2> {log}
		"""
