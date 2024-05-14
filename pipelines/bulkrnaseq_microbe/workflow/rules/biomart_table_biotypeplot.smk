# biomart_table_biotypeplot.rule
# creates a table with biotype info for each gene

rule biomart_table_biotypeplot:
	input:
		"results/counts/{}.finalCounts.txt".format(config["count_table_basename"])
	output:
		"results/geneinfofile.txt"
	params:
		R_path=config["R_path"],
		biomart_db=config["biomart_db"]
	log:
		"logs/biotype_plot/makeTable.log"
	resources:
		mem_mb=5000,
		hours=1
	shell:
		"""
		/bin/echo Job biomart_table_biotypeplot runs on host $(hostname) and has ID "$SLURM_JOB_ID"

		# srun \
		
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/BiomartTableBiotypeplot.R \
		{params.biomart_db} \
		{input} \
		{output} \
		2> {log}
		"""
