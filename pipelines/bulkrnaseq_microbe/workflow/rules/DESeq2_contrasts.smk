# DESeq2_contrasts.rule
# Runs all possible pairwise contrasts among experimental groups defined in expDesign table
# variable contrasts is set up in Snakefile

rule DESeq2_contrasts:
	input:
		dds = "results/DESeq/{}.dds.rds".format(config["count_table_basename"]),
		rlog_in = "results/DESeq/{}.rld.rds".format(config["count_table_basename"])
	output:
		de_tables = expand("results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt", contrasts = contrasts),
		design = "results/ShinyApp/data/expGroups.txt",
		rlog_out = "results/DESeq/Contrasts/rlog.txt"
	params:
		R_path = config["R_path"],
		expDesign = config["expDesignFile"],
		cutoff = config["fdr_threshold"]
	resources:
		mem_mb=1000,
		hours=2
	log:
		"logs/Deseq/deseq.log"
	run:
		for c in contrasts:
			c = c.split(".")
			shell('/bin/echo DESeq.contrasts runs on host `hostname` and has job ID "$SLURM_JOB_ID"; ')

			shell('''
			# srun \
			
			{params.R_path}/Rscript --vanilla \
			workflow/scripts/RunContrasts.R \
			{input.dds} \
			Group \
			{c[0]} \
			{c[1]} \
			{input.rlog_in} \
			{output.rlog_out} \
			{output.design} \
			{params.cutoff} \
			2>> {log}
			''')
