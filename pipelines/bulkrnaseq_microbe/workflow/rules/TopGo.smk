# TopGo.rule
# Runs TopGo using Weight01-Fisher, Weight01-KS and Classic-Fisher. Genes will be ranked based on Weight01.Fisher.
# If neither of the 3 P-values is < 0.05, the term will not be output

rule TopGo:
	input:
		"results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt"
	output:
		topgo_tables = "results/TopGo/{contrasts}.topGoResults.txt",
		genes_go = "results/TopGo/GenesInGOTerms/{contrasts}.GenesInGoTerms.rds"
	params:
		R_path = config["R_path"],
		fdr_threshold = float(config["fdr_threshold"]),
		topgo_db = config["topgo_db"],
		nodes = int(config["nodes"])
	resources:
		mem_mb=5000,
		hours=2
	log:
		"logs/TopGo/{contrasts}.log"
	shell:
		"""
		/bin/echo TopGo runs on host `hostname` and has job ID "$SLURM_JOB_ID"

		# srun \
		
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/TopGo.R \
		{input} \
		{params.fdr_threshold} \
		{params.topgo_db} \
		{params.nodes} \
		{output.topgo_tables} \
		{output.genes_go} \
		2> {log}

		mkdir -p results/TopGo/Plots

		mv -f results/DESeq/Contrasts/{wildcards.contrasts}*.pdf results/TopGo/Plots/
		"""
