# gsea.smk
# Runs gsea.r doing GSEA of GO terms and KEGG pathways. 

rule gsea:
	input:
		dbTranslationTable = "workflow/scripts/dbTranslationTable.txt",
		deResults = "results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt"
	output:
		gseKEGG_table = "results/GSEA/KEGG/{contrasts}/{contrasts}.gseKEGG.txt",
		gseKEGG_plot = "results/GSEA/KEGG/{contrasts}/{contrasts}.gseKEGG.pdf"
	params:
		R_path = config["R_path"],
		fdr_threshold = float(config["fdr_threshold"]),
		topgo_db = config["topgo_db"],
		kegg_organism = config["kegg_organism"],
		nodes = int(config["nodes"])
	resources:
		mem_mb=5000,
		hours=2
	log:
		"logs/gsea/{contrasts}.log"
	shell:
		"""
		/bin/echo GSEA runs on host `hostname` and has job ID "$SLURM_JOB_ID"

		# srun \
		
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/gsea.R \
		{input.deResults} \
		{input.dbTranslationTable} \
		{params.fdr_threshold} \
		"{params.topgo_db}" \
		"{params.kegg_organism}" \
		{params.nodes} \
		{output.gseKEGG_table} \
		{output.gseKEGG_plot} \
		2> {log}
		"""
		
#	container: 
#		config["R_container"]
#	script:
#		"../scripts/gsea.R"
