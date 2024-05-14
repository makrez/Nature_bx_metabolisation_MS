# html_customer_report.rule

rule html_customer_report:
	input:
		pca="results/DESeq/{}_PCA-500MostVariableGenes.pdf".format(config["count_table_basename"]),
		mapStats="results/mappingStats.txt"
	output:
		"results/customer_files/RNASeq_analysis_report.html"
	params:
		R_path=config["R_path"],
		project_name=config["project_name"],
		group=config["group"],
		email=config["email"],
		design=config["expDesignFile"],
		fastqc_version=config["fastqc_version"],
		rseqc_version=config["rseqc_version"],
		hisat2_version=config["hisat2_version"],
		featureCounts_version=config["featureCounts_version"],
		gtf=config["gtf"],
		contrasts=config["includeContrasts"]

	resources:
		mem_mb=1000,
		hours=1

	run:
		gtf = str(params.gtf)
		gtf = gtf.split("/")[-1]
		gtf = gtf.replace('.gtf', '')

		command = "input = 'workflow/scripts/CustomerReport.rmd', "
		command = command + "output_file = '$PWD/{}', ".format(output)
		command = command + "params = list(Project = '{}', ".format(params.project_name)
		command = command + "Group = '{}', ".format(params.group)
		command = command + "fastqc_version = '{}', ".format(params.fastqc_version)
		command = command + "rseqc_version = '{}', ".format(params.rseqc_version)
		command = command + "hisat_version = '{}', ".format(params.hisat2_version)
		command = command + "featureCounts_version = '{}', ".format(params.featureCounts_version)
		command = command + "annotation = '{}', ".format(gtf)
		command = command + "contrasts = '{}', ".format(params.contrasts)
		command = command + "design = '{}', ".format(params.design)
		command = command + "mapStats = '{}', ".format(input.mapStats)
		command = command + "pca = '$PWD/{}', ".format(input.pca)
		command = command + "Analyst = '{}')".format(params.email)


		print(command)

		shell("""
		echo xvfb-run -a {params.R_path}/Rscript -e \"rmarkdown::render(knit_root_dir = \'$PWD\', {command})\"
		xvfb-run -a {params.R_path}/Rscript -e \"rmarkdown::render(knit_root_dir = \'$PWD\', {command})\"
		""")
