# Files to be backedup is based on what we agreed on here:
# https://projects.bioinformatics.unibe.ch/projects/ibu-best-practices/wiki/data-backup-storage-and-archiving

# big files are hard-linked, small files are copied to the BACKUP folder.


rule backup:
	input:
		deseq="results/DESeq/{}.dds.rds".format(config["count_table_basename"]),
		internal="results/internal_report/RNASeq_internal_report.html",
		customer="results/customer_files/extended/{}.counts.txt".format(config["count_table_basename"]),
		shiny="results/ShinyApp/data/topGoResults.rds" if(config["includeContrasts"]=="yes" and len(tablelist) > 0) else list()
	output:
		"results/BACKUP/pipeline/Snakefile"
	params:
		design = config["expDesignFile"],
		snakefile = "workflow/Snakefile",
		rules = "rules",
		config = "config"
	resources:
		mem_mb=1000,
		hours=1
	run:
		shell('''
			CWD=`pwd`
			DESEQ=`basename {input.deseq}`
			mkdir -p results/BACKUP
			mkdir -p results/BACKUP/DESeq
			mkdir -p results/BACKUP/pipeline
			ln -f $CWD/{input.deseq} $CWD/results/BACKUP/$DESEQ
			rsync -av results/internal_report results/BACKUP/
			rsync -av {params.design} results/BACKUP/pipeline/
			rsync -av workflow/rules results/BACKUP/pipeline/
			rsync -av config results/BACKUP/pipeline/
			rsync -av logs results/BACKUP/pipeline/
			rsync -av results/customer_files results/BACKUP/
		''')

		if not input.shiny == list():
			shell('rsync -av results/ShinyApp results/BACKUP/')

		shell('rsync -av {params.snakefile} results/BACKUP/pipeline/')
