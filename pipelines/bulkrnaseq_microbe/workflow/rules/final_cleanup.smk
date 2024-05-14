# final_cleanup
# puts all rseqc results in a folder
# puts all slurm logs in a folder
# writes a final "all done" file

rule final_cleanup:
	input:
		report="results/internal_report/RNASeq_internal_report.html",
		customer_report="results/customer_files/RNASeq_analysis_report.html",
		file_check="results/file_check.txt",
		de_results="results/ShinyApp/data/topTables.rds" if(config["includeContrasts"]!="none" and len(tablelist) > 0) else "results/DESeq/{}.dds.rds".format(config["count_table_basename"]),
		counts="results/customer_files/extended/{}.counts.txt".format(config["count_table_basename"]),
		backup="results/BACKUP/pipeline/Snakefile",
		unmapped=expand("results/unmapped/{sample}.unmapped.fa", sample = SampleNames)
	output:
		"results/done.txt"
	threads:
		1
	resources:
		mem_mb=100,
		hours=1
	shell:
		"""
		mkdir -p logs/slurm_logs
		if [ `ls -1 slurm*.out 2>/dev/null | wc -l ` -gt 0 ]; then mv slurm*.out logs/slurm_logs/; fi;
		/bin/echo All done! > {output}
		"""
