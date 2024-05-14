# customer_files.rule
# sets up a directory that contains all files for the customer

rule customer_files:
	input:
		de_tables = expand("results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt", contrasts = contrasts) if config["includeContrasts"]!="none" else list(),
		topgo_tables = expand("results/TopGo/{contrasts}.topGoResults.txt", contrasts = contrasts) if config["includeContrasts"]=="full" else list(),
		gseKEGG_table = expand("results/GSEA/KEGG/{contrasts}/{contrasts}.gseKEGG.txt", contrasts = contrasts) if config["includeContrasts"]!="none" else list(),
		final="results/counts/{}.finalCounts.txt".format(config["count_table_basename"])
	output:
		counts="results/customer_files/extended/{}.counts.txt".format(config["count_table_basename"])
	params:
		contrasts=config["includeContrasts"]
	resources:
		mem_mb=1000,
		hours=1
	run:
		if input.de_tables == list():
			de_files = ""
		else:
			de_files = "--de {}".format(input.de_tables)

		shell('''
		if [ "{params.contrasts}" != "none" ]
		then
		# Add DE results and GSEA results
			mkdir -p results/customer_files/documentation;
			mkdir -p results/customer_files/de_results;
			mkdir -p results/customer_files/gsea_analysis;
			rsync workflow/scripts/RNAseq_AnalysisStepsAndOutput.pdf results/customer_files/documentation/;
			rsync workflow/scripts/GSEA_AnalysisStepsAndOutput.pdf results/customer_files/documentation/;

			xvfb-run -a /software/bin/python workflow/scripts/RNAseq_txt2excel_v2.py {de_files};
			rsync results/DESeq/Contrasts/*.xlsx results/customer_files/de_results/;
			rsync --ignore-missing-args -a --exclude KEGGpathways/ results/GSEA/ results/customer_files/gsea_analysis/;
		fi
		
		if [ "{params.contrasts}" = "full" ]
		then
		# Add TopGO results
			mkdir -p results/customer_files/go_analysis;
			mkdir -p results/customer_files/go_analysis/Plots;
			rsync workflow/scripts/GOEnrichment_AnalysisStepsAndOutput.pdf results/customer_files/documentation/;

			rsync results/TopGo/*.txt results/customer_files/go_analysis/;
			rsync --ignore-missing-args results/TopGo/Plots/*pdf results/customer_files/go_analysis/Plots/;
		fi
		

		cp {input.final} {output.counts};
		''')
