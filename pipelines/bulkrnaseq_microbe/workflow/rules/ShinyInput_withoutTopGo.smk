# ShinyInput.rule
# Creates all files needed by Shiny app


# Reformat expanded input so that it can be passed as single argument to R script:
def ReformatInput(input_list):
	temp = str(input_list)
	temp = temp.replace(" ", ",")
	return temp


# tablelist = []
# for contrast in contrasts:
# 	table = "results/TopGo/{}.topGoResults.txt".format(contrast)
# 	size = os.stat(table).st_size
# 	print(table)
# 	print(size)
# 	if size > 34:
# 		tablelist.append(table)

# if len(tablelist) > 0:
rule ShinyInput:
	input:
		rlog = "results/DESeq/Contrasts/rlog.txt",
		DESeqTables = expand("results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt", contrasts=contrasts),
		design="results/ShinyApp/data/expGroups.txt"
	output:
		toptables = "results/ShinyApp/data/topTables.rds"
	params:
		R_path = config["R_path"]
	resources:
		mem_mb=5000,
		hours=1
	run:
		DETables = ReformatInput(input.DESeqTables)

		shell('''
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/MakeShinyInput.R \
		{input.rlog} \
		NULL \
		NULL \
		{DETables} \
		NULL \
		{output.toptables}

		cp workflow/scripts/server_withoutTopGo.R results/ShinyApp/
		mv results/ShinyApp/server_withoutTopGo.R results/ShinyApp/server.R

		cp workflow/scripts/ui_withoutTopGo.R results/ShinyApp/
		mv results/ShinyApp/ui_withoutTopGo.R results/ShinyApp/ui.R
		''')