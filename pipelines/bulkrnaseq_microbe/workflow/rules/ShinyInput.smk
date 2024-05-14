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
		topGoTables = expand("results/TopGo/{contrasts}.topGoResults.txt", contrasts=contrasts),
		DESeqTables = expand("results/DESeq/Contrasts/{contrasts}.DEResults.original.rlog.txt", contrasts=contrasts),
		GenesInGoTerms = expand("results/TopGo/GenesInGOTerms/{contrasts}.GenesInGoTerms.rds", contrasts=contrasts),
		design="results/ShinyApp/data/expGroups.txt"
	output:
		topgo = "results/ShinyApp/data/topGoResults.rds",
		toptables = "results/ShinyApp/data/topTables.rds"
	params:
		R_path = config["R_path"]
	resources:
		mem_mb=5000,
		hours=1
	run:
		tablelist = []
		for contrast in contrasts:
			table = "results/TopGo/{}.topGoResults.txt".format(contrast)
			size = os.stat(table).st_size
			print(table)
			print(size)
			if size > 34:
				tablelist.append(table)

		if len(tablelist) > 0:
			topGoTables = ReformatInput(input.topGoTables)
			DETables = ReformatInput(input.DESeqTables)
			GenesInGoTerms = ReformatInput(input.GenesInGoTerms)

			shell('''
			{params.R_path}/Rscript --vanilla \
			workflow/scripts/MakeShinyInput.R \
			{input.rlog} \
			{topGoTables} \
			{output.topgo} \
			{DETables} \
			{GenesInGoTerms} \
			{output.toptables}

			cp workflow/scripts/server.R results/ShinyApp/

			cp workflow/scripts/ui.R results/ShinyApp/
			''')
		else: 
			shell('''
			touch {output.toptables}
			touch {output.topgo}

			''')
