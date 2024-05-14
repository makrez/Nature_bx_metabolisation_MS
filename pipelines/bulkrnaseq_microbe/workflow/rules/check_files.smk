# check_files.rule
# outputs the dictionary of files per sample for double-checking

rule check_files:
	input:
		list()
	output:
		"results/file_check.txt"
	resources:
		mem_mb=1000,
		hours=1
	run:
		outfile=open(str(output), 'w')
		outfile.write("The following fastq files are detected for each sample:\n")
		for x in FilesPerSample.items():
			outfile.write(str(x) + "\n")
		outfile.close()
