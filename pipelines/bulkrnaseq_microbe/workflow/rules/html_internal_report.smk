# html_internal_report.rule
# creates a report from the output of the various QC tools
# A table of concatenated mapping statistics will be written to mappingStats.txt. This is then used for customer report

# If filenames are expanded directly with two wildcards, the order can be inconsistent. These two functions will prevent this:
def fastqFiles(FileBasename, mates):
	if len(mates)>1:
		fileList=expand(["results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/per_base_quality.png", "results/fastqc/{FileBasename}"+mates[1]+"_fastqc/Images/per_base_quality.png"], FileBasename=FileBasename)
		fileList=sorted(fileList)
	else:
		fileList=expand("results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/per_base_quality.png", FileBasename=FileBasename)
	return fileList

def adapterFiles(FileBasename, mates):
	if len(mates)>1:
		fileList=expand(["results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/adapter_content.png", "results/fastqc/{FileBasename}"+mates[1]+"_fastqc/Images/adapter_content.png"], FileBasename=FileBasename)
		fileList=sorted(fileList)
	else:
		fileList=expand("results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/adapter_content.png", FileBasename=FileBasename)
	return fileList

def gcFiles(FileBasename, mates):
	if len(mates)>1:
		fileList=expand(["results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/per_sequence_gc_content.png", "results/fastqc/{FileBasename}"+mates[1]+"_fastqc/Images/per_sequence_gc_content.png"], FileBasename=FileBasename)
		fileList=sorted(fileList)
	else:
		fileList=expand("results/fastqc/{FileBasename}"+mates[0]+"_fastqc/Images/per_sequence_gc_content.png", FileBasename=FileBasename)
	return fileList



rule html_internal_report:
	input:
		fastqc=fastqFiles(FileBasename, mates),
		adapter=adapterFiles(FileBasename, mates),
		counts="results/counts/{}.finalCounts.noHeader.txt".format(config["count_table_basename"]),
		gc=gcFiles(FileBasename, mates),
		stats=expand("logs/hisat2/{Sample}.log", Sample=SampleNames),
		pca="results/DESeq/{}_PCA-500MostVariableGenes.pdf".format(config["count_table_basename"]),
		dist=expand("results/rseqc/{sample}.readDistributionStats.txt", sample=SampleNames),
		biotype=expand("results/biotype_plot/{sample}.CountsPerBiotype.pdf", sample=SampleNames),
		inner_dist=expand("results/rseqc/rseqc.{sample}.inner_distance_plot.pdf", sample=SampleNames) if(len(mates)==2) else list(),
		gene_body=expand("results/rseqc/rseqc.{sample}.geneBodyCoverage.curves.pdf", sample=SampleNames) if(config["includeGeneCov"]=="yes") else list(),
		sat=expand("results/rseqc/rseqc.{sample}.junctionSaturation_plot.pdf", sample=SampleNames)
	output:
		report="results/internal_report/RNASeq_internal_report.html"
	params:
		project_name=config["project_name"],
		output_dir="results/internal_report",
		group=config["group"],
		email=config["email"],
		samples=expand("{samples}", samples=SampleNames),
		fastqc_names=expand("{names}", names=FileBasename)
	resources:
		mem_mb=1000,
		hours=1
	run:
		fastqc="-q {}".format(input.fastqc)
		adapter="-a {}".format(input.adapter)
		gc="-g {}".format(input.gc)
		counts="--counts {}".format(input.counts)
		stats="--stat {}".format(input.stats)
		dist="-d {}".format(input.dist)
		inner_dist="-i {}".format(input.inner_dist)
		sat="-t {}".format(input.sat)
		pca="-z {}".format(input.pca)
		biotype="-y {}".format(input.biotype)

		if input.inner_dist==list():
			inner_dist=""
		else:
			inner_dist="-i {}".format(input.inner_dist)

		if input.gene_body==list():
			gene_body=""
		else:
			gene_body="-b {}".format(input.gene_body)

		if len(mates)==1:
			SE="-u"
		else:
			SE=""

		shell('''
		/software/bin/python \
		workflow/scripts/RNASeq_internal_report.py \
		-p {params.project_name} \
		-o {params.output_dir} \
		-c {params.group} \
		-e {params.email} \
		-n {params.samples} \
		-f {params.fastqc_names} \
		{fastqc} \
		{counts} \
		{stats} \
		{adapter} \
		{gc} \
		{dist} \
		{inner_dist} \
		{sat} \
		{gene_body} \
		{biotype} \
		{pca} \
		{SE}
		''')
