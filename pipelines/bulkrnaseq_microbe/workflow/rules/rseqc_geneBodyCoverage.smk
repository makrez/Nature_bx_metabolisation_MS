# rseqc_geneBodyCoverage.rule
# runs the following script from RSeQC:
# geneBody_coverage.py  --> number of reads along gene body


rule rseqc_geneBodyCoverage:
	input:
		"results/coord_sorted/{sample}.coordSorted.bam"
	output:
		bodyCov="results/rseqc/rseqc.{sample}.geneBodyCoverage.curves.pdf",
		rscript="results/rseqc/rseqc.{sample}.geneBodyCoverage.r",
		coverage_table="results/rseqc/rseqc.{sample}.geneBodyCoverage.txt"
	params:
		bed=config["bed"],
		R_path=config["R_path"],
		rseqc=config["rseqc_version"]
	# envmodules:
	# 	"UHTS/Quality_control/RSeQC/2.6.4"
	log:
		"logs/rseqc/{sample}.geneBodyCoverage.log"
	resources:
		mem_mb=10000,
		hours=24
	run:
		basename="results/rseqc/rseqc.{}".format(wildcards.sample)
		shell('''
		/bin/echo Job geneBodyCoverage.{wildcards.sample} with Job ID "$SLURM_JOB_ID" runs on host `hostname`

		export PATH={params.R_path}:$PATH

		module add UHTS/Quality_control/RSeQC/2.6.4

		# srun \
		
		geneBody_coverage.py \
		-i {input} \
		-r {params.bed} \
		-o {basename} \
		2> {log}
		''')
