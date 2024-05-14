# rseqc_inner_distance.rule
# runs the following scripts from RSeQC, only if data are PE:
# inner_distance.py	--> produces plot with inner distances

# several parameters currently set to default
# nbr read pairs used to estimate inner distance. default=1'000'000; minimum mapping quality default=30
# nbr bins and range for plotting

rule rseqc_inner_distance:
	input:
		"results/coord_sorted/{sample}.coordSorted.bam"
	output:
		innerDist="results/rseqc/rseqc.{sample}.inner_distance_plot.pdf"
	params:
		bed=config["bed"],
		R_path=config["R_path"],
		version=config["rseqc_version"]
	# envmodules:
	# 	"UHTS/Quality_control/RSeQC/2.6.4"
	resources:
		mem_mb=1000,
		hours=1
	log:
		"logs/rseqc/{sample}.log"
	shell:
		"""
		/bin/echo Job rseqc_inner_distance.{wildcards.sample} runs on host `hostname` and has ID "$SLURM_JOB_ID"

		export PATH={params.R_path}:$PATH

		module add UHTS/Quality_control/RSeQC/{params.version}

		# srun \
		
		inner_distance.py \
		-i {input} \
		-o results/rseqc/rseqc.{wildcards.sample} \
		-r {params.bed} \
		2>> {log}
		"""
