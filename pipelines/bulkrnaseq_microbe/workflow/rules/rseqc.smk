# rseqc.rule
# runs the following scripts from RSeQC:
# read_distribution.py	--> produces table with read distribution
# junction_saturation.py --> saturation report for splice junctions


# several parameters currently set to default
# junction_saturation.py / RPKM_saturation.py: start (default 5), end (100) and step size (5) for the resampling (argument -l, -u, -s)
# junction_saturation.py: Minimum intron size (bp). default=50 (-m)
# junction_saturation.py: Minimum number of supportting reads to call a junction. default=1 (-v)

# strandedness. Default is none (for unstranded)

strand = strand

if strand==0:
	d=""
else:
	if len(mates)==2:
		if strand==1:
			d="-d '1++,1--,2+-,2-+' "
		elif strand==2:
			d="-d '1+-,1-+,2++,2--' "
	else:
		if strand==1:
			d="-d '++,--' "
		elif strand==2:
			d="-d '+-,-+' "

# basename="results/rseqc/rseqc.{}".format(wildcards.sample)

rule rseqc:
	input:
		"results/coord_sorted/{sample}.coordSorted.bam"
	output:
		readDist="results/rseqc/{sample}.readDistributionStats.txt",
		juncSat="results/rseqc/rseqc.{sample}.junctionSaturation_plot.pdf"
	params:
		bed=config["bed"],
		R_path=config["R_path"],
		version=config["rseqc_version"],
		# strand=config["strand"]
	# envmodules:
	# 	"UHTS/Quality_control/RSeQC/2.6.4"
	resources:
		mem_mb=4000,
		hours=2
	log:
		"logs/rseqc/{sample}.log"
	shell:
		"""
		/bin/echo Job rseqc.{wildcards.sample} runs on host `hostname` and has ID "$SLURM_JOB_ID"

		export PATH={params.R_path}:$PATH

		module add UHTS/Quality_control/RSeQC/{params.version}

		# srun \

		read_distribution.py \
		-i {input} \
		-r {params.bed} \
		> {output.readDist} \
		2> {log}

		# srun \
		
		junction_saturation.py \
		-i {input} \
		-o results/rseqc/rseqc.{wildcards.sample} \
		-r {params.bed} \
		2>> {log}
		"""
