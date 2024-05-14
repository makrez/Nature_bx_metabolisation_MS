# fastqc.rule
# run fastqc on each individual fastq file

rule fastqc:
	input:
		config["fastqFolder"]+"{FileBasename}{mates}."+config["extension"]
	output:
		"results/fastqc/{FileBasename}{mates}_fastqc.html",
		"results/fastqc/{FileBasename}{mates}_fastqc/Images/per_base_quality.png",
		"results/fastqc/{FileBasename}{mates}_fastqc/Images/adapter_content.png",
		"results/fastqc/{FileBasename}{mates}_fastqc/Images/per_sequence_gc_content.png"
	# envmodules:
    #     "UHTS/Quality_control/fastqc/0.11.5"
	# container:
	# 	"docker://biocontainers/fastqc:v0.11.9_cv8"
	params:
		version=config["fastqc_version"]
	threads:
		1
	resources:
		mem_mb=1000,
		hours=1
	shell:
		"""
		module add UHTS/Quality_control/fastqc/{params.version}
		mkdir -p results/fastqc
		/bin/echo Job Fastqc.{wildcards.FileBasename}{wildcards.mates} runs on host `hostname` and has ID "$SLURM_JOB_ID"

		fastqc -o results/fastqc --extract {input}
		"""
