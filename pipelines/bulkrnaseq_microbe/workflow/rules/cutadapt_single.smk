# cutadapt.rule
# Removes poly-A tails and low quality bases
# Illumina adapters are removed at the demultiplexing stage

# Quality trimming = cut low quality read ends. --nextseq-trim is equivalent to -q 20 but will ignore qualities for Gs
# -a "A{100}" will remove stretches of A's from end of R1. Will also work when there are sequencing errors within the
# poly-A tail. Length of the poly-A stretch can be shorter than 100 and it will still be trimmed. https://cutadapt.readthedocs.io/en/stable/recipes.html#trim-poly-a-tails
# Option --revcomp is currently available only for single-end data.If given, Cutadapt searches both the read and its reverse complement for adapters.

rule cutadapt:
	input:
		reads = cutadapt_in + "{FileBasename}" + mates[0] + "." + config["extension"],
		checkfile = list() if config["kit"] == 'Truseq' else "results/extract_umi/" + "{FileBasename}" + ".extract_umi_done.txt"
	output:
		reads = "results/cutadapt/" + "{FileBasename}" + mates[0] + "." + config["extension"]
	# params:
	# 	cutadapt_version = config["cutadapt_version"]
	threads:
		4
	resources:
		mem_mb = 16000,
		hours = 12
	log:
		"logs/cutadapt/{FileBasename}.log"
	container: 
		config["cutadapt_container"]
	shell:
		'''
		/bin/echo Job cutadapt.{wildcards.FileBasename} runs on host `hostname` and has ID "$SLURM_JOB_ID" 
		
		cutadapt \
		--cores={threads} \
		-a "A{{100}}" \
		--revcomp \
		--nextseq-trim=20 \
		--minimum-length=20 \
		-o {output.reads} \
		{input.reads} \
		2>> {log}
		'''
