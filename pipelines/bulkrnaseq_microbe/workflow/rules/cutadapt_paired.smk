# cutadapt.rule
# Removes poly-A tails and low quality bases
# Illumina adapters are removed at the demultiplexing stage
# both reads must have a minimum length of 20bp after trimming to be retained

# Quality trimming = cut low quality read ends. --nextseq-trim is equivalent to -q 20 but will ignore qualities for Gs
# -a "A{100}" -G "T{100}" will remove stretches of A's from end of R1 and stretches of T's from beginning of R2. Will also work when there are sequencing errors within the
# poly-A tail. Length of the poly-A stretch can be shorter than 100 and it will still be trimmed. https://cutadapt.readthedocs.io/en/stable/recipes.html#trim-poly-a-tails

rule cutadapt:
	input:
		read1 = cutadapt_in + "{FileBasename}" + mates[0] + "." + config["extension"],
		read2 = cutadapt_in + "{FileBasename}" + mates[1] + "." + config["extension"],
		checkfile = list() if config["kit"] == 'Truseq' else "results/extract_umi/" + "{FileBasename}" + ".extract_umi_done.txt"
	output:
		read1 = "results/cutadapt/" + "{FileBasename}" + mates[0] + "." + config["extension"],
		read2 = "results/cutadapt/" + "{FileBasename}" + mates[1] + "." + config["extension"]
	# params:
	# 	cutadapt_version = config["cutadapt_version"]
	threads:
		4
	resources:
		mem_mb = 16000,
		hours = 12
	log:
		"logs/cutadapt/{FileBasename}.log"
	container: config["cutadapt_container"]
	shell:
		'''
		/bin/echo Job cutadapt.{wildcards.FileBasename} runs on host `hostname` and has ID "$SLURM_JOB_ID" 
		
		cutadapt \
		--cores={threads} \
		-a "A{{100}}" \
		-a "T{{100}}" \
		-G "A{{100}}" \
		-G "T{{100}}" \
		--nextseq-trim=20 \
		--minimum-length=20 \
		--pair-filter=any \
		-o {output.read1} \
		-p {output.read2} \
		{input.read1} \
		{input.read2} \
		2>> {log}
		'''
