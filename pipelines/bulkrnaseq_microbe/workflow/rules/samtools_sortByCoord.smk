# samtools_sortByCoord.rule
# sorts bam file by coordinates and index. This file will be kept. It can be used in IGV

# ensure that subtracting memory for overhead does not lead to negative values
# for values below 3000, 50% will be used for overhead. This may not be enough

sort_mem = convert_mem(config["sort_mem"])
threads = int(config["samtools_sort_threads"])

mem_per_core = sort_mem / threads
if (mem_per_core > 3000):
	sort_mem_core = int(mem_per_core - 3000)
else:
	sort_mem_core = int(mem_per_core - 0.5 * mem_per_core)

rule samtools_sortByCoord:
	input:
		bam="results/mapped_reads/{sample}.bam",
		check_file="logs/hisat2/{sample}.checkfile.txt"
	output:
		"results/coord_sorted/{sample}.coordSorted.bam"
	# envmodules:
	# 	"UHTS/Analysis/samtools/1.10"
	params:
		version=config["samtools_version"],
		memory=sort_mem_core
	threads:
		int(config["samtools_sort_threads"])
	resources:
		mem_mb=convert_mem(config["sort_mem"]),
		hours=int(config["sort_time"])
	log:
		"logs/sortByCoord/{sample}.log"
	shell:
		"""
		module add UHTS/Analysis/samtools/{params.version}

		/bin/echo Job coordSort.{wildcards.sample} runs on host $(hostname) and has ID "$SLURM_JOB_ID" and uses memory {params.memory}

		# Make sure there are no temp files from previous runs of the pipeline:
		rm -f {wildcards.sample}.coord*.bam;
		
		# srun \
		
		samtools sort \
		-@ {threads} \
		-m {params.memory}M \
		-T {wildcards.sample}.coord \
		-o {output} \
		{input.bam} \
		2> {log}

		# srun \
		
		samtools index {output}
		"""
