# samtools_extractUnmapped.rule
# writes unmapped reads to a separate fasta file

rule samtools_extractUnmapped:
	input:
		"results/coord_sorted/{sample}.coordSorted.bam"
	output:
		"results/unmapped/{sample}.unmapped.fa"
	# envmodules:
	# 	"UHTS/Analysis/samtools/1.10"
	params:
		version=config["samtools_version"]
	resources:
		mem_mb=8000,
		hours=3
	log:
		"logs/unmapped/{sample}.log"
	shell:
		"""
		module add UHTS/Analysis/samtools/{params.version}
		/bin/echo Job extractUnmapped.{wildcards.sample} runs on host $(hostname) and has ID "$SLURM_JOB_ID";
		# srun \
		
		samtools view -b -f 4 {input} \
		| samtools fasta - > {output} 2> {log}
		"""
