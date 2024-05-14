# hisat2.rule
# runs hisat2
# hisat sometimes aborted, in particular with out of time limit, without snakemake noticing. The pipeline would then complete without error using an incomplete bam file. check_file was added to avoid this.

def comma_join(files):
	return ",".join(files) if isinstance(files,list) else files

def checkMate2(wildcards):
	if(len(mates)==2):
		mate2=expand("results/cutadapt/"+"{unit}"+mates[1]+"."+config["extension"], unit=FilesPerSample[wildcards.sample])
	else:
		mate2=""
	return mate2

rule hisat2:
	input:
		mate1=lambda wildcards: expand("results/cutadapt/"+"{unit}"+mates[0]+"."+config["extension"], unit=FilesPerSample[wildcards.sample])
	output:
		bam=temp("results/mapped_reads/{sample}.bam"),
		log="logs/hisat2/{sample}.log",
		check_file="logs/hisat2/{sample}.checkfile.txt"
	params:
		ref=config["reference"],
		mate2=checkMate2,
		hisat2_version=config["hisat2_version"],
		samtools_version=config["samtools_version"],
		strand=strand
	# envmodules:
	# 	"UHTS/Aligner/hisat/2.2.1",
	# 	"UHTS/Analysis/samtools/1.10"
	log:
		"logs/hisat2/{sample}.log"
	threads:
		int(config["hisat2_threads"])
	resources:
		mem_mb=convert_mem(config["hisat2_mem"]),
		hours=int(config["hisat2_time"])
	run:
		mate1Files=comma_join(input.mate1)

		# for paired-end data:
		if len(params.mate2)>0:
			mate2Files=comma_join(params.mate2)
			inputfiles="-1 {} -2 {}".format(mate1Files, mate2Files)
			if params.strand==0:
		# The default setting of HiSat2 is unstranded
				strandness=""
			else:
				if params.strand==1:
					strandness="--rna-strandness FR"
				else:
					strandness="--rna-strandness RF"

		# for single-end data
		else:
			inputfiles="-U {}".format(mate1Files)
			if params.strand==0:
		# The default setting of HiSat2 is unstranded
				strandness=""
			else:
				if params.strand==1:
					strandness="--rna-strandness F"
				else:
					strandness="--rna-strandness R"


		shell('/bin/echo Job Hisat2.{wildcards.sample} runs on host `hostname` and has ID "$SLURM_JOB_ID"')

		shell('''
		module add UHTS/Aligner/hisat/{params.hisat2_version}
		module add UHTS/Analysis/samtools/{params.samtools_version}

		# srun \

		hisat2 -x {params.ref} {inputfiles} {strandness} -p {threads} \
		2> {log} \
		| samtools view -hbS - > {output.bam} 2>> {log}
		''')

		shell('/bin/echo "I am now really done with hisat2" > {output.check_file}')
