# featureCounts.rule
# runs featureCounts on bam file
# Chimeric reads will be excluded (-C). Will count only reads that are uniquely mapped and can be assigned unambiguously to one feature
# in principle, featureCounts wants name sorted bam but it can re-sort on the fly at little computational cost
# featureCounts output table MUST contain Geneid in column 1 and counts for each sample starting in column 7

if(len(mates)==2):
	p="-p"
else:
	p=""

rule featureCounts:
	input:
		bam = expand(featureCounts_in, sample=SampleNames),
		stats = expand(dedup_stats_in, sample=SampleNames)
	output:
		outfeatCount="results/counts/{}.featCount.txt".format(config["count_table_basename"]),
		final="results/counts/{}.finalCounts.txt".format(config["count_table_basename"]),
		noheader="results/counts/{}.finalCounts.noHeader.txt".format(config["count_table_basename"])
	params:
		version=config["featureCounts_version"],
		gtf=config["gtf"],
		strand=strand,
		feature_type=config["feature_type"],
		meta_feature=config["meta_feature"],
		min_mapQ=int(config["min_mapQ"])
	# envmodules:
	# 	"UHTS/Analysis/subread/2.0.1"
	threads:
		int(config["featureCounts_threads"])
	resources:
		mem_mb=convert_mem(config["featureCounts_mem"]),
		hours=int(config["featureCounts_time"])
	log:
		"logs/count/{}.log".format(config["count_table_basename"])
	shell:
		"""
		/bin/echo featureCounts runs on host `hostname` and has job ID "$SLURM_JOB_ID"

		module add UHTS/Analysis/subread/{params.version}

		# srun \
		
		featureCounts \
		{p} -C \
		-a {params.gtf} \
		-s {params.strand} \
		-T {threads} \
		-g {params.meta_feature} \
		-t {params.feature_type} \
		-Q {params.min_mapQ} \
		-o {output.outfeatCount} \
		{input.bam} \
		2> {log}

		tail -n+2 {output.outfeatCount} \
		| cut -f1,7- \
		| sed -e 's/results\/coord_sorted\///g' -e 's/\.coordSorted\.bam//g' \
		| sed -e 's/results\/dedup\///g' -e 's/\.deduplicated\.bam//g' \
		> {output.final}

		tail -n+3 {output.outfeatCount} \
		| cut -f1,7- > {output.noheader}''
		"""
