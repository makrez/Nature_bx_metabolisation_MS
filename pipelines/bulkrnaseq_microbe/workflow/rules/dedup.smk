# runs umi_tools deduplication
# hard-wired to paired-end data

if len(mates) == 2:
	dedup_pair = '--paired'
	discard_chimeric = '--chimeric-pairs=discard'
else:
	dedup_pair = ''
	discard_chimeric = ''

rule dedup:
	input:
		"results/coord_sorted/{sample}.coordSorted.bam"
	output:
		bam = "results/dedup/{sample}.deduplicated.bam",
		stats = "results/dedup/{sample}.deduplication.stats_per_umi.tsv"
	params:
		stats = "results/dedup/{sample}.deduplication.stats",
		umitools_version = config["umitools_version"]
	threads:
		1
	resources:
		mem_mb = 32000,
		hours = 24
	log:
		"logs/dedup/{sample}.log"
#	container:
#		"docker://quay.io/biocontainers/umi_tools:1.1.2--py36hc5360cc_0"
	shell:
		'''
		/bin/echo Job dedup.{wildcards.sample} runs on host `hostname` and has ID "$SLURM_JOB_ID"
		
		module add UHTS/Analysis/UMI-tools/{params.umitools_version};
		
		 umi_tools dedup \
		 -I {input} \
		 {discard_chimeric} \
		 --temp-dir=./ \
		 {dedup_pair} \
		 --output-stats {params.stats} \
		 -S {output.bam} \
		 --no-sort-output \
		 --log {log} \
		 2>> {log}
		 '''
