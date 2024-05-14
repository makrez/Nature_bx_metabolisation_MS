# parse_map_stats.rule
# Takes various stats from hisat2 log files and featureCounts summary and writes them to a table for the customer report

if len(mates)==1:
	SE="-u"
else:
	SE=""

rule parse_map_stats:
	input:
		counts="results/counts/{}.finalCounts.noHeader.txt".format(config["count_table_basename"]),
		stats=expand("logs/hisat2/{sample}.log", sample=SampleNames)
	output:
		"results/mappingStats.txt"
	params:
		samples=expand("{sample}", sample=SampleNames)
	resources:
		mem_mb=1000,
		hours=1
	shell:
		"""
		/software/bin/python workflow/scripts/ParseMappingStats.py \
		--names {params.samples} \
		--stat {input.stats} \
		--counts {input.counts} \
		--output {output} \
		{SE}
		"""
