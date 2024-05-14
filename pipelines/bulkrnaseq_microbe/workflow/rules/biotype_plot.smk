# biotype_plot.rule
# Creates a plot of counts per biotype

rule biotype_plot:
	input:
		counts="results/counts/{}.finalCounts.txt".format(config["count_table_basename"]),
		table="results/geneinfofile.txt"
	output:
		"results/biotype_plot/{sample}.CountsPerBiotype.pdf"
	params:
		R_path=config["R_path"],
	log:
		"logs/biotype_plot/{sample}.log"
	resources:
		mem_mb=5000,
		hours=1
	shell:
		"""
		/bin/echo Job biotype_plot.{wildcards.sample} runs on host $(hostname) and has ID "$SLURM_JOB_ID"

		# srun \
		
		{params.R_path}/Rscript --vanilla \
		workflow/scripts/biotypePlot.R \
		{input.counts} \
		{input.table} \
		{wildcards.sample} \
		{output} \
		2> {log}
		"""
