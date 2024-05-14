rule extract_genes:
  input:
    SCORES = 'results/03_kmerscores/score9.fasta',
    CSV = "results/03_kmerscores/split_aa.csv",
    ANNOT = "../genome_assembly_pipeline/results/4_prokka/{sample}/{sample}.ffn"

  output:
    CSV = "results/04_gene_tables/{sample}.csv"
  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['python_threads'])

  resources:
    mem_mb = int(config['python_mem_mb']),
    runtime = 60*int(config['python_hours']),

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate biopython ;"
    "  python workflow/scripts/search_kmers.py "
    "  --annotated_fasta {input.ANNOT} "
    "  --kmer_file {input.SCORES} "
    "  --output_file {output.CSV} ;"
