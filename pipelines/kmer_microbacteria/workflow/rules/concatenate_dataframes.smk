rule concatenate:
  input:
    expand("results/01_kmercounting/merged/{sample}_kmers.txt", \
           sample = samples)
  output:
    CSV = "results/02_kmertable/kmer_counts.csv",

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['short_commands_threads'])

  resources:
    mem_mb = int(config['short_commands_mb']),
    runtime = 60*int(config['short_commands_hours']),

  log:
    "results/logs/02_kmertable/concatenate.log"

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate python_pandas ;"
    "  python workflow/scripts/join_kmers.py "
    "  --file_path results/01_kmercounting/merged/ {log} > {output.CSV} ;"
