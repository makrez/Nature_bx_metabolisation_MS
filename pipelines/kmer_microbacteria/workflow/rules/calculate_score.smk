rule calculate_score:
  input:
    "results/02_kmertable/split_aa"
  output:
    CSV = "results/03_kmerscores/split_aa.csv",
  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['kmer_score_threads'])

  resources:
    mem_mb = int(config['kmer_score_mem_mb']),
    runtime = 60*int(config['kmer_score_hours']),

  log:
    "results/logs/03_kmerscores/kmer_scoring.log"

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate python_pandas ;"
    "  python workflow/scripts/score_kmers_mpc.py "
    "  --csv_file_path results/02_kmertable/ "
    "  --phenotype_file phenotype.csv "
    "  --output_path results/03_kmerscores/ "
    "  --log_file {log} ;"
