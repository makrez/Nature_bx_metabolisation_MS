rule kmer_count:
  input:
    config['DataFolder'] + "{sample}" + config['complete_extension'],

  output:
    COUNTS = "results/01_kmercounting/{sample}_kmer_counts.kmc_pre",

  log:
    LOG="results/logs/01_kmercounting/{sample}.log",
    EXECUTION_LOG_F="results/logs/01_kmercounting/{sample}_execution_f.json",

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",
    K_VALUE = int(config['k_value']),
    tmp_dir = config['kmer_tempdir'] + "tmp_{sample}"

  threads:
    int(config['kmc']['kmc_threads'])

  resources:
    mem_mb = int(config['kmc']['kmc_mem_mb']),
    runtime = 60*int(config['kmc']['kmc_hours']),

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate kmc_3.1.1 ; newgrp p776;"
    " mkdir -p {params.tmp_dir} ;" #results/01_kmercounting/tmp_{wildcards.sample} ;"
    "  kmc " #kmc -m1 /share/data/D-ecoli-10K-ilmn/z1.fastq out1 tmp
    "  -k{params.K_VALUE} "
    "  -m16  "
    "  -fm "
    "  -ci1 "
    "  -t{threads} "
    "  -j{log.EXECUTION_LOG_F} "
    "  {input} "
    "  results/01_kmercounting/{wildcards.sample}_kmer_counts "
    "  {params.tmp_dir} ;"

rule kmer_dump:
  input:
    REVERSE = "results/01_kmercounting/{sample}_kmer_counts.kmc_pre"

  output:
    COUNTS = "results/01_kmercounting/merged/{sample}_kmers.txt"

  params:
    conda_profile = "/mnt/apps/centos7/Conda/miniconda3/etc/profile.d/conda.sh",

  threads:
    int(config['kmc']['kmc_threads'])

  resources:
    mem_mb = int(config['kmc']['kmc_mem_mb']),
    runtime = 60*int(config['kmc']['kmc_hours']),

  shell:
    " set +u ;"
    " source {params.conda_profile} ;"
    " conda activate kmc_3.1.1 ;"
    " mkdir -p results/01_kmercounting/tmp_{wildcards.sample} ;"
    "  kmc_dump " #kmc_tools simple out1 out2 union merge
    "  -t{threads} "
    "   results/01_kmercounting/{wildcards.sample}_kmer_counts "
    "   {output.COUNTS} ;"
