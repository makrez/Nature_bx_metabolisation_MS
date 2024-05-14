rule split_csv:
  input:
    "results/02_kmertable/kmer_counts.csv"
  output:
    CSV = "results/02_kmertable/split_aa",

  threads:
    int(config['short_commands_threads'])

  resources:
    mem_mb = int(config['short_commands_mb']),
    runtime = 60*int(config['short_commands_hours']),

  log:
    "results/logs/02_kmertable/concatenate.log"

  shell:
    "  workflow/scripts/equal_split_file.sh "
    "  -f {input} "
    "  -n 100000 "
    "  -o results/02_kmertable/ ;"
