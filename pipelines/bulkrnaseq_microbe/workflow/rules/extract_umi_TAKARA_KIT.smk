# extract_umi.rules
# Expects UMI in first 8bp of read2
# The following 3bp are adapter sequence. However, Fastqc shows that generally the base quality is low for bases 9-15 --> we discard 7bp following the UMI
# Container: sudo singularity pull --name "umitools" docker://parisepigenetics/umitools:latest

rule extract_umi:
  input:
    read1 = config["fastqFolder"] + "{FileBasename}" + mates[0] + "." + config["extension"],
    read2 = config["fastqFolder"] + "{FileBasename}" + mates[1] + "." + config["extension"]
  output:
    read1 = "results/extract_umi/" + "{FileBasename}" + mates[0] + "." + config["extension"],
    read2 = "results/extract_umi/" + "{FileBasename}" + mates[1] + "." + config["extension"],
    checkfile = "results/extract_umi/" + "{FileBasename}" + ".extract_umi_done.txt"
  params:
    umitools_version = config["umitools_version"]
  threads:
    1
  resources:
    mem_mb = 16000,
    hours = 12
  log:
    "logs/extract_umi/{FileBasename}.log"
#  container:
#    "docker://quay.io/biocontainers/umi_tools:1.1.2--py36hc5360cc_0"
  shell:
    '''
    /bin/echo Job extract_umi.{wildcards.FileBasename} runs on host `hostname` and has ID "$SLURM_JOB_ID"
    
    module add UHTS/Analysis/UMI-tools/{params.umitools_version};

    umi_tools extract \
    -I {input.read1} \
    --read2-in={input.read2} \
    --extract-method=regex \
    --bc-pattern2="(?P<umi_1>.{{8}})(?P<discard_1>.{{7}}).*" \
    -S {output.read1} \
    --read2-out={output.read2} \
    --log {log} \
    2>> {log}

    touch {output.checkfile}
    '''
