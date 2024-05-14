# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#
#/ Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#/
#/ OPTIONS
#/
#/ EXAMPLES
#/

#{{{ CL arguments
while getopts ":f:r:i:b:d:" opt; do
  case ${opt} in
    f )
	  FOLDER_KMERS=$OPTARG
	  ;;
    r )
	  RESULTS_DIR=$OPTARG
	  ;;
    i )
    SINGULARITY_IMAGE=$OPTARG
    ;;
    b )
    SINGULARITY_BIND_DIR=$OPTARG
    ;;
    d )
    PFAM_DATABASE=$OPTARG
    ;;
    \? )
    echo "Invalid option: $OPTARG" 1>&2
	  ;;
    : )
	  echo "Invalid option: $OPTARG requires an argument" 1>&2
	  ;;
  esac
done
shift $((OPTIND -1))
#}}}

# Bash settings

set -o errexit # abort on nonzero exitstatus
set -o nounset # abort on unbound variable
set -o pipefail # dont hide errors within pipes


#{{{ Variables
readonly script_name=$(basename "${0}")
readonly script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
IFS=$'\t\n'   # Split on newlines and tabs (but not on spaces)
#}}}

main() {

# create results directory
mkdir -p $RESULTS_DIR;

# get list of kmers
KMER_FILES=$(ls ${FOLDER_KMERS} | grep 'fasta' | \
sed "s,^,${FOLDER_KMERS},g")

# loop over each kmer_file
for i in $KMER_FILES; do \
  KMER_NAME=$(echo $i | cut -d/ -f2 | cut -d. -f1)
  mkdir -p ${RESULTS_DIR}/${KMER_NAME};
  echo $i;
  (
  cd ${RESULTS_DIR}/${KMER_NAME};
  echo "Analysing kmer-file:" $i;

  singularity exec --bind $SINGULARITY_BIND_DIR $SINGULARITY_IMAGE \
  hmmer2go getorf -i ../../$i -o ${KMER_NAME}_orfs_aa.faa \
  --verbose || true;

  singularity exec --bind $SINGULARITY_BIND_DIR $SINGULARITY_IMAGE \
   hmmer2go run -i ${KMER_NAME}_orfs_aa.faa -d ../../$PFAM_DATABASE \
    -o ${KMER_NAME}_genes_orf_Pfam-A.tblout -n 4 || true; #kmer_1012_genes_

  singularity exec --bind $SINGULARITY_BIND_DIR $SINGULARITY_IMAGE \
   hmmer2go fetchmap -o pfam2go || true;

  singularity exec --bind $SINGULARITY_BIND_DIR $SINGULARITY_IMAGE \
   hmmer2go mapterms \
   -i ${KMER_NAME}_genes_orf_Pfam-A.tblout \
   -p pfam2go \
   -o ${KMER_NAME}_genes_orfs_Pfam-A_GO.tsv \
   --map || true;

   singularity exec --bind $SINGULARITY_BIND_DIR $SINGULARITY_IMAGE \
    hmmer2go map2gaf \
    -i ${KMER_NAME}_genes_orfs_Pfam-A_GO_GOterm_mapping.tsv \
    -o ${KMER_NAME}_orfs_Pfam-A_GO_GOterm_mapping.gaf \
    -s 'Bacteria' || true;
    ); done;

echo $SINGULARITY_IMAGE;

}
#{{{ Helper functions

#myfunc(){

#}

#}}}

main "${@}"
