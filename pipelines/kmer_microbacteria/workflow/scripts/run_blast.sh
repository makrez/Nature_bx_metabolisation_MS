# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#
#/ Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#/ 
#/ OPTIONS
#/ d: database, /path/to/db/db_prefix
#/ q: query, path/to/query.fasta
#/ o: output, path/to/output.m8
#/ l: log, path/to/log.txt
#/ EXAMPLES
#/ ./run_blast.sh -d /path/to/db/db_prefix -q query.fasta -o output.m8 -l log.txt

#{{{ CL arguments
while getopts ":d:q:o:l:" opt; do
  case ${opt} in
    d )
	  database=$OPTARG
	  ;;
    q )
	  query=$OPTARG
	  ;;
    o )
	  output=$OPTARG
	  ;;
    l )
	  log=$OPTARG
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

	blastn -query ${query} -task blastn	-db ${database} -outfmt 6 \
			-out ${output} -num_threads 12 &> ${log};
}
#{{{ Helper functions

#myfunc(){
#}

#}}}

main "${@}"
