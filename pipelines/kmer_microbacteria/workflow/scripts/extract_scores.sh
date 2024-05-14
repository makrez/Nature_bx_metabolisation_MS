# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#
#/ Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#/ 
#/ OPTIONS
#/
#/ EXAMPLES
#/  

#{{{ CL arguments
while getopts ":t:s:o:" opt; do
  case ${opt} in
    t )
	  target_folder=$OPTARG
	  ;;
    s )
	  high_score=$OPTARG
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

		echo $script_dir;
		echo $PWD;
	(cd $target_folder;
	echo "extracting gt 0"
	grep -v ',-' *.csv | grep -v 'kmer' > scores_ge_0.csv;
	if [[ $high_score -gt 0 ]]
	then 
	echo "extracting gt 1"
	grep -v ',0' scores_ge_0.csv \
			| grep -v 'kmer' > scores_ge_1.csv;
	fi
	if [[ $high_score -gt 1 ]]
	then 
	echo "extracting gt 2"
	grep -v ',1'  scores_ge_1.csv \
		   | grep -v 'kmer' > scores_ge_2.csv;
	fi
	if [[ $high_score -gt 2 ]]
	then 
	echo "extracting gt 3"
	grep -v ',2' scores_ge_2.csv \
			| grep -v 'kmer' > scores_ge_3.csv;
	fi
	if [[ $high_score -gt 3 ]]
	then 
	echo "extracting gt 4"
	grep -v ',3'  scores_ge_3.csv \
			| grep -v 'kmer' > scores_ge_4.csv;
	fi
	if [[ $high_score -gt 4 ]]
	then 
	echo "extracting gt 5"
	grep -v ',4' scores_ge_4.csv \
			| grep -v 'kmer' > scores_ge_5.csv;
	fi
	if [[ $high_score -gt 5 ]]
	then 
	echo "extracting gt 6"
	grep -v ',5' *scores_ge_5.csv \
			| grep -v 'kmer' > scores_ge_6.csv;
	fi
	if [[ $high_score -gt 6 ]]
	then 
	echo "extracting gt 7"
	grep -v ',6' scores_ge_6.csv \
			| grep -v 'kmer' > scores_ge_7.csv;
	fi
	)
}
#{{{ Helper functions

#myfunc(){
#}

#}}}

main "${@}"
