#!/usr/bin/bash
# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#
#/ Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#/ 
#/ OPTIONS
#/
#/ EXAMPLES
#/  

#{{{ CL arguments
while getopts ":f:n:o:" opt; do
  case ${opt} in
    f )
	  input_file=$OPTARG
	  ;;
	n )
	  split_number=$OPTARG
	  ;;
	o )
	  output_folder=$OPTARG
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

tail -n +2 ${input_file} | \
		split -l ${split_number} - "${output_folder}"/split_
for file in "${output_folder}"/split_*
  do
  head -n 1 ${input_file} > tmp_file
  cat "$file" >> tmp_file
  mv -f tmp_file "$file"
  done
}
#{{{ Helper functions

#myfunc(){
#}

#}}}

main "${@}"
