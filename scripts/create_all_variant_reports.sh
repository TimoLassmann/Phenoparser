#!/usr/bin/env bash

# source common
. $SNGSCRIPTS/common.sh

pwd=$(pwd)

usage() {
cat <<EOF
usage: $0  -g <gemini_database> -d <phenoparser database> -t <report template>
EOF
  exit 1;
}

main() {

  PATIENT_ID=
  GEMINI_DATABASEPATH=
  DATABASEPATH=
  TEMPLATE=

  while getopts g:d:t:  opt
  do
      case ${opt} in
	  g) GEMINI_DATABASEPATH=${OPTARG};;
	  d) DATABASEPATH=${OPTARG};;
	  t) TEMPLATE=${OPTARG};;
	  *) usage;;
      esac
  done

  if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
  if [ "${DATABASEPATH}" = "" ]; then usage; fi
  if [ "${TEMPLATE}" = "" ]; then usage; fi

  samples_array=( $(gemini query -q "select name from samples" $GEMINI_DATABASEPATH) )

  for i in ${samples_array[@]}; do
      step "Creating report for $i"
      try create_variant_report.sh -i $i -g $GEMINI_DATABASEPATH -d $DATABASEPATH -t $TEMPLATE -o $SNGHPOBO
      next 
  done

  echo "DONE!!! Hurrah ";
}

main "$@"
