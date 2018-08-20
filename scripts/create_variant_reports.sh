#!/usr/bin/env bash

# source common
. $SNGSCRIPTS/common.sh

pwd=$(pwd)

usage() {
cat <<EOF
usage: $0 -i <patient id> -g <gemini_database> -d <phenotype database> -t <report template> -l <logfile>
EOF
  exit 1;
}

main() {

  PATIENT_ID=
  GEMINI_DATABASEPATH=
  DATABASEPATH=
  TEMPLATE=
  LOG_STEPS=

  while getopts i:g:d:t:l: opt
  do
      case ${opt} in
          i) PATIENT_ID=${OPTARG};;
	  g) GEMINI_DATABASEPATH=${OPTARG};;
	  d) PHENO_DATABASEPATH=${OPTARG};;
	  t) TEMPLATE=${OPTARG};;
	  l) LOG_STEPS=${OPTARG};;
	  *) usage;;
      esac
  done

  if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
  if [ "${PHENO_DATABASEPATH}" = "" ]; then usage; fi
  if [ "${TEMPLATE}" = "" ]; then usage; fi
  if [ "${LOG_STEPS}" = "" ]; then usage; fi
  if [ ! -f "$GEMINI_DATABASEPATH" ]; then 
   echo "Your GEMINI database $GEMINI_DATABASEPATH does not exist" > $LOG_STEPS 2>&1;
   exit 1;
  fi
  if [ ! -f "$PHENO_DATABASEPATH" ]; then
   echo "Your phenotype database $PHENO_DATABASEPATH does not exist" > $LOG_STEPS 2>&1;
   exit 1;
  fi

  samples_array=();
  if [[ "$PATIENT_ID" != "" ]]; then
   patients=();
   patientID_arr=(${PATIENT_ID//,/ });
   for i in ${patientID_arr[@]}; do
    i_nospace="$(echo -e "${i}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')";
    patients+=("\""$i_nospace"\"");
   done;
   patients_str=$(IFS=, ; echo "${patients[*]}");
   samples_array=( $(gemini query -q "select name from samples where name in ($patients_str)" $GEMINI_DATABASEPATH) )
  else
   samples_array=( $(gemini query -q "select name from samples" $GEMINI_DATABASEPATH) )
  fi > $LOG_STEPS 2>&1;

  if [ ${#samples_array[@]} -eq 0 ]; then
   echo "No samples found in $GEMINI_DATABASEPATH named $PATIENT_ID" > $LOG_STEPS 2>&1;
   try exit 1; 
  fi

  for i in ${samples_array[@]}; do
      step "Creating report for $i"
      timestamp=$(date +"%m%d%y")
      
      reportname="Report_"$i"_"$timestamp".Rmd";
      
      PATHNAME=$(dirname $TEMPLATE);
      echo "$PATHNAME";
      
      cat $TEMPLATE \
      | sed -e "s=VARPATIENT_ID=$i=g" \
      | sed -e "s=VARGEMINI_DATABASEPATH=$GEMINI_DATABASEPATH=g" \
      | sed -e "s=VARPHENO_DATABASEPATH=$PHENO_DATABASEPATH=g" \
      | sed -e "s=VARPATHTOHPOOBO=$SNGHPOBO=g" \
      > $reportname
      
      R -e "rmarkdown::render('$reportname')"
      
      echo "$timestamp $reportname";
      
      next 
  done > $LOG_STEPS 2>&1;

  echo "$0 DONE!";
}

main "$@"
