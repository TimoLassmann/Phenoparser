#!/usr/bin/env bash

# source common
. $SNGSCRIPTS/common.sh

pwd=$(pwd)

usage() {
cat <<EOF
usage: $0  -i <patient> -g <gemini_database> -d <phenoparser database> -t <template>
EOF
  exit 1;
}

main() {
  PATIENT_ID=
  GEMINI_DATABASEPATH=
  DATABASEPATH=
  TEMPLATE=

  while getopts i:g:d:t:  opt
  do
  case ${opt} in
i) PATIENT_ID=${OPTARG};;
g) GEMINI_DATABASEPATH=${OPTARG};;
d) DATABASEPATH=${OPTARG};;
t) TEMPLATE=${OPTARG};;
*) usage;;
esac
done

if [ "${PATIENT_ID}" = "" ]; then usage; fi
if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
if [ "${DATABASEPATH}" = "" ]; then usage; fi

timestamp=$(date +"%m%d%y")

reportname="Report_"$PATIENT_ID"_"$timestamp".Rmd";

PATHNAME=$(dirname $TEMPLATE);
echo "$PATHNAME";

cat $TEMPLATE \
| sed -e "s=VARPATIENT_ID=$PATIENT_ID=g" \
| sed -e "s=VARGEMINI_DATABASEPATH=$GEMINI_DATABASEPATH=g" \
| sed -e "s=VAROMIM_DATABASEPATH=$DATABASEPATH=g" \
| sed -e "s=VARPHENOLYZER_DATABASEPATH=$PHENOLYZER_DATABASEPATH=g" \
| sed -e "s=VARPATHTOHPOOBO=$SNGHPOBO=g" \
> $reportname

R -e "rmarkdown::render('$reportname')"

echo "$timestamp $reportname";

echo "DONE."
}

main "$@"
