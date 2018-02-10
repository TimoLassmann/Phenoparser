#!/usr/bin/env bash

pwd=$(pwd)

usage() {
cat <<EOF
usage: $0  -i <inputdir> -d <HPOdatabase> -l <logfile>
EOF
exit 1;
}

main() {
LOG_STEPS=
OUT_DATABASE=
INDIR=

while getopts d:i:l: opt
do
  case ${opt} in
      d) OUT_DATABASE=${OPTARG};;
      i) INDIR=${OPTARG};;
      l) LOG_STEPS=${OPTARG};;
      *) usage;;
  esac
done

if [ "${OUT_DATABASE}" = "" ]; then usage; fi
if [ "${INDIR}" = "" ]; then usage; fi

for file in $(find $INDIR -name *.vcf  -and ! -name "*.d.vcf"  -and ! -name "*combined*"  -type f); do  
  make_phenotype_tables $file $INDIR;
done  

# not needed as docker not involved here
# cleanup_docker

echo "DONE."
}

make_phenotype_tables() {
local SAMPLENAME=
local basefilename=$(basename $1)
local INDIR=$(basename $2)

step "Working on $basefilename";    
try get_sampleID_from_vcf_file SAMPLENAME $1;
next

BASEDIR=$(dirname "$1") 
# this needs to be checked as on Ubuntu and MacOS the following code strips off the end of the directory path and means that the files cannot be found.
#BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

local WORKINGDIR=$pwd/$INDIR

#
# Look for HPO 
#

HPOname="$SAMPLENAME"
HPOname+="_hpo.txt"

OMIMname="$SAMPLENAME"
OMIMname+="_omim.txt"

outname="$SAMPLENAME"
outname+="_term_list.txt"

phenolyzeroutput="phenolyzer_"
phenolyzeroutput+="$SAMPLENAME"

phenoparser termlist --id $SAMPLENAME --db $OUT_DATABASE $WORKINGDIR/$HPOname  $WORKINGDIR/$OMIMname -o $WORKINGDIR/$outname

step "Retrieve HPO info"
try perl $PLBIN/disease_annotation.pl ./$INDIR/$outname -f -p -ph -logistic -out ./$INDIR/$phenolyzeroutput/$SAMPLENAME -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
try phenoparser readphe  --id $SAMPLENAME --pheno ./$INDIR/$phenolyzeroutput/$SAMPLENAME  --db $OUT_DATABASE

#try docker run -v $WORKINGDIR:/data -u `stat -c "%u:%g" $WORKINGDIR` seqnextgen_phenolyzer perl /src/phenolyzer/disease_annotation.pl /data/$outname -f -p -ph -logistic -out /data/$phenolyzeroutput/$SAMPLENAME  -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
#try phenoparser readphe  --id $SAMPLENAME --pheno $WORKINGDIR/$SAMPLENAME_phenolyzer/hpo.seed_gene_list  --db $OUT_DATABASE

next
}

main "$@"
