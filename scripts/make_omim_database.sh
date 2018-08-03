#!/usr/bin/env bash

# source common
. $SNGSCRIPTS/common.sh

pwd=$(pwd)

usage() {
cat <<EOF
usage: $0  -d <local sql > -l <logfile> -k <omimkey>
EOF
 exit 1;
}

main() {
 LOG_STEPS=
 INDIR=$pwd/$SNGTMP
 DATABASE=
 OMIMKEY=

 while getopts l:d:k: opt
 do
     case ${opt} in
	 l) LOG_STEPS=${OPTARG};;
	 d) DATABASE=${OPTARG};;
	 k) OMIMKEY=${OPTARG};;
	 *) usage;;
     esac
 done

 if [ "${DATABASE}" = "" ]; then usage; fi
 if [ "${OMIMKEY}" = "" ]; then usage; fi
 if [ "${LOG_STEPS}" = "" ]; then usage; fi

 #
 #    Extract phenotype information - store in flatfiles... 
 #
 for file in $(find $INDIR -name *.vcf  -and ! -name "*.d.vcf"  -and ! -name "*combined*"  -type f); do  
     make_phenotype_tables $file;
 done > $LOG_STEPS 2>&1;             

 echo "$0 DONE!";

}

make_phenotype_tables() {
 local SAMPLENAME=
 local basefilename=$(basename $1)

 step "Working on $basefilename";
 try get_sampleID_from_vcf_file SAMPLENAME $1;
 next

 BASEDIR=$(dirname "$1") 
 # this needs to be checked as on Ubuntu and MacOS the following code strips off the end of the directory path and means that the files cannot be found.
 # BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

 local WORKINGDIR=$pwd/$SNGTMP
 local already_processed=
 newname="$SAMPLENAME"
 newname+="_omim.txt"

 #
 # Look for omim 
 #
 if [ -f $BASEDIR/$newname ]; then
     step "Retrieving OMIM info"
     echo " $SAMPLENAME $BASEDIR/$newname $DATABASE";
     $SNGPPBIN/phenoparser insert --id $SAMPLENAME --pheno $BASEDIR/$newname --key $OMIMKEY --db $DATABASE
     next
 fi
}

main "$@";
