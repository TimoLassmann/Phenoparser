#!/usr/bin/env bash

export PATH=$PATH:~/gemini/anaconda/bin
 export PATH=$PATH:~/bin
 . /etc/init.d/functions
 
 step() {
     echo -n "STEP: $@"
 
     STEP_OK=0
     [[ -w /tmp ]] && echo $STEP_OK > /tmp/step.$$
 }
 
 try() {
     # Check for `-b' argument to run command in the background.
     local BG=
 
     [[ $1 == -b ]] && { BG=1; shift; }
     [[ $1 == -- ]] && {       shift; }
 
     # Run the command.
     if [[ -z $BG ]]; then
         "$@"
     else
         "$@" &
     fi
 
     # Check if command failed and update $STEP_OK if so.
     local EXIT_CODE=$?
 
     if [[ $EXIT_CODE -ne 0 ]]; then
         STEP_OK=$EXIT_CODE
         [[ -w /tmp ]] && echo $STEP_OK > /tmp/step.$$
 
         if [[ -n $LOG_STEPS ]]; then
             local FILE=$(readlink -m "${BASH_SOURCE[1]}")
             local LINE=${BASH_LINENO[0]}
 
             echo "$FILE: line $LINE: Command \`$*' failed with exit code $EXIT_CODE." >> "$LOG_STEPS"
         fi
     fi
 
     return $EXIT_CODE
 }
 
 next() {
     [[ -f /tmp/step.$$ ]] && { STEP_OK=$(< /tmp/step.$$); rm -f /tmp/step.$$; }
     [[ $STEP_OK -eq 0 ]]  && echo_success || echo_failure
     echo
 
     return $STEP_OK
 }
 
 
 function file_exists ()
 {
     if ! [ "$1" ]
     then
         echo "file exists function needs an input";
         return 1;
     fi
     if [[ -f "$1" ]]; then
         return 0;
     else
         return 1;
     fi;
 }
 
 
   
          pwd=$(pwd)

          function usage()
          {
cat <<EOF
usage: $0  -i <inputdir> -d <HPOdatabase> -l <logfile>
EOF
              exit 1;
          }

          main() {
              export PATH=$PATH:~/gemini/anaconda/bin
              export PATH=$PATH:~/bin
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
            
              cleanup_docker
              
              echo "DONE."
          }
          
          function make_phenotype_tables()
          {
              local SAMPLENAME=
              local basefilename=$(basename $1)
              local INDIR=$(basename $2)
         
              step "Working on $basefilename";    
              try get_sampleID_from_vcf_file SAMPLENAME $1;
              next
         
         
              BASEDIR=$(dirname "$1") 
              BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)
         
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
              try perl ~/programs/phenolyzer/disease_annotation.pl ./$INDIR/$outname -f -p -ph -logistic -out ./$INDIR/$phenolyzeroutput/$SAMPLENAME -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
              try phenoparser readphe  --id $SAMPLENAME --pheno ./$INDIR/$phenolyzeroutput/$SAMPLENAME  --db $OUT_DATABASE


              #try docker run -v $WORKINGDIR:/data -u `stat -c "%u:%g" $WORKINGDIR` seqnextgen_phenolyzer perl /src/phenolyzer/disease_annotation.pl /data/$outname -f -p -ph -logistic -out /data/$phenolyzeroutput/$SAMPLENAME  -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
              #try phenoparser readphe  --id $SAMPLENAME --pheno $WORKINGDIR/$SAMPLENAME_phenolyzer/hpo.seed_gene_list  --db $OUT_DATABASE
              

              next
             
              
          }


          function cleanup_docker () {
              list=$(docker ps -a -f status=exited | grep seqnextgen | cut -f1 -d' ')
              if [[ ! $list ]]; then
                  echo "No docker containers found"
              else
                  docker rm -v $list
              fi
              return 0;
          }

          function get_sampleID_from_vcf_file()
          {
              local  __resultvar=$1
              local  myresult=$( bcftools query -l $2)
              eval $__resultvar="'$myresult'"
          }

          main "$@"
