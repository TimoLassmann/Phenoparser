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
usage: $0  -i <inputdir> -d <local sql > -l <logfile> -k <omimkey>
EOF
         exit 1;
     }


     main() {
         export PATH=$PATH:~/gemini/anaconda/bin
         export PATH=$PATH:~/bin
         LOG_STEPS=
         INDIR=
         DATABASE=
         OMIMKEY=
    
         while getopts i:l:d:k: opt
         do
             case ${opt} in
                 i) INDIR=${OPTARG};;
                 l) LOG_STEPS=${OPTARG};;
                 d) DATABASE=${OPTARG};;
                 k) OMIMKEY=${OPTARG};;
                 *) usage;;
             esac
         done

         if [ "${INDIR}" = "" ]; then usage; fi
         if [ "${DATABASE}" = "" ]; then usage; fi
         if [ "${OMIMKEY}" = "" ]; then usage; fi

    
         #
         #    Extract phenotyp information - store in flatfiles... 
         #
         for file in $(find $INDIR -name *.vcf  -and ! -name "*.d.vcf"  -and ! -name "*combined*"  -type f); do  
             make_phenotype_tables $file;             
         done  

         echo "DONE!!! Hurrah ";

     }



     function make_phenotype_tables()
     {
         local SAMPLENAME=
         local basefilename=$(basename $1)

         step "Working on $basefilename";
    
         try get_sampleID_from_vcf_file SAMPLENAME $1;
         next

    

         BASEDIR=$(dirname "$1") 
         BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)


         local WORKINGDIR=$pwd/tmp

         local already_processed=
         newname="$SAMPLENAME"
         newname+="_omim.txt"
    
    
         #
         # Look for omim 
         #
         if [ -f $BASEDIR/$newname ]; then
        
             step "Retrieving OMIM info"
             echo " $SAMPLENAME  $BASEDIR/$newname  $DATABASE";
             phenoparser insert --id  $SAMPLENAME --pheno  $BASEDIR/$newname --key $OMIMKEY --db $DATABASE
             next
         fi

     }


     function get_sampleID_from_vcf_file()
     {
         local  __resultvar=$1
         local  myresult=$( bcftools query -l $2)
         eval $__resultvar="'$myresult'"
     }


     main "$@";
