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
     usage: $0  -p <platform> -d <database> -l <logfile>
     EOF
              exit 1;
          }

          main(){
              export PATH=$PATH:~/gemini/anaconda/bin
              export PATH=$PATH:~/bin
         
              LOG_STEPS=
              GEMINI_DATABASE=
              PLATFORM=
         
              NUM_THREADS=6
         
         
              while getopts p:d:t:l: opt
              do
                  case ${opt} in
                      d) GEMINI_DATABASE=${OPTARG};;
                      p) PLATFORM=${OPTARG};;
                      t) NUM_THREADS=${OPTARG};;
                      l) LOG_STEPS=${OPTARG};;
                 
                      *) usage;;
                  esac
              done
         
              if [ "${GEMINI_DATABASE}" = "" ]; then usage; fi
              if [ "${PLATFORM}" = "" ]; then usage; fi
         
              samplelist=$(cat sample_info.txt |  grep $PLATFORM | cut -f 1);
         
              step "Merge all vcf files"    
              try bcftools merge ${samplelist[*]} > $pwd/tmp/combined_$PLATFORM.vcf
              next
         
              step "Run VT and VEP on combined"
              pipeline $pwd/tmp/combined_$PLATFORM.vcf;
              next
         
              step "Load into gemini";
              try load_into_gemini  $pwd/tmp combined_$PLATFORM.d.n.vep.vcf.gz
              next 

         
              echo "Done!"
          }

function pipeline()
{
    if ! [ "$1" ]
    then
        echo "Pipeline function needs an input vcf file";
        return 1;
    fi

    local OUTNAME=
    local INNAME=
    
    local SAMPLENAME=
    local WORKINGDIR=$pwd/tmp
    
    local basefilename=$(basename $1)
    
    

    # Test if sample is in existing database... 

    

    step "Working on $file $basefilename";
    try get_sampleID_from_vcf_file SAMPLENAME $WORKINGDIR/$basefilename;
    try echo "Got sample name: $SAMPLENAME";
    if [ "$SAMPLENAME" = "Unknown" ]; then
        echo "ERROR: Unknown sample!!!";
        exit 1; 
    fi
    
    
    next

    # Start pipeline with correct filename
    INNAME=$basefilename
    
    step "Run vt decompose";
    try run_vt_decompose OUTNAME $WORKINGDIR  $INNAME
    next

    
    # Swap output / inoput name 
    INNAME=$OUTNAME
    OUTNAME=
    
    echo "$INNAME in $OUTNAME  ";

    step "Run vt normalize";
    try run_vt_normalize OUTNAME $WORKINGDIR  $INNAME
    next

    # Swap output / inoput name 
    INNAME=$OUTNAME
    OUTNAME=
    
    echo "$INNAME in $OUTNAME  ";

    step "Run vep";
    try run_vep OUTNAME $WORKINGDIR  $INNAME
    next

    # Swap output / inoput name 
    INNAME=$OUTNAME
    OUTNAME=

    
    index_vcf $WORKINGDIR $INNAME

    #bcf combine al vcf ...
    
    step "Docker cleanup"
    try cleanup_docker
    next 

    return 0;
}

function cleanup_docker () {
    list=$(docker ps -a -f status=exited | grep seqnextgen | cut -f1 -d' ')
    if [[ ! $list ]]; then
        echo "No docker containers found".
    else
        docker rm -v $list
    fi
    return 0;
}

function load_into_gemini () {
    if ! [ "$1" ]
    then
        echo "index needs and input file";
        return 1;
    fi

    gemini load  --passonly   -v $1/$2 -t VEP --tempdir /home/user/tmp/ --cores $NUM_THREADS $GEMINI_DATABASE
    return 0;
}

function run_vep () {
    if ! [ "$1" ]
    then
        echo "run_vep needs a working directory";
        return 1;
    fi

    if ! [ "$2" ]
    then
        echo "run_vep needs input vcf file";
        return 1;
    fi

    local  __resultname=$1
    local myresultname=$(basename "$3" | cut -d. -f1).d.n.vep.vcf

    docker run -v $2:/data seqnextgen_vep perl /src/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl -i /data/$3 -o /data/$myresultname --vcf --fork $NUM_THREADS --offline --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --assembly GRCh37 --dir_cache /root/.vep

    eval $__resultname="'$myresultname'"
}



function run_vt_normalize () {
    if ! [ "$1" ]
    then
        echo "run_vt_normalize needs a working directory";
        return 1;
    fi

    if ! [ "$2" ]
    then
        echo "run_vt_normalize needs input vcf file";
        return 1;
    fi

    local  __resultname=$1
    local myresultname=$(basename "$3" | cut -d. -f1).d.n.vcf
    if file_exists $2/$myresultname; then
        echo "$myresultname exists.";
    elif file_exists $2/$myresultname.gz; then
        echo "$myresultname exists.";
    else
        echo "$myresultname does not exists.";
        docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt normalize -r /genome/hg19/hg19.fa.gz -o /data/$myresultname /data/$3
    fi
    
    eval $__resultname="'$myresultname'"
}


function run_vt_decompose () {
    if ! [ "$1" ]
    then
        echo "run_vt_decompose needs a working directory";
        return 1;
    fi

    if ! [ "$2" ]
    then
        echo "run_vt_decompose needs input vcf file";
        return 1;
    fi

    local  __resultname=$1
    local myresultname=$(basename "$3" | cut -d. -f1).d.vcf

    if file_exists $2/$myresultname; then
        echo "$myresultname exists.";
    else
        echo "$myresultname does not exists.";
        docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt decompose -s /data/$3 -o /data/$myresultname
    fi
    eval $__resultname="'$myresultname'"
}

function get_sampleID_from_vcf_file()
{
    local  __resultvar=$1
    local  myresult=$( bcftools query -l $2)
    eval $__resultvar="'$myresult'"
}

function index_vcf() {
    if ! [ "$1" ]
    then
        echo "index needs and input file";
        return 1;
    fi
    step "Index vcf file"
    try bgzip -f -c $1/$2 > $1/$2.gz

    try tabix -p vcf -f $1/$2.gz;
    
    try grabix index $1/$2.gz;
    next
    
    return 0;
}

main "$@";
