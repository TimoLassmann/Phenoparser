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
usage: $0  -i <inputdir>  -l <logfile>
EOF
    exit 1;
}

main() {
    <<pathvar>>
    LOG_STEPS=
    
    INDIR=

    NUM_THREADS=6

    
    while getopts d:i:t:l: opt
    do
        case ${opt} in
             
            i) INDIR=${OPTARG};;
            t) NUM_THREADS=${OPTARG};;
            l) LOG_STEPS=${OPTARG};;
       
            *) usage;;
        esac
    done


    if [ "${INDIR}" = "" ]; then usage; fi

    
    make_clean_copy_of_input;

    # run vt on input files... 
    for file in $(find tmp -name *.vcf -type f); do
        if ! [[ $file =~ ".d.vcf" ]]; then
            echo "$file running";            
            vt_pipeline $file;
        fi;
    done   
    #
    #    Extract phenotyp infotmation - store in flatfiles... 
    #
    #for file in $(find tmp -name *.vcf -type f); do  
    #    make_phenotype_tables $file;             
    #done  

    if file_exists sample_info.txt; then 
        rm sample_info.txt
    fi
    
    for file in $(find tmp -name *.d.n.vcf.gz -type f); do
        sanity_check_vcf_files $file;
    done
    exit;    
    echo "DONE!!! Hurrah ";

}

function make_clean_copy_of_input ()
{
    
    local WORKINGDIR=$pwd/tmp/ 
    step "Create working directory";
    try make_working_directory $WORKINGDIR
    next


    # copy all relevant files into tmp directory.
    # !!!!!rename!!!!! samples if name in vcf is Unknown!!! 
    
    step "Make local copies of vcf and annotation files"
    next
    for file in $(find $INDIR -name *.vcf -type f); do
   
        local basefilename=$(basename $file)
        local SAMPLENAME=

        step "Working with:  $basefilename";
        if file_exists $WORKINGDIR/$basefilename; then
            echo "Warning: $WORKINGDIR/$basefilename exists.";
        else
       
            try cp $file $WORKINGDIR
            try chmod 700 $WORKINGDIR/$basefilename
            local SAMPLENAME=
            get_sampleID_from_vcf_file SAMPLENAME $WORKINGDIR/$basefilename;
       
            if [ "$SAMPLENAME" = "Unknown" ]; then
                local BASEDIR=$(dirname "$file") 
                local BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)
                local DIRNAME=$(basename $BASEDIR)
                SAMPLENAME="$DIRNAME";
                echo "WARNING: Unknown sample! Will rename to: $SAMPLENAME ! ";

                echo "$SAMPLENAME" > $WORKINGDIR/samplenames.tmp;
                try bcftools reheader  -s $WORKINGDIR/samplenames.tmp $WORKINGDIR/$basefilename  > $WORKINGDIR/tmp.vcf
                try cp $WORKINGDIR/tmp.vcf $WORKINGDIR/$basefilename
                try rm $WORKINGDIR/samplenames.tmp
                try rm $WORKINGDIR/tmp.vcf
           
            else
                echo "$SAMPLENAME";
            fi
       
            BASEDIR=$(dirname "$file") 
            BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

            echo "$SAMPLENAME";
            #
            # Look for omim 
            #
            if [ -f $BASEDIR/pheno/omim.txt ]; then
                newname="$SAMPLENAME"
                newname+="_omim.txt"
           
                try cp $BASEDIR/pheno/omim.txt $WORKINGDIR/$newname
                try chmod 700 $WORKINGDIR/$newname
            fi
            #
            # Look for hpo file...
            #
            if [ -f $BASEDIR/pheno/hpo.txt ]; then
                newname="$SAMPLENAME"
                newname+="_hpo.txt"
           
                try cp $BASEDIR/pheno/hpo.txt $WORKINGDIR/$newname
                try chmod 700 $WORKINGDIR/$newname
            fi
        fi
        next
    done
}

function sanity_check_vcf_files ()
{
    if ! [ "$1" ]
    then
        echo "Sanity Check needs an input vcf file";
        return 1;
    fi
    local basefilename=$(basename $1)

    BASEDIR=$(dirname "$1") 
    BASEDIR=$(echo $BASEDIR | rev | cut -d'/' -f2- | rev)

    echo "SANITY: $1";
    local SAMPLENAME=

    get_sampleID_from_vcf_file SAMPLENAME $1;
    
    DIRNAME=$(basename $BASEDIR)

    SOURCE=$(bcftools view -h $1 |  grep "^##source=" | sed 's/^##source=//' | perl -pe 's/[\n,\t, ]+/_/g' );

    REF=$(bcftools view -h $1 | grep "^##reference=" | sed 's/^##reference=//' | perl -pe 's/[\n,\t, ]+/_/g' );

    
    printf "%s\t%s\t%s\t%s\n"  "$1" "$SAMPLENAME" "$SOURCE" "$REF" >> sample_info.txt;

}

function vt_pipeline ()
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

    INNAME=$OUTNAME
    OUTNAME=
    
    step "Indexing..."
    if file_exists $WORKINGDIR/$INNAME; then 
        try bgzip $WORKINGDIR/$INNAME
        try tabix -p vcf -f $WORKINGDIR/$INNAME.gz;
        try grabix index $WORKINGDIR/$INNAME.gz;        
    else
        echo "$WORKINGDIR/$INNAME.gz exists".
    fi
    
    next
    
    step "Docker cleanup"
    try cleanup_docker
    next 
    
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

function make_working_directory()
{
    if ! [ "$1" ]
    then
        echo "mkdir function needs an input";
        return 1;
    fi
    echo "$1";
    mkdir -p $1;
    chmod 700 $1;
    return 0;
}

function get_sampleID_from_vcf_file()
{
    local  __resultvar=$1
    local  myresult=$( bcftools query -l $2)
    eval $__resultvar="'$myresult'"
}




main "$@";
