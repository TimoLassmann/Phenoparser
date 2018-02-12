#!/usr/bin/env bash

cleanup_docker () {
    list=$( $SNGDOCKERCMD ps -a -f status=exited | grep seqnextgen | cut -f1 -d' ' )
    if [[ ! $list ]]; then
        echo "No docker containers found".
    else
        $SNGDOCKERCMD rm -v $list
    fi
    return 0;
}

run_vt_normalize () {
    if ! [ "$2" ]
    then
        echo "run_vt_normalize needs a working directory";
        return 1;
    fi

    if ! [ "$3" ]
    then
        echo "run_vt_normalize needs input vcf file";
        return 1;
    fi

    local  __resultname=$1
    local myresultname=$(basename "$3" | cut -d. -f1).d.n.vcf
    if file_exists $2/$myresultname; then
        echo "$myresultname exists already.";
    elif file_exists $2/$myresultname.gz; then
        echo "$myresultname exists already.";
    else
        echo "$myresultname does not exists. Running normalise";
        $SNGDOCKERCMD run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt normalize -r /genome/hg19/hg19.fa.gz -o /data/$myresultname /data/$3
    fi

    eval $__resultname="'$myresultname'"
}

run_vt_decompose () {
    if ! [ "$2" ]
    then
        echo "run_vt_decompose needs a working directory";
        return 1;
    fi

    if ! [ "$3" ]
    then
        echo "run_vt_decompose needs input vcf file";
        return 1;
    fi

    local  __resultname=$1
    local myresultname=$(basename "$3" | cut -d. -f1).d.vcf

    if file_exists $2/$myresultname; then
        echo "$myresultname exists already.";
    else
        echo "$myresultname does not exists. Running decompose";
        $SNGDOCKERCMD run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt decompose -s /data/$3 -o /data/$myresultname
    fi

    eval $__resultname="'$myresultname'"
}

index_vcf() {

    if ! [ "$1" ]
    then
        echo "index needs and input file";
        return 1;
    fi
    step "Index vcf file"
    try $SNGHTSLIB/bgzip -f -c $1/$2 > $1/$2.gz
    try $SNGHTSLIB/tabix -p vcf -f $1/$2.gz;
    try $SNGGRABIX/grabix index $1/$2.gz;
    next
    
    return 0;
}

get_sampleID_from_vcf_file(){
    local  __resultvar=$1
    local  myresult=$( $SNGBCFTOOLS/bcftools query -l $2)
    eval $__resultvar="'$myresult'"
}


vt_pipeline (){
    if ! [ "$1" ]
    then
        echo "Pipeline function needs an input vcf file";
        return 1;
    fi

    local __innamevar=$2

    local OUTNAME=
    local INNAMEVAR=
    
    local SAMPLENAME=
    local WORKINGDIR=$3
    
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
    INNAMEVAR=$basefilename
    
    step "Run vt decompose";
    try run_vt_decompose OUTNAME $WORKINGDIR  $INNAMEVAR
    next
    
    # Swap output / input name 
    INNAMEVAR=$OUTNAME
    OUTNAME=
    
    echo "$INNAMEVAR in $OUTNAME  ";

    step "Run vt normalize";
    try run_vt_normalize OUTNAME $WORKINGDIR  $INNAMEVAR
    next

    # Swap output / input name 
    INNAMEVAR=$OUTNAME
    OUTNAME=

    echo "$INNAMEVAR in $OUTNAME  ";
    eval $__innamevar="'$INNAMEVAR'";
}


# thanks to http://tech.franzone.blog/2008/08/25/bash-script-init-style-status-message/
# Column number to place the status message
RES_COL=60
# Command to move out to the configured column number
MOVE_TO_COL="echo -en \\033[${RES_COL}G"
# Command to set the color to SUCCESS (Green)
SETCOLOR_SUCCESS="echo -en \\033[1;32m"
# Command to set the color to FAILED (Red)
SETCOLOR_FAILURE="echo -en \\033[1;31m"
# Command to set the color back to normal
SETCOLOR_NORMAL="echo -en \\033[0;39m"

# Function to print the SUCCESS status
echo_success() {
  $MOVE_TO_COL
  echo -n "["
  $SETCOLOR_SUCCESS
  echo -n $"  OK  "
  $SETCOLOR_NORMAL
  echo -n "]"
  echo -ne "\r"
  return 0
}

# Function to print the FAILED status message
echo_failure() {
  $MOVE_TO_COL
  echo -n "["
  $SETCOLOR_FAILURE
  echo -n $"FAILED"
  $SETCOLOR_NORMAL
  echo -n "]"
  echo -ne "\r"
  return 1
}

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


file_exists () {
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
