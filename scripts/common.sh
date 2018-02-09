function cleanup_docker () {
    list=$(docker ps -a -f status=exited | grep seqnextgen | cut -f1 -d' ')
    if [[ ! $list ]]; then
        echo "No docker containers found".
    else
        docker rm -v $list
    fi
    return 0;
}

function run_vt_normalize () {
    #if ! [ "$1" ]
    if ! [ "$2" ]
    then
        echo "run_vt_normalize needs a working directory";
        return 1;
    fi

    #if ! [ "$2" ]
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
        #sudo docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt normalize -r /genome/hg19/hg19.fa.gz -o /data/$myresultname /data/$3
        sudo docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt normalize -r /genome/hg19/hg19.fa.gz -o /data/$myresultname /data/$3
    fi

    eval $__resultname="'$myresultname'"
}

function run_vt_decompose () {
    #if ! [ "$1" ]
    if ! [ "$2" ]
    then
        echo "run_vt_decompose needs a working directory";
        return 1;
    fi

    #if ! [ "$2" ]
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
        #docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt decompose -s /data/$3 -o /data/$myresultname
        sudo docker run -v $2:/data -u `stat -c "%u:%g" $2` seqnextgen_vt vt decompose -s /data/$3 -o /data/$myresultname
    fi

    eval $__resultname="'$myresultname'"
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

function get_sampleID_from_vcf_file()
{
    local  __resultvar=$1
    local  myresult=$( bcftools query -l $2)
    eval $__resultvar="'$myresult'"
}


function vt_pipeline ()
{
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


