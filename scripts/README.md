Introduction
============

The purpose of this pipeline is to discover causative variants in
patients suffering from rare disease. Currently we only focus on exome
data from various platforms.

Software installation
=====================

The pipeline applies the following steps to patient VCF files and accompanying annotation:

|  **Step**                  |  **Software**                                                                         |  **Purpose**                                                                                                                         | 
| ---------------------------| --------------------------------------------------------------------------------------| -------------------------------------------------------------------------------------------------------------------------------------|
|  **Data cleaning**         |                                                                                       |  Get all vcf files into a comparable state
|  Decompose                 |  vt ( <https://genome.sph.umich.edu/wiki/Vt> )                                        |  Decompose multiallelic variants
|  Normalize                 |  vt                                                                                   |  Standardize vcf format (see [Tan paper](https://academic.oup.com/bioinformatics/article/31/13/2202/196142) )
|  =====                     |                                                                                       |
|  **Phenotype data**        |                                                                                       |  Parse phenotype data
|  Create gene panels        |  phenoparser ( <https://github.com/TimoLassmann/Phenoparser> )                        |  A simple C program to query omim and run phenolyzer, storing all results in a database
|                            |  requires: phenolyzer ( <http://phenolyzer.wglab.org/> )                              |  Creates gene lists (panels) from input HPO terms
|                            |  requires: sqlite ( <https://www.sqlite.org/download.html> )                          |  A single database is used to store all data 
|  =====                     |                                                                                       |
|  **Variant data**          |                                                                                       |  Annotate and store/manage variant data
|  Annotation                |  Variant Effect Predictor ( <https://www.ensembl.org/info/docs/tools/vep/index.html> )|  Add annotation prior to import into Gemini
|  Storage/management        |  GEMINI ( <https://gemini.readthedocs.io/en/latest/> )                                |  Database holding variants and annotation; can be queried using sql statements
|                            |  requires: sqlite ( <https://www.sqlite.org/download.html> )                          |  We generate a separate database for each sequencing technology (IonTorrent, Illumina and Solid)
|  =====                     |                                                                                       |
|  **Reporting**             |                                                                                       |  
|  Generate patient reports  |  Shell script and R ( <https://www.r-project.org/> )                                  |  Extracts variant information and combines it with phenotype data into one report of ranked candidate causal variants
|                            |  Requires availability of Human Phenotype Ontology txt file                           |  

In addition to the scripts in this distribution (detailed below), the highlighted software components in the table above need to be installed for the pipeline to work.


Notes for particular software
-----------------------------

**Variant Effect Predictor**

VEP installs are linked to Ensembl release versions. So that HGVS expressions are properly created in the reports, be sure to download a compatible (i.e. same Ensembl release) genomic reference sequence file.
Include the path to this file in the sng_config file using the SNGVEPREF variable.
For example:
Using VEP https://github.com/Ensembl/ensembl-tools/archive/release/84.zip
Use Assembly ftp://ftp.ensembl.org/pub/grch37/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

In sng_config
SNGVEPBIN=/path/to/installed/vep
SNGVEPREF=/path/to/bgzipped/Homo_sapiens.GRCh37.dna.primary_assembly.fa

GEMINI solely supports human genetic variation mapped to build 37 (aka hg19) of the human genome so be sure to obtain GRCh37 Ensembl data.

**R**

R is used to perform data manipulation and generate the reports. In order for the template to run properly the following R packages need to be installed
dplyr
rlang
stringr
knitr
DT
tidyverse
knitr
kableExtra
ontologyIndex


Data Organisation
=================

Create a directory to store the raw patient data.
For each patient, create one directory named after the patient ID (preferably the one used in the generation of the VCF file).
In this directory, create two sub-directories: vcf and pheno. The latter contains text files with the phenotype / disease annotation.

    /home/user/patient_data/D12-3456/
    ├── pheno
    │   ├── hpo.txt
    │   └── omim.txt
    └── vcf
        └── D12-3456.vcf

Code
====

The pipeline consists of various scripts that process the input data and run the software components above.

SNG Config File
---------------

The code ships with a sample config file. This is essentially a list of environment variables denoting the location of the relevant software components

The `$PATH` variable has to contain the gemini path:

    PATH=$PATH:~/gemini/anaconda/bin

and all shell script have to he reachable:

    PATH=$PATH:~/local_install/scripts

###
# change all docker references to local installs
# remove reference to docker tidy up
# add HPO obo file location to sng_config and alter create_all_variant_reports.sh to use this environment variable rather that the -o option
Convenience library {#sec-4-2}
-------------------

This file contains some functions I regularly use:

    export PATH=$PATH:~/gemini/anaconda/bin
    export PATH=$PATH:~/bin

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

Pre-processing {#sec-4-3}
--------------

This script will copy and normalize vcf files. Run it from within a
directory you create for the analysis (e.g.
20180108$_{\text{Patient}}$$_{\text{analysis}}$).

It takes two option (run without option will print a help message):

-i &lt;directory holding data of multiple patients&gt; -l &lt;log file
to keep track of what’s happening&gt;

Example usage:

    ../scripts/run_pre_merge.sh -i ~/patient_data -l pre_merge_log

And here is the script itself:

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

The main function calling all other functions. You may need to adjust
the path variables below to fit your installation.

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

Before doing anything make a copy or the patient data. No analysis files
should end up in the same directory as the raw input data.

NOTE: in our case there are discrepancies between the patient ID used in
the vcf files (sample IDs) and the naming of the vcf and other files.
Because of this I decided to consistently use the names in the vcf
files.

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

This unfortunately named function extracts some extra information from
vcf files. This information is used later to determine whether Illumina
/ Solid / Ion torrent was used.

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

This is the pipeline to run the vt decompose and normalize. These steps
are necessary to load the data into gemini. Some annotations are also
added to variants. I also index the resulting files (if I recall
correctly there were some issues with the bcftools version - make sure
you path is set correctly).

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

This functions shuts down docker instances….

    function cleanup_docker () {
        list=$(docker ps -a -f status=exited | grep seqnextgen | cut -f1 -d' ')
        if [[ ! $list ]]; then
        echo "No docker containers found".
        else
        docker rm -v $list
        fi
        return 0;
    }

Not used any more…

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

VT normalize step run in docker…

NOTE: needs to be updated if the genome in another place than genome
hg19.

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

Runs decompose in docker…

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

Creates a new directory.

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

extracts sample ID from vcf file - very handy.

Needs bcftools installed.

    function get_sampleID_from_vcf_file()
    {
        local  __resultvar=$1
        local  myresult=$( bcftools query -l $2)
        eval $__resultvar="'$myresult'"
    }




    main "$@";

Build gemini database(s): {#sec-4-4}
-------------------------

This script combines vcf files from different sequencing technology and
loads them into gemini. I re-run vt decompose and normalization for good
measure but don’t think this is strictly necessary.

NOTE: the script will search (grep!) through the
sample$_{\text{info}}$.txt file created by
run$_{\text{pre}}$$_{\text{merge}}$.sh to identify all vcf files from
one technology. Currently we have:

-   “Life” for life technologies - SOLID

-   “Torrent” for Ion Torrent data

-   “GATK” for Illumia files

i.e. the base caller is used to identify the technology as older vcf
files are devoid of and usable meta-data.

There are three parameters:

-p “technology” from the options above -d “output database” -l log file

Example usage:

    ../scripts/build_gemini_db.sh -p GATK -d GATK.db -l gatk.log 
    ../scripts/build_gemini_db.sh -p Life -d Life.db -l gatk.log 
    ../scripts/build_gemini_db.sh -p Torrent -d Torrent.db -l gatk.log

And here is the actual script:

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

The vt pipeline

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

The command to load data into gemini. Note in the past this was done
within docker. However, this caused lots of problems because of the huge
(50GB+) annotation files that gemini requires. In future we may want to
move away from docker completely and provide an analysis image. This may
also align better with using Pawsey and other cloud resources.

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

Make Omim database {#sec-4-5}
------------------

The script below runs phenoparser written by myself to extract
information from the OMIM database. All information is stored in a sql
lite database to make sure the analysis is reproducible even if OMIM
changes. This also allows us to re-run previously un-diagnosed cases if
there is a major update.

To access OMIM you need an OMIM key which you can request online.

The options are:

-i &lt;directory holding data of multiple patients&gt; -d &lt;output
database&gt; -l &lt;log file&gt; -k &lt;OMIM key&gt;

Example usage:

    ../scripts/make_omim_database.sh  -i ~/patient_data -d omim.db -l omim.log -k
    <KEY>

And here is the actual script:

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

Phenolyzer {#sec-4-6}
----------

This script takes HPO terms **and Disease terms** for each patient and
queries Phenolyzer. For usability I now store the resuls in the same
database as the OMIM information (simply use -d &lt;same database name
as above&gt;.

The options are:

-i &lt;directory holding data of multiple patients&gt; -d &lt;output
database&gt; -l &lt;log file&gt;

Example usage:

    ../scripts/make_phenolyzer_database.sh -i <directory where the copied vcf files
    are> -d hpo.txt -l phenolyzer.log

And here is the actual script:

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

Patient reports {#sec-4-7}
---------------

I use RMarkown and knitr to create per-patient reports. A template
contains special variables that are replaced by patient details. This
gives us a lot of flexibility for the future: we can have reports for
research including more variants etc…

### Create variant report {#sec-4-7-1}

This script brings all the data together into one report. It extracts
variants together with their annotation from gemini and overlays in
silico gene panels from OMIM and Phenolyzer.

Then script works by copying a report template, replacing placeholder
variables with patient data and then running the template in R. The
output is a html file containing information about the analysis as well
as a tab separated file with the variant table.

The options are:

-i &lt;patient id&gt; -g &lt;gemini database&gt; -o &lt;omim
database&gt; -p &lt;phenolyzer database&gt; -t &lt;report template&gt;

Example usage:

    ../scripts/create_variant_report.sh -i <patientID> -g Torrent.db -o omim.db -p
    hpo.txt -t ../scripts/report_master_template.Rmd

And here is the actual script:

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
    usage: $0  -i <patient> -g <gemini_database> -d <phenoparser database> -t <template>
    EOF
          exit 1;
          }

          main() {
          export PATH=$PATH:~/gemini/anaconda/bin
          export PATH=$PATH:~/bin

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
          | sed -e "s=VARPATHTOHPOOBO=$PATHNAME=g" \
          > $reportname

          R -e "rmarkdown::render('$reportname')"


          echo "$timestamp $reportname";



          echo "DONE."
          }

          main "$@"

### Master Template {#sec-4-7-2}

    ---
    output:
        html_document:
        keep_md: true
    ---

    # Patient report: VARPATIENT_ID

    Sequencing platform: VARGEMINI_DATABASEPATH

    <style>
         .main-container { width: 1600px; max-width:1600px;}
    </style>

    ```{r setup, warning = FALSE, message = FALSE, include=FALSE}
              # Load the packages into R
              library(dplyr)
              library(rlang)
              library(stringr)
              library(knitr)
              library(DT)
              library(tidyverse)
              library(knitr)
              library(kableExtra)
              library(ontologyIndex)
    ```

    ```{r loadOBO, warning = FALSE, message = FALSE, include=FALSE}
          ontology <- get_ontology("VARPATHTOHPOOBO/hp.obo")
    ```
    ## Extract variants and annotation from gemini

    In this step the following thresholds are used:

    1. Maximum allele frequency in any population: 0.01
    2. Minimum CADD score 15 OR impact severity HIGH

    ```{r Gemini, include=FALSE}

              try(system("gemini query --header -q 'select chrom, start, end, gene, impact, impact_severity,
               cadd_scaled,polyphen_score,sift_score, clinvar_sig,max_aaf_all,sub_type,num_het, num_hom_alt, exon, codon_change, aa_change,aa_length,
              gts.VARPATIENT_ID, 
              gt_depths.VARPATIENT_ID,
              gt_quals.VARPATIENT_ID,
              gt_ref_depths.VARPATIENT_ID,
              gt_alt_depths.VARPATIENT_ID 
              from variants
               where 
              max_aaf_all < 0.01 AND 
              (cadd_scaled >= 15 OR impact_severity == \"HIGH\")
               ' VARGEMINI_DATABASEPATH --gt-filter \"gt_depths.VARPATIENT_ID > 0\" > raw_var_table_VARPATIENT_ID.tsv", intern = TRUE, ignore.stderr = TRUE))
    ```


    ```{r postGemini, include=FALSE}

              df = read_tsv("raw_var_table_VARPATIENT_ID.tsv")
              df = df %>% unite(pos, chrom, start, end)

              df = rename(df, Severity = impact_severity)

              df = rename(df, GTS = "gts.VARPATIENT_ID")
              df = rename(df, Depth = "gt_depths.VARPATIENT_ID")
              df = rename(df, CallQ = "gt_quals.VARPATIENT_ID")
              df = rename(df, RefD = "gt_ref_depths.VARPATIENT_ID")
              df = rename(df, AltD = "gt_alt_depths.VARPATIENT_ID")

              df$cadd_scaled = gsub("None",15.555, df$cadd_scaled)
              df$cadd_scaled = as.numeric(df$cadd_scaled)

              df = df %>% mutate(max_aaf_all = sprintf("%0.1e", max_aaf_all))
              df = df %>% mutate(CallQ = sprintf("%0.0f",CallQ))
              df$CallQ = as.numeric(df$CallQ)



    ```



    ```{r Export Phenotype information,echo=FALSE, include=FALSE}

              try(system("phenoparser  panel  --id VARPATIENT_ID --db VAROMIM_DATABASEPATH --out VARPATIENT_ID"))


    ```

    ## HPO terms and/or suspected diseases 

    The table contains all HPO and/or suspected diseases considered when ranking the variants. 



    ```{r make term table, echo=FALSE,include=TRUE}

              info = file.info("VARPATIENT_ID_terms.csv")
              if(info[1,1] != 0){
              terms = suppressMessages(read_csv("VARPATIENT_ID_terms.csv",col_names = FALSE))
              colnames(terms) = c("Patient ID","Term")
              terms$Description = ontology$name[terms$Term];
              terms <- terms %>% filter(!str_detect(Term, "PS"))  

              #datatable(terms,rownames=TRUE,filter = 'top',options = list(pageLength=50,autoWidth = TRUE,columnDefs = list(list(width = '5px', targets = "_all"))))
              kable(terms,"html") %>% kable_styling(bootstrap_options = "striped", full_width = F,position = "left")

              }
         ```

    ```{r Add OMIM Phenotype information to df , include=FALSE}

              info = file.info("VARPATIENT_ID_omim.csv")
              if(info[1,1] != 0){
            omim = read_csv("VARPATIENT_ID_omim.csv",col_names = FALSE)
            omim_genes = select(omim, X4,X6,X5,X2)
            colnames(omim_genes) = c("OmimPanel","gene","Inheritance","PhenotypicSeries")
            omim_genes = distinct(omim_genes, OmimPanel, gene, Inheritance,PhenotypicSeries)
              }else{
            omim_genes = tibble(OmimPanel = character(), gene = character(),Inheritance = character( ),PhenotypicSeries = character())
              }

              #df = df %>% left_join(omim_genes, by = "gene")

         ```


    ```{r Add Phenolyzer information to df , include=FALSE}

              info = file.info("VARPATIENT_ID_phenolyzer.csv")
              if(info[1,1] != 0){
            hpo = read_csv("VARPATIENT_ID_phenolyzer.csv",col_names = FALSE)
            hpo_genes = select(hpo, X4,X2,X5)
            colnames(hpo_genes) = c("PhenolyzerPanel","gene","evidence")
            hpo_genes = distinct(hpo_genes, PhenolyzerPanel, gene, evidence)

              }else{
             hpo_genes  = tibble(PhenolyzerPanel = numeric(), gene = character(), evidence = character()) 
              }

              #df = df %>% left_join(hpo_genes, by = "gene")


         ```

    ```{r joining... , include=FALSE} 
           pheno <- omim_genes %>% full_join(hpo_genes, by = "gene")

          pheno <-pheno %>%  group_by(gene) %>% dplyr::summarise(OmimPanel=paste(unique(OmimPanel), collapse="</br>"), PhenotypicSeries=paste(unique(PhenotypicSeries), collapse = "</br>"),Inheritance=paste(unique(Inheritance), collapse = "</br>"),PhenolyzerPanel=max(PhenolyzerPanel),evidence = paste(unique(evidence), collapse = "</br>")) %>% arrange(desc(PhenolyzerPanel))


           pheno <- pheno %>% replace_na(list(PhenolyzerPanel = 0.0, OmimPanel = "NA"))

           df = df %>% left_join(pheno, by = "gene")

         ```
    ## Ranking of variants

    The table is sorted by variants in genes associated with the patients disease
    phenotype and then by descending CADD score. In addition, all variants occuring
    in more than 5 patients are given a low priority. Variants with a call quality of less than 10 are discarded.

    ```{r Sorting,include=FALSE}

              #df <- df %>% group_by(OmimPanel,PhenolyzerPanel)  %>% arrange(OmimPanel)  %>%  arrange(desc(PhenolyzerPanel))  %>%  arrange(desc(cadd_scaled))

              df <- filter(df, CallQ > 10)
              df <- arrange(df,((num_het+ num_hom_alt) > 5) , desc(str_length(OmimPanel)> 1), desc(PhenolyzerPanel),desc(cadd_scaled)) #OmimPanel ,desc(PhenolyzerPanel), desc(cadd_scaled))

              df <- add_column(df , Rank = 1:nrow(df),.before = 1)

              df  <- df  %>% select(exon,codon_change,aa_change,aa_length,evidence,impact,polyphen_score,sift_score, clinvar_sig,max_aaf_all,sub_type,PhenotypicSeries,Rank, pos, GTS, gene, Severity ,cadd_scaled, PhenolyzerPanel,OmimPanel,Inheritance,CallQ, Depth,  RefD,AltD, everything())


         ```



    ```{r Write to tsv, include=FALSE}

              write.table(df, file = "VARPATIENT_ID_report_V2.tsv", append = FALSE, quote = FALSE, sep = "\t",na = "NA",row.names = FALSE,col.names = TRUE) 

    ```
    ## List of candidate variants 

    Legend: 

    ```{r Option table,echo=FALSE, include=TRUE}

         friends_data <- data_frame(
           Column = c("Rank", "pos", "GTS", "PhenolyzerPanel","OmimPanel","CallQ","Depth", "RefD","AltD","cadd_scaled"),
           Description = c("rank of variant after sorting the list according to the criteria above",
                 "genomic coordinates",
                 "genotpye",
                 "contains a score for genes returned by phenolyzer",
                 "contain genes associated with the patients disease phenotype",
                 "Phred scaled call quality", 
                 "Read depth",
                 "Depth of reference allele",
                 "Depth of alternate allele",
                 "Scaled CADD score - NOTE: for visualization purposes I give all indels a score of 15.555 and color the cell blue."
               ),                       
         )          
         kable(friends_data,"html") %>% kable_styling(bootstrap_options = "striped", full_width = F,position = "left")
    ```

    ```{r, echo=FALSE,include=TRUE}

         df <- add_column(df , Info = "<font size=\"+2\" color=\"green\">&oplus;</font>",.before = 1)





           datatable(df,rownames=FALSE,filter = 'top',escape = FALSE, 
           options = list(
           dom = 'Blfrtip',
           pageLength=25,
           autoWidth = FALSE,
           columnDefs = list(
                  list(visible=FALSE, targets= c(1:12)),
                  list(orderable = FALSE, className = 'details-control', targets = 0)
              )
           ),
           callback = JS("
           table.column(1).nodes().to$().css({cursor: 'pointer'});
           var format = function(d) {
              return '<h2>Consequence of variant on protein:</h2><table>' + 
              '<tr>' + 
              '<th>Exon</th>' + 
              '<th>Codon Change</th>' + 
              '<th>Amino Acid Change</th>' + 
              '<th>Position In Protein</th>' +
              '</tr><tr>' + 
              '<td>' + d[1] +  '</td>' + 
              '<td>' + d[2] +  '</td>' + 
              '<td>' + d[3] +  '</td>' + 
              '<td>' + d[4] +  '</td>' + 
              '</tr></table><br><hr><h2>Evidence used by phenolyzer:</h2>' + 
             '</div>' + '<div style=\"background-color:#eee; padding: .5em;\">' +
              d[5] + '</div><br><hr>' + 
              '<h2>Evidence used by OMIM:</h2><table>' + 
              '<tr>' + 
              '<th>OMIM disease</th>' + 
              '<th>Inheritance</th>' + 
              '<th>Phenotypic series</th>' + 
              '</tr><tr>' + 
              '<td>' + d[20] +  '</td>' + 
              '<td>' + d[21] +  '</td>' + 
              '<td><a href=\"https://www.omim.org/phenotypicSeries/' + d[12] + '\">' +  d[12] +  '</a></td>' +  
              '</tr></table><br><hr>' + 
              '<h2>More information:</h2><table>' + 
              '<tr>' + 
              '<th>Impact</th>' + 
              '<th>Polyphen_score</th>' + 
              '<th>Sift score</th>' + 
              '<th>Clinvar Significance</th>' +
              '<th>Max. allele frequency in any population</th>' +
              '<th>Sub Type</th>' +
              '</tr><tr>' + 
              '<td>' + d[6] +  '</td>' + 
              '<td>' + d[7] +  '</td>' + 
              '<td>' + d[8] +  '</td>' + 
              '<td>' + d[9] +  '</td>' + 
             '<td>' + d[10] +  '</td>' +
              '<td>' + d[11] +  '</td>' +
              '</tr></table>';
           };
           table.on('click', 'td.details-control', function() {
           var td = $(this), row = table.row(td.closest('tr'));
           if (row.child.isShown()) {
              row.child.hide();
              td.html('<font size=\"+2\" color=\"green\">&oplus;</font>');
           } else {
              row.child(format(row.data())).show();
              td.html('<font size=\"+2\" color=\"red\">&CircleMinus;</font>');
           }
           });"
           )) %>%            formatStyle(
           'Severity',
           backgroundColor = styleEqual(c('HIGH','MED','LOW') , c( 'lightpink', 'lightgreen', 'lightblue')
           )
           ) %>% formatStyle(
           'cadd_scaled',
           backgroundColor = styleInterval(c(15.554,15.556), c('white', 'lightblue', 'white')
           )
           ) %>% formatStyle(
           'num_het',
           backgroundColor = styleInterval(c(4.9), c('white', 'lightpink')
           )
           ) %>% formatStyle(
           'num_hom_alt',
           backgroundColor = styleInterval(c(4.9), c('white', 'lightpink')
           )
           ) 








         ```
    End.

### Create all reports {#sec-4-7-3}

A simple script to generate reports for all patients found in a
database.

The options are:

-g &lt;gemini database&gt; -o &lt;omim database&gt; -p &lt;phenolyzer
database&gt;

Example usage:

    ../scripts/create_all_variant_reports.sh -g Torrent.db -o omim.db -p hpo.txt

And here is the actual script:

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
    usage: $0  -g <gemini_database> -d <phenoparser database> -t <report template>
    EOF
          exit 1;
          }

          main() {
        export PATH=$PATH:~/gemini/anaconda/bin
        export PATH=$PATH:~/bin


          PATIENT_ID=
          GEMINI_DATABASEPATH=
          DATABASEPATH=
          TEMPLATE=


          while getopts g:d:t:  opt
          do
              case ${opt} in
              g) GEMINI_DATABASEPATH=${OPTARG};;
              d) DATABASEPATH=${OPTARG};;
              t) TEMPLATE=${OPTARG};;
              *) usage;;
              esac
          done


          if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
          if [ "${DATABASEPATH}" = "" ]; then usage; fi
          if [ "${TEMPLATE}" = "" ]; then usage; fi

          samples_array=$(gemini query -q "select name from samples" $GEMINI_DATABASEPATH)


          for i in ${samples_array[@]}; do
              step "Creating report for $i"
              try create_variant_report.sh  -i $i    -g $GEMINI_DATABASEPATH -d $DATABASEPATH  -t $TEMPLATE 
              next 
          done

          echo "DONE!!! Hurrah ";
          }

          main "$@"
