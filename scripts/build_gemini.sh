#!/usr/bin/env bash

  pwd=$(pwd)

  usage(){
cat <<EOF
usage: $0  -p <platform> -d <database> -l <logfile>
EOF
      exit 1;
  }

main(){
 
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
 
      # converting to an array
      samplelist=( $(cat sample_info.txt | grep $PLATFORM | cut -f 1) );
 
      step "Merge all vcf files"
      # RF checking to see if there's only one vcf file for this technology
      if [ ${#samplelist[@]} == 1 ]; then
          echo " - Only one sample";
          cat "${samplelist[@]}" > $pwd/tmp/combined_$PLATFORM.vcf;
      else
          echo " - Merging samples";
          try bcftools merge "${samplelist[@]}" > $pwd/tmp/combined_$PLATFORM.vcf
      fi
      next
 
      step "Run VT and VEP on combined"
      pipeline $pwd/tmp/combined_$PLATFORM.vcf;
      next
 
      step "Load into gemini";
      try load_into_gemini  $pwd/tmp combined_$PLATFORM.d.n.vep.vcf.gz
      next 

      step "Extract from gemini";
      try extract_from_gemini
      next
 
      echo "Done!"
  }

pipeline(){
if ! [ "$1" ]
then
echo "Pipeline function needs an input vcf file";
return 1;
fi

local OUTNAME=
local INNAME=

local WORKINGDIR=$pwd/tmp

# run vt_pipeline
vt_pipeline $1 INNAME $WORKINGDIR;

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

load_into_gemini () {
if ! [ "$1" ]
then
echo "index needs and input file";
return 1;
fi

$GEMINIBIN/gemini load --passonly -v $1/$2 -t VEP --tempdir $GEMINI_TMP --cores $NUM_THREADS $GEMINI_DATABASE
return 0;
}

extract_from_gemini () {

    $GEMINIBIN/gemini query --header --show-samples --sample-delim ";" -q "select variant_id, num_het, num_hom_alt from variants" $GEMINI_DATABASE > sample_variants_$PLATFORM.txt
    return 0;
}

function run_vep () {
if ! [ "$2" ]
then
echo "run_vep needs a working directory";
return 1;
fi

if ! [ "$3" ]
then
echo "run_vep needs input vcf file";
return 1;
fi

local  __resultname=$1
local myresultname=$(basename "$3" | cut -d. -f1).d.n.vep.vcf

$DOCKERCMD/docker run -v $2:/data seqnextgen_vep perl /src/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl -i /data/$3 -o /data/$myresultname --vcf --fork $NUM_THREADS --offline --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --assembly GRCh37 --dir_cache /root/.vep

eval $__resultname="'$myresultname'"
}

main "$@";
