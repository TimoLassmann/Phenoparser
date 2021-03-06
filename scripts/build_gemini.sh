#!/usr/bin/env bash

# source common
. $SNGSCRIPTS/common.sh

  pwd=$(pwd)

  usage(){
cat <<EOF
usage: $0  -p <platform> -d <database> -l <logfile> -t <threads>
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
      if [ "${LOG_STEPS}" = "" ]; then usage; fi
 
      # converting to an array
      # the use of sort -u allows pre-merged vcf files to be used
      samplelist=( $(cat sample_info.txt | grep $PLATFORM | cut -f 1 | sort -u) );
 
      step "Merge all vcf files" > $LOG_STEPS 2>&1;
      # RF checking to see if there's only one vcf file for this technology
      if [ ${#samplelist[@]} == 0 ]; then
          echo " - No samples detected for this platform ($PLATFORM). Exiting!";
          try exit 1;
      elif [ ${#samplelist[@]} == 1 ]; then
          echo " - Only one sample file";
          try cat "${samplelist[@]}" > $pwd/$SNGTMP/combined_$PLATFORM.vcf;
      else
          echo " - Merging samples";
          try $SNGBCFTOOLS/bcftools merge "${samplelist[@]}" > $pwd/$SNGTMP/combined_$PLATFORM.vcf
      fi >> $LOG_STEPS 2>&1;
      next >> $LOG_STEPS 2>&1;
 
      step "Run VT and VEP on combined" >> $LOG_STEPS 2>&1;
      pipeline $pwd/$SNGTMP/combined_$PLATFORM.vcf >> $LOG_STEPS 2>&1;
      next >> $LOG_STEPS 2>&1;
 
      step "Load into gemini" >> $LOG_STEPS 2>&1;
      try load_into_gemini  $pwd/$SNGTMP combined_$PLATFORM.d.n.vep.vcf.gz >> $LOG_STEPS 2>&1;
      next >> $LOG_STEPS 2>&1;

      echo "$0 Done!"
  }

pipeline(){
if ! [ "$1" ]
then
echo "Pipeline function needs an input vcf file";
return 1;
fi

local OUTNAME=
local INNAME=

local WORKINGDIR=$pwd/$SNGTMP

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
#step "Docker cleanup"
#try cleanup_docker
#next 

return 0;
}

load_into_gemini () {
if ! [ "$1" ]
then
echo "gemini needs a working directory";
return 1;
fi

if ! [ "$2" ]
then
echo "gemini needs and input file";
return 1;
fi

$SNGGEMINIBIN/gemini load --passonly -v $1/$2 -t VEP --tempdir ./$SNGGEMINI_TMP --cores $NUM_THREADS $GEMINI_DATABASE
return 0;
}

run_vep () {
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

#$SNGDOCKERCMD run -v $2:/data seqnextgen_vep perl /src/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl -i /data/$3 -o /data/$myresultname --vcf --fork $NUM_THREADS --offline --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --assembly GRCh37 --dir_cache /root/.vep
#$SNGDOCKERCMD run -v $2:/data seqnextgen_vep2 perl /src/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl -i /data/$3 -o /data/$myresultname --vcf --fork $NUM_THREADS --offline --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --fasta /src/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,HGVSc,HGVSp,HGVS_OFFSET --assembly GRCh37 --dir_cache /root/.vep
if [ ! -d "$SNGVEPTMP" ]; then
 echo "VEP cache directory $SNGVEPTMP does not exist. See documentation.";
fi
perl $SNGVEPBIN/vep -i $2/$3 -o $2/$myresultname --vcf --fork $NUM_THREADS --offline --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --fasta $SNGVEPREF --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,HGVSc,HGVSp,HGVS_OFFSET --assembly GRCh37 --dir_cache $SNGVEPTMP

eval $__resultname="'$myresultname'"
}

main "$@";
