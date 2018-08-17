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
|  Decompose                 |  [vt](<https://genome.sph.umich.edu/wiki/Vt>)                                         |  Decompose multiallelic variants
|  Normalize                 |  vt                                                                                   |  Standardize vcf format (see [Tan paper](https://academic.oup.com/bioinformatics/article/31/13/2202/196142) )
|  =====                     |                                                                                       |
|  **Phenotype data**        |                                                                                       |  Parse phenotype data
|  Create gene panels        |  [phenoparser](<https://github.com/TimoLassmann/Phenoparser>)                         |  A simple C program to query omim and run phenolyzer, storing all results in a database
|                            |  requires: [phenolyzer](<http://phenolyzer.wglab.org/>)                               |  Creates gene lists (panels) from input HPO terms
|                            |  requires: [sqlite](<https://www.sqlite.org/download.html>)                           |  A single database is used to store all data 
|  =====                     |                                                                                       |
|  **Variant data**          |                                                                                       |  Annotate and store/manage variant data
|  Annotation                |  [Variant Effect Predictor](<https://www.ensembl.org/info/docs/tools/vep/index.html>) |  Add annotation prior to import into Gemini
|  Storage/management        |  [GEMINI](<https://gemini.readthedocs.io/en/latest/>)                                 |  Database holding variants and annotation; can be queried using sql statements
|                            |  requires: [sqlite](<https://www.sqlite.org/download.html>)                           |  We generate a separate database for each sequencing technology (IonTorrent, Illumina and Solid)
|                            |  requires: perl                                                                       |
|  =====                     |                                                                                       |
|  **Reporting**             |                                                                                       |  
|  Generate patient reports  |  Shell script and [R](<https://www.r-project.org/>)                                   |  Extracts variant information and combines it with phenotype data into one report of ranked candidate causal variants
|                            |  requires: [Human Phenotype Ontology](<https://github.com/obophenotype/human-phenotype-ontology>) hp.obo file                           |  

In addition to the scripts in this distribution (detailed below), the highlighted software components in the table above need to be installed for the pipeline to work.


Notes for particular software
-----------------------------

**Variant Effect Predictor**

VEP installs are linked to Ensembl release versions. So that HGVS expressions are properly created in the reports, be sure to download a compatible (i.e. same Ensembl release) genomic reference sequence file.
Include the path to this file in the sng_config file using the SNGVEPREF variable.

For example:<br/>
Using VEP [version 84](<https://github.com/Ensembl/ensembl-tools/archive/release/84.zip>)<br/>
Use Assembly [version 84](<https://bit.ly/2uQa9oB>)<br/>

In sng_config:<br/>
    `SNGVEPBIN=/path/to/installed/vep`<br/>
    `SNGVEPREF=/path/to/bgzipped/Homo_sapiens.GRCh37.dna.primary_assembly.fa`<br/>

Also, if the pipeline is to be used by multiple people on a given server then the VEP cache directory should be accessible to all. By default this directory is created in the installing user's `$HOME/.vep`. However, this will prevent others from accessing it. Instead, use the `--dir` option as [here](<https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#cacheopt>) to specify a communal directory for this data and refer to this in your config file using `SNGVEPTMP`.<br/>

In sng_config:<br/>
    `SNGVEPTMP=/path/to/installed/vep/cache`<br/>

**GEMINI**

GEMINI solely supports human genetic variation mapped to build 37 (aka hg19) of the human genome so be sure to obtain GRCh37 Ensembl data.

**R**

R is used to perform data manipulation and generate the reports. In order for the template to run properly the following R packages need to be installed:

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
The pipeline can handle individual patient VCF files or a single pre-merged VCF file and accompanying phenotype data.

Individual VCF
--------------

For each patient, create one directory named after the patient ID (preferably the one used in the generation of the VCF file).
In this directory, create two sub-directories: vcf and pheno. The latter contains text files with the phenotype / disease annotation.

    /home/user/patient_data/D12-3456/
    ├── pheno
    │   ├── hpo.txt
    │   └── omim.txt
    └── vcf
        └── D12-3456.vcf

Merged VCF
----------

For a merged VCF file the structure is slightly different. Still create two sub-directories as above: vcf and pheno, then place your merged VCF file in the vcf folder.
For each patient create [patientID]_hpo.txt and [patientID]_omim.txt files in the pheno folder where this data is available.
The key is to make sure that the [patientID] in these file names is identical to the corresponding sample names used in the VCF file (use: `bcftools query -l file.vcf` to obtain the sample names in your VCF).

    /home/user/patient_data/
    ├── pheno
    │   ├── D12-3456_hpo.txt
    │   └── D12-3456_omim.txt
    │   ├── D98-7654_hpo.txt
    │   └── D98-7654_omim.txt
    └── vcf
        └── merged.vcf

The analysis of this data should be performed in a separate directory. the first part of the pipeline copies all raw data to this directory to leave the original data untouched.

In this document this analysis directory is referred to as the "patient_analysis" directory.

Code
====

The pipeline consists of various scripts that process the input data and run the software components above.

SNG Config File
---------------

The code ships with a sample config file. This is essentially a list of environment variables denoting the location of the relevant software components and data containing or temporary directories

```bash
# this is the directory to which all patient data is copied and where the pipeline writes temporary files
SNGTMP="tmp"
# this is the directory where the gemini binary lives
SNGGEMINIBIN="....gemini/tools/bin/"
# a temp directory that is made in your patient_analysis folder to store temporary Gemini files
SNGGEMINI_TMP="geminitmp"
# this is the subdirectory of the phenoparser code that contains the pipeline scripts
SNGSCRIPTS="..../Phenoparser/scripts/"
# file path to the bcftools root directory where the bcftools binary lives
SNGBCFTOOLS="/usr/local/src/bcftools-1.6/"
# file path to the htslib root directory where bgzip and tabix binaries live
SNGHTSLIB="/usr/local/src/htslib-1.6/"
# file path to the grabix root directory where the grabix binary lives
SNGGRABIX="/usr/local/src/grabix/grabix-master/"
# this is the directory where the vt binary lives
SNGVTBIN="/usr/local/src/vt/"
# this is the directory where the vt reference sequence lives - should be the reference the variants were aligned to
SNGVTREF="/usr/local/src/vt/hg19/hg19.fa.gz"
# this is the directory where the VEP binary lives
SNGVEPBIN="/usr/local/src/vep/"
# this is the directory where VEP has installed its cache and plugins - see note above
SNGVEPTMP="/home/vep_cache"
# this is the directory where the bgzipped genome assembly file for VEP lives
SNGVEPREF="/path/to/bgzipped/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
# this is the phenolyzer root directory where the disease_annotation.pl perl script lives 
SNGPLBIN="/usr/local/src/Phenolyzer/phenolyzer/"
# this is the directory where the phenoparser binary lives
SNGPPBIN="..../Phenoparser/src/"
# this is the location of the human phenotype ontology file
SNGHPOBO="/path/to/hp.obo"
```

The above locations are made available in your environment with some added to your `$PATH` by sourcing this file at the start of the pipeline

Common functions
----------------

The script `common.sh` contains a number of functions used across the pipeline including running vt and indexing
vcf files along with some utility functions. Each of the other scripts in this collection source this file.

Pre-processing
--------------

The script `run_pre_merge.sh` initially copies all variant and phenotype to your analysis directory so as to ensure that the original files
remain untouched. It then normalizes and indexes each vcf file, extracting the variant caller used.
It should be run from within the directory you create for the analysis (e.g. 20180108_patient_analysis).

It takes two options (run without option will print a help message):

```
-i <directory holding data of multiple patients>
-l <log file to keep track of what’s happening>
```

Example usage:

```
run_pre_merge.sh -i ~/patient_data -l pre_merge_log
```

Build gemini database(s)
------------------------

This script `build_gemini.sh` combines vcf files from different variant callers and loads them into GEMINI.
vt decompose and normalise is then run on these files for good measure.

The script will search (grep!) through the `sample_info.txt` file created by `run_pre_merge.sh` to identify all vcf files from a given caller.
For example we used:

-   “Life” for life technologies - SOLID
-   “Torrent” for Ion Torrent data
-   “GATK” for Illumia files

i.e. the base caller is used to identify the technology as older vcf files are devoid of usable meta-data.

There are three parameters:

```
-p <"technology" like the options above>
-d <output GEMINI database name>
-l <log file>
```

Example usage:

```
build_gemini_db.sh -p GATK -d GATK.db -l gatk.log 
build_gemini_db.sh -p Life -d Life.db -l gatk.log 
build_gemini_db.sh -p Torrent -d Torrent.db -l gatk.log
```

Make Omim database
------------------

The script `make_omim_database.sh` runs [phenoparser]( <https://github.com/TimoLassmann/Phenoparser> ) to extract information from the OMIM database based on any supplied OMIM phenotypic terms for each patient.
It also pulls down Phenotypic Series records - collections of entries with overlapping clinical manifestations. In all cases, disease associated genes and annotation are retrieved and stored in an SQLite database
to make sure the analysis is reproducible even if OMIM changes. This also allows you to re-run previously un-diagnosed cases if there is a major update.

To access OMIM you need an OMIM key which you can request online at ( <https://www.omim.org/api> ). Applications can take a few days to be approved.

The options to the script are:

```
-d <output database where all results will be stored>
-l <log file to keep track of what's happening>
-k <OMIM key>
```

Example usage:

```
make_omim_database.sh -d omim.db -l omim.log -k <KEY>
```

Phenolyzer
----------

The script `make_extended_phenolyzer_database.sh` takes in any available HPO terms **and Disease terms** for each patient and
queries Phenolyzer. All data is stored in the same SQLite database as the omim data.

The options are:

```
-d <output database where all data will be stored - needs to be the same database name as above>
-l <log file to keep track of what's happening>
```

Example usage:

```
make_extended_phenolyzer_database.sh -d omim.db -l phenolyzer.log
```

Patient reports
---------------

### Script

The pipeline uses RMarkown and knitr to create per-patient reports. A template contains special variables that are replaced by patient details.

Either a single report can be generated for a list of patient IDs or reports for all patients can be created.

The script `create_variant_reports.sh` brings all the data together into one report. It extracts variants together with their annotation from gemini and overlays in
silico gene panels from OMIM and Phenolyzer.  The script works by copying a report template, replacing placeholder variables with patient data and then running the template in R. The
output is a html file containing information about the analysis as well as a tab separated file with the variant table, which is ranked in order of importance.

The options are:

```
-i <patient id - comma separated list or absent to produce reports for all patients>
-g <gemini database path for a particular platform>
-d <phenotype database path>
-t <report template>
-l <logfile to keep track of what is happening>
```

Example usage:

```
# specific patients
create_variant_reports.sh -i "patient1,patient2" -g Torrent.db -d omim.db -t report_master_template.Rmd -l reports.log

# all patients
create_variant_reports.sh -g Torrent.db -d omim.db -t report_master_template.Rmd -l reports.log
```

### Template

The template `report_V2_template.Rmd` combines variants and disease annotation from the various sources generated in the pipeline.
It then orders the subsequent interactive table so as to place the most likely causative variants at the top of the list.
Additional annotation is presented such as HGVS expressions, links out to UniProt and ClinVar plus other useful references.


Complete script
---------------

A complete script to run the entire pipeline can be created, which might look something like the following.

```bash
#!/usr/bin/env bash

# source the config file
. ./sng_config

# make a clean copy of the data and merge
run_pre_merge.sh -i ~/data/patient_data -l pre_merge_log

# create gemini databases for the three platforms
build_gemini.sh -p GATK -d GATK.db -l gatk.log
build_gemini.sh -p Torrent -d Torrent.db -l torrent.log
build_gemini.sh -p Life -d Life.db -l life.log

# make an omim database - uses phenoparser
make_omim_database.sh -d omim.db -l omim.log -k <OMIM_KEY>

# add phenolyzer data to omim database - uses phenoparser
make_extended_phenolyzer_database.sh -d omim.db -l phenolyzer.log

# create a patient report
##create_variant_reports.sh -i "patient1,patient2" -g Torrent.db -d omim.db -t $SNGSCRIPTS/report_V2_template.Rmd -l reports.log

# or all patient reports
create_variant_reports.sh -g Torrent.db -d omim.db -t $SNGSCRIPTS/report_V2_template.Rmd -l reports.log
create_variant_reports.sh -g GATK.db -d omim.db -t $SNGSCRIPTS/report_V2_template.Rmd -l reports.log
create_variant_reports.sh -g Life.db -d omim.db -t $SNGSCRIPTS/report_V2_template.Rmd -l reports.log
```
