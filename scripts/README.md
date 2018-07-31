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
|                            |  requires: perl                                                                       |
|  =====                     |                                                                                       |
|  **Reporting**             |                                                                                       |  
|  Generate patient reports  |  Shell script and R ( <https://www.r-project.org/> )                                  |  Extracts variant information and combines it with phenotype data into one report of ranked candidate causal variants
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
For each patient, create one directory named after the patient ID (preferably the one used in the generation of the VCF file).
In this directory, create two sub-directories: vcf and pheno. The latter contains text files with the phenotype / disease annotation.

    /home/user/patient_data/D12-3456/
    ├── pheno
    │   ├── hpo.txt
    │   └── omim.txt
    └── vcf
        └── D12-3456.vcf

The analysis of this data should be performed in a separate directory. the first part of the pipeline copies all raw data to this directory to leave the original data untouched.

In this document this analysis directory is referred to as the "patient_analysis" directory.

Code
====

The pipeline consists of various scripts that process the input data and run the software components above.

SNG Config File
---------------

The code ships with a sample config file. This is essentially a list of environment variables denoting the location of the relevant software components

```bash
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
# this is the directory where VEP can store temporary files
SNGVEPTMP="veptmp"
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

### Master Template {#sec-4-7-2}

### Create all reports {#sec-4-7-3}

A simple script to generate reports for all patients found in a
database.

The options are:

-g &lt;gemini database&gt; -o &lt;omim database&gt; -p &lt;phenolyzer
database&gt;

Example usage:

    ../scripts/create_all_variant_reports.sh -g Torrent.db -o omim.db -p hpo.txt

Complete script
---------------

A complete script to run the entire pipeline can be created, which might look something like the following.

```bash
#!/usr/bin/env bash

# source the config file
. /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/sng_config
# make a clean copy of the data and merge
#run_pre_merge.sh -i /home/richard/data/patient_data -l pre_merge_log
# create gemini databases for the three platforms
build_gemini.sh -p GATK -d GATK.db -l gatk.log
build_gemini.sh -p Torrent -d Torrent.db -l torrent.log
build_gemini.sh -p Life -d Life.db -l life.log
# make an omim database - uses phenoparser
make_omim_database.sh -i tmp/ -d omim.db -l omim.log -k Wqy5lssmS7uWGdpyy8H9zw
# add phenolyzer data to omim database - uses phenoparser
make_extended_phenolyzer_database.sh -i tmp -d omim.db -l phenolyzer.log
# create a patient report
#create_variant_report.sh -i 07_D14-0866 -g Torrent.db -d omim.db -t /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/report_V2_template.Rmd -o /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/hp.obo
# or all patient reports
create_all_variant_reports.sh -g Torrent.db -d omim.db -t /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/report_V2_template.Rmd -o /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/hp.obo
create_all_variant_reports.sh -g GATK.db -d omim.db -t /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/report_V2_template.Rmd -o /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/hp.obo
create_all_variant_reports.sh -g Life.db -d omim.db -t /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/report_V2_template.Rmd -o /home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/hp.obo
```
