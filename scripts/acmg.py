#!/usr/bin/env python3

import os
import getopt
import requests
import requests.exceptions
import sys
import json
import pandas as pd
import multiprocessing
import collections
import math
import re
import pickle
import dill
import numpy as np
import sqlite3
from joblib import Parallel, delayed
from collections import defaultdict
from nested_dict import nested_dict
from pybiomart import Server
from sqlite3 import Error
from time import sleep
from gene_disease import restClient, geneClient, getGeneList

high_impact = ['frameshift_variant','splice_acceptor_variant','splice_donor_variant','start_lost','stop_gained','stop_lost','initiator_codon_variant','initiator_codon_variant','rare_amino_acid_variant','chromosomal_deletion']
med_impact = ['missense_variant','inframe_insertion','inframe_deletion','coding_sequence_variant','disruptive_inframe_deletion','disruptive_inframe_insertion','5_prime_UTR_truncation','3_prime_UTR_truncation','splice_region_variant','mature_miRNA_variant','regulatory_region_variant','TF_binding_site_variant','regulatory_region_ablation','regulatory_region_amplification','TFBS_ablation','TFBS_amplification']
low_impact = ['stop_retained_variant','synonymous_variant','5_prime_UTR_variant','3_prime_UTR_variant','intron_variant','coding_sequence_variant','upstream_gene_variant','downstream_gene_variant','intergenic_variant','intragenic_variant','gene_variant','transcript_variant,exon_variant','5_prime_UTR_premature_start_codon_gain_variant','start_retained_variant','conserved_intron_variant','nc_transcript_variant','NMD_transcript_variant','incomplete_terminal_codon_variant','non_coding_exon_variant','transcript_ablation','transcript_amplification','feature_elongation','feature_truncation']

domain_sources = ['TIGRFAM_domain','SMART_domains','PROSITE_profiles','Prints_domain','Pfam_domain']
ps1data = collections.namedtuple('ps1data', 'Index ps1 ps1_caution pm5 pm5_caution domains')
ppdata = collections.namedtuple('ppdata', 'Index pp2 pp2_caution pp3 pp3_caution')

regex = re.compile("\.\d+:[nc].+", re.IGNORECASE)
ppsift_min = 0.7
cadd_min = 15
consens_min = 2
af = 0.001

server = Server(host='grch37.ensembl.org',use_cache=True)
dataset_snp = (server.marts['ENSEMBL_MART_SNP'].datasets['hsapiens_snp'])
dataset_ens = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])

gene_data = []

##########
# functions
##########

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return None

def get_genes(conn):
    """
    Get the genes:disease data out of the database
    :param conn: the Connection object
    :return: dictionary of gene:impact:disease
    """
    cur = conn.cursor()
    cur.execute("select gene,impact,disease from impact_disease")
    rows = cur.fetchall()
    allrows = nested_dict()
    #for r in rows:
    #    allrows.append(r[0])
    
    for row in rows:
        if row[0] not in allrows:
            allrows[row[0]]

        if row[1] not in allrows[row[0]]['impacts']:
            allrows[row[0]]['impacts'][row[1]] = []

        if row[2] not in allrows[row[0]]['impacts'][row[1]]:
            allrows[row[0]]['impacts'][row[1]].append(row[2])
            
    allrows_dict = allrows.to_dict()        
    return allrows_dict

def getdata(row):    
    """ get the ps1 data
    :param row: tuple representing a row of a file
    :return: ps1data namedtuple
    """
    ps1 = 0
    pm5 = 0
    ps1_caution = ""
    pm5_caution = ""

    regex = re.compile("^.+:p.", re.IGNORECASE)
    regex2 = re.compile("Ter\d+$", re.IGNORECASE)
    domains = []
    domains_src = []

    if isinstance(row.vep_hgvsc, float):
        return ps1data(Index=row.Index,ps1=ps1,ps1_caution="unusable HGVSc",pm5=pm5,pm5_caution=pm5_caution,domains=list(set(domains)))
    if isinstance(row.vep_hgvsp, float):
        return ps1data(Index=row.Index,ps1=ps1,ps1_caution="unusable HGVSp",pm5=pm5,pm5_caution=pm5_caution,domains=list(set(domains)))
    if (row.vep_hgvsc == ""):
        return ps1data(Index=row.Index,ps1=ps1,ps1_caution="unusable HGVSc",pm5=pm5,pm5_caution=pm5_caution,domains=list(set(domains)))
    if (row.vep_hgvsp == ""):
        return ps1data(Index=row.Index,ps1=ps1,ps1_caution="unusable HGVSp",pm5=pm5,pm5_caution=pm5_caution,domains=list(set(domains)))
    
    ensembl_client = restClient()
    #print("{}".format(row.vep_hgvsc))
    url = 'https://grch37.rest.ensembl.org/vep/human/hgvs/'+row.vep_hgvsc
    # return type is JSON so bring the data back to an associative array
    json_ensembl = ensembl_client.perform_rest_action(
            url=url,
            params={'content-type':'application/json','domains':1},
            resource='ensembl')

    if "error" in json_ensembl:
        return ps1data(Index=row.Index,ps1=ps1,pm5=pm5,pm5_caution=pm5_caution,ps1_caution="ensembl failed with error code {} for {}".format(json_ensembl["error"],row.vep_hgvsc),domains=list(set(domains)))
    else:
        trim_hgvsp = regex.sub("", row.vep_hgvsp)
        trim_hgvsp = regex2.sub("",trim_hgvsp)
        # this goes through each object in the json
        # there may only be one
        for data in json_ensembl:
            # look for the transcript consequences as this contains the domain
            # need to know this information for pm1
            if "transcript_consequences" in data:
                for tran_con in data["transcript_consequences"]:
                    if "domains" in tran_con:
                        for dom in tran_con["domains"]:
                            if (dom["db"] in domain_sources):
                                domains.append(dom["name"])

            try:
                results1 = dataset_snp.query(
                        attributes=['clinical_significance', 'allele', 'ensembl_peptide_allele'],
                        filters={'chr_name': data['seq_region_name'],
                                 'start' : data['start'],
                                 'end' : data['end']},
                        use_attr_names = True)

                results = results1.replace(np.nan, '', regex=True)
                
                for res_row in results.itertuples():
                    # don't care about blank aa_changes in either 
                    if((res_row.ensembl_peptide_allele == "") | (row.aa_change == "")):
                        continue
                    
                    if ("pathogenic" in res_row.clinical_significance.split(',')):
                        if res_row.ensembl_peptide_allele == row.aa_change:
                            ps1 = 1
                            ps1_caution += res_row.allele + " " + row.codon_change
                        else:
                            pm5 = 1
                            pm5_caution += res_row.ensembl_peptide_allele

            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ConnectTimeout) as e:
                pm5 = 1
                if e.response is None:
                    pm5_caution += "Connection error"
                else:
                    pm5_caution += e.response.reason
             
    return ps1data(Index=row.Index,ps1=ps1,pm5=pm5,pm5_caution=pm5_caution,ps1_caution=ps1_caution,domains=list(set(domains)))

def getppdata(row):
    """ get the pp data
    :param row: tuple representing a row of a file
    :return: ppdata namedtuple
    """
    pp2 = 0
    pp2_caution = ""
    pp3 = 0
    pp3_caution = ""
    trimhgvsc = ""
    
    # assess pp2 first
    if (not isinstance(row.vep_hgvsc, float)) & (row.vep_hgvsc != ""):
        trimhgvsc = regex.sub("",row.vep_hgvsc)

        print("{}".format(trimhgvsc))

        try:
            results1 = dataset_ens.query(
                attributes = ['ensembl_gene_id','ensembl_transcript_id','variation_name', 'clinical_significance', 'synonymous_status'],
                filters = {'link_ensembl_transcript_stable_id' : trimhgvsc }, use_attr_names=True
            )

            # remove NaN
            results = results1.replace(np.nan, '', regex=True)
            sigcon = results.loc[(results.clinical_significance != "") & (results.synonymous_status == "missense_variant")]
            bcount = sigcon[(sigcon.clinical_significance == "benign") | (sigcon.clinical_significance == "likely_benign")]
            bperc = 0 if len(sigcon) == 0 else len(bcount.index) / len(sigcon)
            
            if (row.gene in gene_data):
                if ("missense_variant" in gene_data[row.gene]['impacts']):
                    if (bperc < 0.4):
                        pp2 = 1

        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ConnectTimeout) as e:
            pp2 = 0
            if e.response is None:
                pp2_caution += "Connection error"
            else:
                pp2_caution += e.response.reason
        
    # now pp3
    consens = 0;
    
    if (row.polyphen_score != "None"):
        if (float(row.polyphen_score) > ppsift_min):
            consens += 1
    if (row.sift_score != "None"):
        if (1 - float(row.sift_score) > ppsift_min):
            consens += 1
    if (row.cadd_scaled != "None"):
        if (float(row.cadd_scaled) > cadd_min):
            consens += 1
    
    if (consens >= consens_min):
        pp3 = 1
        pp3_caution = consens
    
    return ppdata(Index=row.Index,pp2=pp2,pp3=pp3,pp2_caution=pp2_caution,pp3_caution=pp3_caution)


def assignACMG(row):
    """ assign category using acmg criteria
    :param row: tuple representing a row of a file
    :return: acmg category
    """
    #Pathogenic
    #  (i) 1 Very strong (PVS1) AND
    if (row.pvs1 == 1):
    #    (a) ≥1 Strong (PS1–PS4) OR
        if (row.ps1 == 1):
            return("pathogenic")
    #    (b) ≥2 Moderate (PM1–PM6) OR
        elif ((row.pm1 + row.pm2 + row.pm4 + row.pm5) >= 2):
            return("pathogenic")
    #    (c) 1 Moderate (PM1–PM6) and 1 supporting (PP1–PP5) OR
    # note that this could be 0 for all pm but 2 for pp
    # however this is pretty much what (d) looks at anyway
        elif ((row.pm1 + row.pm2 + row.pm4 + row.pm5 + row.pp2 + row.pp3) >= 2):
            return("pathogenic")
    #    (d) ≥2 Supporting (PP1–PP5)
        elif ((row.pp2 + row.pp3) >= 2):
            return("pathogenic")

    #  (ii) ≥2 Strong (PS1–PS4) OR
    # Can't do this one as we only have 1 ps

    #  (iii) 1 Strong (PS1–PS4) AND
    elif (row.ps1 == 1):
    #    (a) ≥3 Moderate (PM1–PM6) OR
        if ((row.pm1 + row.pm2 + row.pm4 + row.pm5) >= 3):
            return("pathogenic")
    #    (b) 2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
    ### used >= 2 as it could be any one of the pms
        elif (((row.pm1 + row.pm2 + row.pm4 + row.pm5) == 2) & ((row.pp2 + row.pp3) >= 2)):
            return("pathogenic")
    #    (c) 1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)
    # Can't do this one as we only have 2 pp

    #Likely pathogenic
    #  (i) 1 Very strong (PVS1) AND 1 moderate (PM1– PM6) OR
    ### used >= 1 as it could be any one of the pms
    if ((row.pvs1 == 1) & (row.pm1 + row.pm2 + row.pm4 + row.pm5 >= 1)):
            return("likely pathogenic")
    #  (ii) 1 Strong (PS1–PS4) AND 1–2 moderate (PM1–PM6) OR
    elif ((row.ps1 == 1) & (row.pm1 + row.pm2 + row.pm4 + row.pm5 >= 1)):
            return("likely pathogenic")
    #  (iii) 1 Strong (PS1–PS4) AND ≥2 supporting (PP1–PP5) OR
    elif ((row.ps1 == 1) & (row.pp2 + row.pp3 >= 2)):
            return("likely pathogenic")
    #  (iv)  ≥3 Moderate (PM1–PM6) OR
    elif (row.pm1 + row.pm2 + row.pm4 + row.pm5 >= 3):
            return("likely pathogenic")
    #  (v) 2 Moderate (PM1–PM6) AND ≥2 supporting (PP1–PP5) OR
    ### used >= 2 as it could be any one of the pms
    elif ((row.pm1 + row.pm2 + row.pm4 + row.pm5 == 2) & (row.pp2 + row.pp3 >= 2)):
            return("likely pathogenic")
    #  (vi) 1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5) 
    # Can't do this one as we only have 2 pp

    return("undefined")

##########
# doacmg 
##########

def doacmg(df,db_file):
    """ oversee the retrieval of data and assignment of acmg category
    :param df: a pandas data frame containing the snp data
    :param db_file: sqlite database containing gene:impact:disease data
    :return: pandas data frame
    """
    # add an acmg column to store the final acmg value
    df["acmg"] = ""

    # add a domains column to complement pfam_domain
    df["domains"] = ""

    # add columns for all the acmg criteria to check
    # add default vals to col
    acmgcols = ['pvs1','ps1','pm1','pm2','pm4','pm5','pp2','pp3']
    for col in acmgcols:
        df[col] = 0
        df[col+"_caution"] = ""

    conn = create_connection(db_file)
    gene_data = get_genes(conn)

    print("got gene data",flush=True)

##########
# pvs1
##########

    for row in df.itertuples(index=True):
        pvs1 = 0
        pvs1_caution = ""
        # if clinvar_sig is pathogenic then just set to 1
        # if row.clinvar_sig contains "pathogenic":
        if "pathogenic" in row.clinvar_sig.split(','):
            pvs1 = 1
            pvs1_caution = ""
        
        # is impact (taken from the most affected transcript) in high_impact
        elif row.impact in high_impact:
            # is this impact a known mechanism of disease in the gene this variant affects
            if row.gene in gene_data:
                if row.impact in gene_data[row.gene]['impacts']:
                    pvs1 = 1;
                    if row.exon is str:
                        exons = row.exon.split('/') 
                        if exons[0] == exons[1]:
                            pvs1_caution = pvs1_caution+"variant in last exon"
        
        df.at[row.Index, 'pvs1'] = pvs1
        df.at[row.Index, 'pvs1_caution'] = pvs1_caution

    print("pvs1 set",flush=True)

##########
# ps1
##########
        
    ps1_information = Parallel(n_jobs=20)(delayed(getdata)(row) for row in df.itertuples(index=True))
    # pickle
    #filename = '/Users/rfrancis_adm/Documents/GeneticsHealth/Timo/Python/acmg-guidelines/acmg_final/ps1_information.pickle'
    #filename = '/home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/ps1_information.pickle'
    #with open(filename, 'wb') as f:
    #    dill.dump(ps1_information, f)
    #with open(filename, 'rb') as f:
    #    ps1_information = dill.load(f)

    for row in ps1_information:
        df.at[row.Index, 'ps1'] = row.ps1
        df.at[row.Index, 'pm5'] = row.pm5
        df.at[row.Index, 'ps1_caution'] = row.ps1_caution
        df.at[row.Index, 'pm5_caution'] = row.pm5_caution
        df.at[row.Index, 'domains'] = ",".join([str(x) for x in row.domains])

##########
# pm1, pm2, pm4
##########

    for row in df.itertuples(index=True):
        pm1 = pm2 = pm4 = 0
        pm1_caution = pm2_caution = pm4_caution = ""
        if (row.pfam_domain != "None") | (row.domains != ""):
            pm1 = 1
        if float(row.max_aaf_all) <= af:
            pm2 = 1
        # surely other things can cause protein length changes
        if row.impact in ['stop_lost','inframe_insertion','inframe_deletion']:
            pm4 = 1
            
        df.at[row.Index, 'pm1'] = pm1
        df.at[row.Index, 'pm1_caution'] = pm1_caution
        df.at[row.Index, 'pm2'] = pm2
        df.at[row.Index, 'pm2_caution'] = pm2_caution
        df.at[row.Index, 'pm4'] = pm4
        df.at[row.Index, 'pm4_caution'] = pm4_caution

    print("got pm1,2,4",flush=True)

##########
# pp2, pp3
##########

    pp_information = Parallel(n_jobs=20)(delayed(getppdata)(row) for row in df.itertuples(index=True))
    # pickle
    #filename = '/Users/rfrancis_adm/Documents/GeneticsHealth/Timo/Python/acmg-guidelines/acmg_final/pp_information.pickle'
    #filename = '/home/richard/pipeline/SeqNextGen_pipeline/Phenoparser/scripts/pp_information.pickle'
    #with open(filename, 'wb') as f:
    #    dill.dump(pp_information, f)
    #with open(filename, 'rb') as f:
    #    pp_information = dill.load(f)

    for row in pp_information:
        df.at[row.Index, 'pp2'] = row.pp2
        df.at[row.Index, 'pp2_caution'] = row.pp2_caution
        df.at[row.Index, 'pp3'] = row.pp3
        df.at[row.Index, 'pp3_caution'] = row.pp3_caution

    print("got pp2,3",flush=True)

##########
# assign acmg
##########

    for row in df.itertuples(index=True):
        acmg = assignACMG(row)
        df.at[row.Index, 'acmg'] = acmg

##########
# output
##########

    return(df)

##########
# MAIN
##########

def main(argv):
    #rawdatafile = "/home/richard/data/alexiadata4acmg/raw_var_table_02_CF_P17-5052911G.tsv"
    rawdatafile = "/Users/rfrancis_adm/Documents/GeneticsHealth/Timo/Python/acmg-guidelines/raw_var_table_02_CF_P17-5052911G.tsv"
    db_file = "GATK_id.db"

    rawdatafile = ''
    db_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["snpdata=","gddb="])
    except getopt.GetoptError:
        print('acmg.py -s <snpdata> -g <genediseasedb>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('acmg.py -s <inputfile> -g <genediseasedb>')
            sys.exit()
        elif opt in ("-s", "--snpdata"):
            rawdatafile = arg
        elif opt in ("-g", "--gddb"):
            db_file = arg

    if ((rawdatafile == "") | ()):
        print("You must provide both an input file and a dabase file\nacmg.py -s <inputfile> -g <genediseasedb>")
        sys.exit()
    else:
        df1 = pd.read_csv(rawdatafile, sep='\t')
        df = df1.replace(np.nan, '', regex=True)
        doacmg(df,db_file)


if __name__ == "__main__":
   main(sys.argv[1:])

