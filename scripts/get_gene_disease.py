#!/usr/bin/env python3

from gene_disease import restClient, geneClient, getGeneList
import multiprocessing as mp
from nested_dict import nested_dict
import sqlite3
from sqlite3 import Error
import collections
import dill

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

def make_table(conn):
    """
    Make the table to store the impact:disease data
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute('''create table if not exists impact_disease( gene TEXT default NULL, impact TEXT default NULL, disease text default NULL )''')
    cur.execute('''delete from impact_disease''')
 
    return None

def get_genes(conn):
    """
    Get the genes out of the gemini database
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("select distinct gene from gene_summary where is_hgnc = 1 and hgnc_id is not null")
    rows = cur.fetchall()
    allrows = []
    for r in rows:
        allrows.append(r[0])
    return allrows

def load_table(conn,row_data):
    """
    Load the retrieved data
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()

    for rowd in row_data:
        for row in rowd:
            if row.error == "":
                cur.execute('''insert into impact_disease(gene,impact,disease) VALUES (?,?,?)''', (row.gene,row.impact,row.disease))
            else:
                print("error from {} => {}".format(row.gene,row.error))

    conn.commit()

def getginfoasync(geneList,omim_key):
    g_information = []
    #pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(100)
    #result_objects = [pool.apply_async(getGeneList, args=(gene, omim_key)) for gene in geneList]
    g_information = pool.starmap_async(getGeneList, [(gene, omim_key) for i, gene in enumerate(geneList)]).get()
    #print(results)
    # result_objects is a list of pool.ApplyResult objects
#    g_information = [r.get()[0] for r in result_objects]
    pool.close()
    #pool.join()

    return g_information

def getginfo(geneList,omim_key):
    g_information = []
    pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(15)
    g_information = [pool.apply(getGeneList, args=(gene,omim_key)) for gene in geneList]
    #g_information = pool.starmap(getGeneList, [(gene, omim_key) for gene in geneList])
    #g_information = Parallel(n_jobs=100)(delayed(getGeneList)(gene,omim_key) for gene in geneList)
    pool.close()

    return g_information

if __name__ == '__main__':
    gene_db = "GATK.db"
    impact_disease_db = "GATK_id.db"
    omim_key = "Wqy5lssmS7uWGdpyy8H9zw"
    geneList = []
     
    # create database connections
    print("connecting to {}".format(gene_db),flush=True)
    gene_conn = create_connection(gene_db)
    print("connecting to {}".format(impact_disease_db),flush=True)
    id_conn = create_connection(impact_disease_db)
    
    # make the table
    print("making impact_disease database",flush=True)
    make_table(id_conn)

    # get the genes from the database
    print("getting genes",flush=True)
    #geneList = ['IRF6','MT-ND4','VWA1']
    geneList = get_genes(gene_conn)

    # get the information
    print("getting gene information",flush=True)
    #g_information = getginfo(geneList,omim_key)
    g_information = getginfoasync(geneList,omim_key)
    # pickle it
    filename = '/home/richard/pipeline/scripts/Phenoparser/scripts/g_information.pickle'
    with open(filename, 'wb') as f:
            dill.dump(g_information, f)
    print("collected information for {} genes".format(len(g_information)),flush=True)

    # load the data
    print("loading gene information",flush=True)
    load_table(id_conn,g_information)

    print("closing connections",flush=True)
    gene_conn.close()
    id_conn.close()


