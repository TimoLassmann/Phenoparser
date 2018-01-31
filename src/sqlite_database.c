
#include "phenoparser.h"

                               



int check_if_db_exists_otherwise_create(struct parameters* param)
{
        char buffer[BUFFER_LEN]; 
        sqlite3 *sqlite_db = NULL;
        int rc;          
        ASSERT(param->local_sqlite_database_name != NULL,"No database.");
        if(!my_file_exists(param->local_sqlite_database_name)){
                DPRINTF1("Creating sqlite database:%s\n", param->local_sqlite_database_name);
                	
                rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
                if(rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                DPRINTF2("Creating table: patient\n");
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE patient(patient_id TEXT NOT NULL, DiseaseSearch TEXT,  unique (patient_id,DiseaseSearch) );");
                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                DPRINTF2("Creating table: diseaseMIM\n");
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE diseaseMIM(DiseaseSearch TEXT, phenotypeMimNumber INT NOT NULL,phenotypeDescription TEXT,phenotypeInheritance TEXT, unique (DiseaseSearch,phenotypeMimNumber, phenotypeInheritance));");
                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                DPRINTF2("Creating table: MIMgene\n");
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE MIMgene(phenotypeMimNumber INT NOT NULL, gene TEXT NOT NULL,unique (phenotypeMimNumber,gene));");
                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                DPRINTF2("Creating table: Phenolyzer\n");
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE phenolyzer(patient_id TEXT NOT NULL,gene TEXT NOT NULL,identifier INT,score REAL,status TEXT NOT NULL, unique (patient_id,gene));"); 

                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                DPRINTF2("Creating table: Phenolyzer Gene Information \n");
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE phenolyzerGeneData(patient_id TEXT NOT NULL,gene TEXT NOT NULL,gene_evidence TEXT NOT NULL, unique (patient_id,gene));"); 

                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

                
                rc = sqlite3_close(sqlite_db);
                if( rc!=SQLITE_OK ){
                        ERROR_MSG("sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));
                }
        }
        return OK;
ERROR:
        if(sqlite_db){
                DPRINTF3("Error in check if db exists handled...");
                rc = sqlite3_close(sqlite_db);
                if( rc!=SQLITE_OK ){
                        fprintf(stderr,"sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));
                }
        }  
        return FAIL;
}

/*
  Some test queries 
  INSERT INTO patient VALUES ('OTTO', 'epidermolysis bullosa','');

  INSERT INTO diseaseMIM VALUES ('epidermolysis bullosa',603576);
  INSERT INTO diseaseMIM VALUES ('epidermolysis bullosa',613815);

  INSERT INTO MIMgene VALUES (603576,'TRPM1');
  INSERT INTO MIMgene VALUES (603576,'MLSN1');
  INSERT INTO MIMgene VALUES (603576,'CSNB1C');

  INSERT INTO MIMgene VALUES (613815,'CYP21A2');
  INSERT INTO MIMgene VALUES (613815,'CYP21');

  INSERT OR IGNORE INTO MIMgene VALUES (613815,'CA21H');


  SELECT * FROM MIMgene  AS a INNER JOIN diseaseMIM AS b ON a.phenotypeMimNumber = b.phenotypeMimNumber;


  SELECT
  patient_id AS ID,
  patient.DiseaseSearch AS DISEASE,
  MIMgene.phenotypeMimNumber AS MIMnumber,
  diseaseMIM.phenotypeDescription AS DESC,
  MIMgene.gene AS GENE
  FROM
  patient 
  INNER JOIN diseaseMIM ON diseaseMIM.DiseaseSearch = patient.DiseaseSearch 
  INNER JOIN MIMgene ON MIMgene.phenotypeMimNumber = diseaseMIM.phenotypeMimNumber 
  WHERE patient_id == "dOTTO";

  SELECT 
  *
  FROM
  diseaseMIM
  INNER JOIN MIMgene ON MIMgene.phenotypeMimNumber = diseaseMIM.phenotypeMimNumber;

*/
