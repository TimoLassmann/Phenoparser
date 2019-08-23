#ifndef PHENOPARSER_HEADER

#define PHENOPARSER_HEADER

#define MAX_ALT_GENE_NAMES 10

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

#include <getopt.h>
#include <stdint.h>
#include <float.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <expat.h>
#include <curl/curl.h>
#include <sqlite3.h>


#include "openssl/conf.h"
#include "openssl/err.h"
#include "openssl/engine.h"
#include "openssl/ssl.h"


#include "tldevel.h"

#include "rbtree.h"


struct parameters{
        char* local_sqlite_database_name;
        char* omimkey;
        char* patient_id;
        char* phenofile;

        char* outfile;
        char** infile;
        int num_infiles;
};

struct OMIM_list{
        struct OMIM** terms; 
        int num_entries;
        int num_malloced;
};

struct OMIM{
        char* mimNumber;
        char* phenotype;
        char* phenotypeMimNumber;
        char* phenotypeMappingKey;
        char* phenotypeInheritance;
        char* phenotypicSeriesNumber;
        char* sequenceID;
        char* chromosome;
        char* chromosomeSymbol;
        char* chromosomeSort;
        char* chromosomeLocationStart;
        char* chromosomeLocationEnd;
        char* transcript;
        char* cytoLocation;
        char* computedCytoLocation;
        char** geneSymbols;
        char* geneInheritance;
};

struct MemoryStruct {
        char *memory;
        size_t size;
};
 
struct ParserStruct {
        int ok;
        size_t tags;
        size_t depth;
        struct MemoryStruct characters;
        struct OMIM_list* ol;
        struct OMIM** omimout;
        int num_omim_out; 

};

struct string_struct{
        char* name;
};



/* user interface */
int print_global_help(int argc, char * argv[]);
struct parameters* get_panel_param(int argc, char * argv[]);

struct parameters* get_termlist_param(int argc, char * argv[]);
struct parameters* get_insert_param(int argc, char * argv[]);
struct parameters* get_readphe_param(int argc, char * argv[]);

void free_param(struct parameters* param);

/* parser.c  */


int phenotype_series_search(struct parameters* param, char* search_term);
int search_and_insert_disease(struct parameters* param, char* search_term, struct rbtree_root* series);

/* database stuff */
char* create_query_string(char* query,int* query_len,const char * format, ...);


int check_if_db_exists_otherwise_create(struct parameters* param);

/* kind of database but defined in main.c */
int action_insert_into_sqlite( struct OMIM_list* ol,struct parameters* param, char* search_term, struct rbtree_root* series);
int enter_term(struct OMIM_list* ol, char* name,char* value );


/* omim list */

struct OMIM_list* init_omim_list(int n);
int clear_omim_list(struct OMIM_list* ol);
int resize_omim_list(struct OMIM_list* ol, int add);
void free_omim_list(struct OMIM_list* ol);

struct OMIM* init_omim_entry(void);

int clear_term(struct OMIM* omim);
void free_omim(struct OMIM* omim);

/* misc */
struct rbtree_root* make_string_tree(void);
int remove_comma_quote(char* in);

#endif
