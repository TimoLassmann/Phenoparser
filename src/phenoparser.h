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

#include "tldevel.h"

struct parameters{
        char* local_sqlite_database_name;
        char* omimkey;
        char* patient_id;
        char* phenofile;
        
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

/* user interface */
int print_global_help(int argc, char * argv[]);
struct parameters* get_panel_param(int argc, char * argv[]);


struct parameters* get_insert_param(int argc, char * argv[]);

void free_param(struct parameters* param);

#endif
