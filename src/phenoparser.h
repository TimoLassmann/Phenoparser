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
        char* mimNumber;//	601282
        char* phenotype;//	Epidermolysis bullosa simplex with muscular dystrophy
        char* phenotypeMimNumber;//	226670
        char* phenotypeMappingKey;//	3
        char* phenotypeInheritance;//	Autosomal recessive
        char* sequenceID;//	7256
        char* chromosome;//	8
        char* chromosomeSymbol;//	8
        char* chromosomeSort;//	552
        char* chromosomeLocationStart;//	143915146
        char* chromosomeLocationEnd;//	143976799
        char* transcript;//	uc064rfy.1
        char* cytoLocation;//	8q24
        char* computedCytoLocation;//	8q24.3
        char** geneSymbols;//	PLEC1, PLEC, PLTN, EBS1, LGMD2Q, EBSOG, EBSPA, EBSMD, EBSND
        char* geneInheritance;//	(null)      
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
