

#include "phenoparser.h"

  
struct OMIM_list* init_omim_list(int n)
{
        struct OMIM_list* ol = NULL;
        int i;
        ASSERT(n != 0,"Requesting zero entries!"); 
        
        MMALLOC(ol,sizeof(struct OMIM_list));
        ol->terms = NULL;
        ol->num_entries = 0;
        ol->num_malloced = n;
        MMALLOC(ol->terms,sizeof(struct OMIM*)* n);
        for(i = 0; i < n;i++){
                ol->terms[i] = NULL;
                RUNP(ol->terms[i] =  init_omim_entry());
        }
        
        
        return ol;
ERROR:
        return NULL;
}

int clear_omim_list(struct OMIM_list* ol)
{
        int i;
        ASSERT(ol != NULL,"No list!");
        ol->num_entries = 0;
        for(i = 0;i < ol->num_malloced;i++){
                RUN(clear_term(ol->terms[i]));
        }
        return OK;
ERROR:
        return FAIL;

}

int resize_omim_list(struct OMIM_list* ol, int add)
{
        int i;
        ASSERT(ol != NULL,"Nolist...");
        
        MREALLOC(ol->terms,sizeof(struct OMIM*)* (ol->num_malloced+add));
        
        for(i = ol->num_malloced; i < ol->num_malloced+add;i++){
                ol->terms[i] = NULL;
                RUNP(ol->terms[i] =  init_omim_entry());
        }
        ol->num_malloced += add;
        
        return OK;
ERROR:
        return FAIL;
}


void free_omim_list(struct OMIM_list* ol)
{
        if(ol){
                int i;
                for(i = 0; i < ol->num_malloced;i++){
                        free_omim(ol->terms[i]);
                }
                MFREE(ol->terms);
                MFREE(ol);
        }
}


        
struct OMIM* init_omim_entry(void)
{
        struct OMIM* omim = NULL;
        
        int i;
        MMALLOC(omim,sizeof(struct OMIM));
        omim->chromosome = NULL;
        omim->chromosomeLocationEnd = NULL;
        omim->chromosomeLocationStart = NULL;
        omim->chromosomeSort = NULL;
        omim->chromosomeSymbol = NULL;
        omim->computedCytoLocation = NULL;
        omim->cytoLocation = NULL;
        omim->geneInheritance = NULL;
        omim->geneSymbols = NULL;
        omim->mimNumber = NULL;
        omim->phenotype = NULL;
        omim->phenotypeInheritance = NULL;
        omim->phenotypeMappingKey = NULL;
        omim->phenotypeMimNumber = NULL;
        omim->phenotypicSeriesNumber = NULL;
        omim->sequenceID = NULL;
        omim->transcript = NULL;
                

        MMALLOC(omim->chromosome,sizeof(char) * 128);//
        MMALLOC(omim->chromosomeLocationEnd,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->chromosomeLocationStart,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->chromosomeSort,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->chromosomeSymbol,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->computedCytoLocation,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->cytoLocation,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->geneInheritance,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->mimNumber,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->phenotype,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->phenotypeInheritance,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->phenotypeMappingKey,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->phenotypeMimNumber,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->phenotypicSeriesNumber,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->sequenceID,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->transcript,sizeof(char) * 128);// = NULL;

        

        omim->geneSymbols = NULL;
        MMALLOC(omim->geneSymbols,sizeof(char*) * MAX_ALT_GENE_NAMES);
        for(i = 0; i < MAX_ALT_GENE_NAMES;i++){
                omim->geneSymbols[i] = NULL;
                MMALLOC(omim->geneSymbols[i] ,sizeof(char) * 128);// = NULL;

        }
        RUN(clear_term(omim));
        return omim;
ERROR:
        free_omim(omim);
        return NULL;
}

int clear_term(struct OMIM* omim)
{
        int i;
        
        ASSERT(omim != NULL,"No term structure!");

        snprintf(omim->chromosome,128,"NA");
        snprintf(omim->chromosomeLocationEnd,128,"NA");
        snprintf(omim->chromosomeLocationStart,128,"NA");
        snprintf(omim->chromosomeSort,128,"NA");
        snprintf(omim->chromosomeSymbol,128,"NA");
        snprintf(omim->computedCytoLocation,128,"NA");
        snprintf(omim->cytoLocation,128,"NA");
        snprintf(omim->geneInheritance,128,"NA");
        snprintf(omim->mimNumber,128,"NA");
        snprintf(omim->phenotype,128,"NA");
        snprintf(omim->phenotypeInheritance,128,"NA");
        snprintf(omim->phenotypeMappingKey,128,"NA");
        snprintf(omim->phenotypeMimNumber,128,"NA");
        snprintf(omim->phenotypicSeriesNumber,128,"NA");
        snprintf(omim->sequenceID,128,"NA");
        snprintf(omim->transcript,128,"NA");

        for(i = 0; i < MAX_ALT_GENE_NAMES;i++){
                snprintf(omim->geneSymbols[i],128,"NA");
        }
        
        
        return OK;
ERROR:
        return FAIL; 
}

void free_omim(struct OMIM* omim)
{
        if(omim){
                int i;
                for(i = 0; i < MAX_ALT_GENE_NAMES;i++){
                        MFREE(omim->geneSymbols[i]);// ,sizeof(char) * 128);// = NULL;
                }
                MFREE(omim->geneSymbols);//,sizeof(char*) * 10);
        
                MFREE(omim->chromosome);
                MFREE(omim->chromosomeLocationEnd);
                MFREE(omim->chromosomeLocationStart);
                MFREE(omim->chromosomeSort);
                MFREE(omim->chromosomeSymbol);
                MFREE(omim->computedCytoLocation);
                MFREE(omim->cytoLocation);
                MFREE(omim->geneInheritance);
                MFREE(omim->mimNumber);
                MFREE(omim->phenotype);
                MFREE(omim->phenotypeInheritance);
                MFREE(omim->phenotypeMappingKey);
                MFREE(omim->phenotypeMimNumber);
                MFREE(omim->phenotypicSeriesNumber);
                MFREE(omim->sequenceID);
                MFREE(omim->transcript);
                MFREE(omim);
        }
}


