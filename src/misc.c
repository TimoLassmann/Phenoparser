
#include "phenoparser.h"

long int compare_name(void* keyA, void* keyB);

static int resolve_default(void* ptr_a,void* ptr_b);

static void* get_name(void* ptr);

void print_string_struct(void* ptr,FILE* out_ptr);
void free_string_struct(void* ptr);

static void* get_name(void* ptr)
{
        struct string_struct* tmp = (struct string_struct*)  ptr;
        return tmp->name;
}

static int resolve_default(void* ptr_a,void* ptr_b)
{
        struct string_struct* tmp = (struct string_struct*)  ptr_b;
        if(tmp){
                if(tmp->name){
                        MFREE(tmp->name);
                }
                MFREE(tmp);
        }
        return 0;
}


void free_string_struct(void* ptr)
{
        struct string_struct* tmp = (struct string_struct*)  ptr;
        if(tmp){
                if(tmp->name){
                        MFREE(tmp->name);
                }
                MFREE(tmp);
        }
}


long int compare_name(void* keyA, void* keyB)
{
        return strcmp(keyA,keyB);
}



void print_string_struct(void* ptr,FILE* out_ptr)
{
        struct string_struct* tmp = (struct string_struct*)  ptr;
        fprintf(out_ptr,"%s\n",tmp->name );
}

struct rbtree_root* make_string_tree(void)
{
        struct rbtree_root* root = NULL;

        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;

        void (*fp_free)(void* ptr) = NULL;
	

        fp_get = &get_name;
        fp_cmp = &compare_name;
        fp_print = &print_string_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_string_struct;
	
        RUNP(root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free));
        return root;
ERROR:
        return NULL;
}


int remove_comma_quote(char* in)
{
        int i;
        int len;
        len = strlen(in);
        for(i = 0; i < len;i++){
                if(in[i] == ','){
                        in[i] = '_';
                }
                if((int)in[i] == 39 || (int)in[i] == 64 || (int)in[i] ==  96){
                        in[i] = '_';
                }
        }
        return OK;
}

