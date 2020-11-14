#include "tldevel.h"
#include "phenoparser.h"

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
