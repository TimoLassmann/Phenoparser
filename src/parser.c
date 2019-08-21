
#include "phenoparser.h"

#ifdef XML_LARGE_SIZE
#  if defined(XML_USE_MSC_EXTENSIONS) && _MSC_VER < 1400
#    define XML_FMT_INT_MOD "I64"
#  else
#    define XML_FMT_INT_MOD "ll"
#  endif
#else
#  define XML_FMT_INT_MOD "l"
#endif

#ifdef XML_UNICODE_WCHAR_T
#  define XML_FMT_STR "ls"
#else
#  define XML_FMT_STR "s"
#endif


static void startElement(void *userData, const XML_Char *name, const XML_Char **atts);
static void characterDataHandler(void *userData, const XML_Char *s, int len);
static void endElement(void *userData, const XML_Char *name);
static size_t parseStreamCallback(void *contents, size_t length, size_t nmemb, void *userp);

static int callback(void *NotUsed, int argc, char **argv, char **azColName);

int search_and_insert_disease(struct parameters* param, char* search_term, struct rbtree_root* series)
{
        char buffer[BUFFER_LEN*10];
        char search_term_tmp[BUFFER_LEN];
        CURL *curl_handle;
        CURLcode res;
        XML_Parser parser;
        struct ParserStruct state;
        int i,c,len;

        int start = 0;
        int step = 10;

        /* Initialize the state structure for parsing. */
        memset(&state, 0, sizeof(struct ParserStruct));
        state.ok = 1;
        state.ol = NULL;

        RUNP(state.ol = init_omim_list(50));

        len = strlen(search_term);

        c = 0;
        for(i = 0; i<= len;i++){
                if(isspace((int) search_term[i])){
                        search_term_tmp[c] = '%';
                        c++;
                        search_term_tmp[c] = '2';
                        c++;
                        search_term_tmp[c] = '0';
                        c++;
                }else{
                        search_term_tmp[c] = search_term[i];
                        c++;
                }
        }
        search_term_tmp[c] = 0;

        state.ol->num_entries = 1;
        while(state.ol->num_entries){
                RUN(clear_omim_list(state.ol));
                state.ok = 1;
                /* Initialize a namespace-aware parser. */
                parser = XML_ParserCreateNS(NULL, '\0');
                //LOG_MSG("STATE:%d", state.ok);
                XML_SetUserData(parser, &state);
                XML_SetElementHandler(parser, startElement, endElement);
                XML_SetCharacterDataHandler(parser, characterDataHandler);
                //LOG_MSG("STATE:%d", state.ok);
                /* Initialize a libcurl handle. */
                curl_global_init(CURL_GLOBAL_ALL ^ CURL_GLOBAL_SSL);
                curl_handle = curl_easy_init();

                snprintf(buffer,BUFFER_LEN*10,"https://api.omim.org/api/entry/search?search=%s%s%s+AND+gm_phenotype_exists:true&include=geneMap&start=%d&limit=%d&apiKey=%s","%22",search_term_tmp,"%22",start,step, param->omimkey);

                LOG_MSG("%s",buffer);
                DPRINTF2("%s",buffer);




                curl_easy_setopt(curl_handle, CURLOPT_URL,buffer);
                curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, parseStreamCallback);
                curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)parser);

                //printf("Depth   Characters   Closing Tag\n");
                //LOG_MSG("STATE:%d", state.ok);
                /* Perform the request and any follow-up parsing. */
                res = curl_easy_perform(curl_handle);
                //LOG_MSG("STATE:%d", state.ok);
                if(res != CURLE_OK) {
                        fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
                } else if(state.ok) {
                        /* Expat requires one final call to finalize parsing. */
                        if(XML_Parse(parser, NULL, 0, 1) == 0) {
                                int error_code = XML_GetErrorCode(parser);
                                fprintf(stderr, "Finalizing parsing failed with error code %d (%s).\n",
                                        error_code, XML_ErrorString(error_code));
                        }
                }


                /* insert results into tables...  */

                RUN(action_insert_into_sqlite(state.ol,param,search_term,series));



                start += step;

                XML_ParserFree(parser);
                curl_easy_cleanup(curl_handle);
                curl_global_cleanup();

                //free(state.characters.memory);
        }




        /* Clean up. */
        free(state.characters.memory);
        free_omim_list(state.ol);

        LOG_MSG("Done");
        return OK;
ERROR:
        return FAIL;
}

int phenotype_series_search(struct parameters* param, char* search_term)
{
        char buffer[BUFFER_LEN*10];
        char search_term_tmp[BUFFER_LEN];
        CURL *curl_handle;        CURLcode res;
        XML_Parser parser;
        struct ParserStruct state;
        int i,c,len;
        int start = 0;
        int step = 10;

        /* Initialize the state structure for parsing. */
        memset(&state, 0, sizeof(struct ParserStruct));
        state.ok = 1;
        state.ol = NULL;

        RUNP(state.ol = init_omim_list(50));

        /* Initialize a namespace-aware parser. */
        //parser = XML_ParserCreateNS(NULL, '\0');
        //XML_SetUserData(parser, &state);
        //XML_SetElementHandler(parser, startElement, endElement);
        //XML_SetCharacterDataHandler(parser, characterDataHandler);

        /* Initialize a libcurl handle. */
        //curl_global_init(CURL_GLOBAL_ALL ^ CURL_GLOBAL_SSL);
        //c/url_handle = curl_easy_init

        len = strlen(search_term);

        c = 0;
        for(i = 0; i<= len;i++){
                if(isspace((int) search_term[i])){
                        search_term_tmp[c] = '%';
                        c++;
                        search_term_tmp[c] = '2';
                        c++;
                        search_term_tmp[c] = '0';
                        c++;
                }else{
                        search_term_tmp[c] = search_term[i];
                        c++;
                }
        }
        search_term_tmp[c] = 0;

        state.ol->num_entries = 1;
        while(state.ol->num_entries){
                RUN(clear_omim_list(state.ol));
                state.ok = 1;
                /* Initialize a namespace-aware parser. */
                parser = XML_ParserCreateNS(NULL, '\0');
                XML_SetUserData(parser, &state);
                XML_SetElementHandler(parser, startElement, endElement);
                XML_SetCharacterDataHandler(parser, characterDataHandler);

                /* Initialize a libcurl handle. */
                curl_global_init(CURL_GLOBAL_ALL ^ CURL_GLOBAL_SSL);
                curl_handle = curl_easy_init();


                LOG_MSG("start search with %d genes.", state.ol->num_entries);
                snprintf(buffer,BUFFER_LEN*10,"http://api.omim.org/api/entry/search?search=phenotypic_series_number:%s%s%s&start=%d&limit=%d&include=geneMap&apiKey=%s","%22",search_term_tmp,"%22",start,step, param->omimkey);
                DPRINTF2("%s",buffer);





                curl_easy_setopt(curl_handle, CURLOPT_URL,buffer);
                curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, parseStreamCallback);
                curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)parser);

                //printf("Depth   Characters   Closing Tag\n");

                /* Perform the request and any follow-up parsing. */
                res = curl_easy_perform(curl_handle);
                if(res != CURLE_OK) {
                        fprintf(stderr, "curl_easy_perform() failed: %s\n",
                                curl_easy_strerror(res));
                } else if(state.ok) {
                        /* Expat requires one final call to finalize parsing. */
                        if(XML_Parse(parser, NULL, 0, 1) == 0) {
                                int error_code = XML_GetErrorCode(parser);
                                fprintf(stderr, "Finalizing parsing failed with error code %d (%s).\n",
                                        error_code, XML_ErrorString(error_code));
                        }
                }
                /* Insest query results into sqlite. Last parameter is null as I
                 * do not want anything extra to be inserted into the rb_tree
                 * holding phenotypic series. */
                RUN(action_insert_into_sqlite(state.ol, param, search_term, NULL));

                /*ol = state.ol;
                  LOG_MSG("Found: %d associations.",ol->num_entries);
                  c = 0;
                  for (i = 0; i <= ol->num_entries;i++){
                  for(j = 0; j < MAX_ALT_GENE_NAMES;j++){
                  if(strcmp(ol->terms[i]->geneSymbols[j],"NA")){
                  fprintf(stdout,"%s %s %s %s\n",ol->terms[i]->phenotype,ol->terms[i]->phenotypeMimNumber, ol->terms[i]->geneSymbols[j], ol->terms[i]->phenotypicSeriesNumber);
                  }
                  }
                  }*/
                LOG_MSG("Found: %d gene symbols.",c);
                start += step;

                //free(state.characters.memory);
                XML_ParserFree(parser);
                curl_easy_cleanup(curl_handle);
                curl_global_cleanup();


        }

        //XML_ParserFree(parser);
        /* Clean up. */
        free(state.characters.memory);
        free_omim_list(state.ol);
        //curl_easy_cleanup(curl_handle);
        //curl_global_cleanup();


        //XML_ParserFree(parser);
        //curl_easy_cleanup(curl_handle);
        //curl_global_cleanup();


        //free_omim_list(state.ol);
        //free(state.characters.memory);
        //XML_ParserFree(parser);
        //curl_easy_cleanup(curl_handle);
        //curl_global_cleanup();


        LOG_MSG("Done");
        return OK;
ERROR:
        return FAIL;
}





static void startElement(void *userData, const XML_Char *name, const XML_Char **atts)
{
        struct ParserStruct *state = (struct ParserStruct *) userData;
        state->tags++;
        state->depth++;

        /* Get a clean slate for reading in character data. */

        free(state->characters.memory);
        state->characters.memory = NULL;
        state->characters.size = 0;
}

static void characterDataHandler(void *userData, const XML_Char *s, int len)
{
        struct ParserStruct *state = (struct ParserStruct *) userData;
        struct MemoryStruct *mem = &state->characters;

        mem->memory = realloc(mem->memory, mem->size + len + 1);
        if(mem->memory == NULL) {
                /* Out of memory. */
                fprintf(stderr, "Not enough memory (realloc returned NULL).\n");
                state->ok = 0;
                return;
        }

        memcpy(&(mem->memory[mem->size]), s, len);
        mem->size += len;
        mem->memory[mem->size] = 0;
}

static void endElement(void *userData, const XML_Char *name)
{
        struct ParserStruct *state = (struct ParserStruct *) userData;
        //int i;

        state->depth--;

        /* THIS IS A HACK! the tags I wanty happen to be at depth 6...  */
        if(state->depth == 6){
                if (state->characters.size){
                        enter_term(state->ol,(char*) name,(char*)state->characters.memory );
                }
        }


        //       printf("%5lu   %10lu   NAME:%s\t%s\n", state->depth, state->characters.size, name, (char*)state->characters.memory);
}


static size_t parseStreamCallback(void *contents, size_t length, size_t nmemb, void *userp)
{
        XML_Parser parser = (XML_Parser) userp;
        size_t real_size = length * nmemb;
        struct ParserStruct *state = (struct ParserStruct *) XML_GetUserData(parser);

        /* Only parse if we're not already in a failure state. */
        if(state->ok && XML_Parse(parser, contents, real_size, 0) == 0) {
                int error_code = XML_GetErrorCode(parser);
                fprintf(stderr, "Parsing response buffer of length %lu failed"
                        " with error code %d (%s).\n",
                        real_size, error_code, XML_ErrorString(error_code));
                state->ok = 0;
        }
        return real_size;
}


int callback(void *NotUsed, int argc, char **argv, char **azColName)
{
        int i;
        static int first_call = 1;
        NotUsed = 0;

        if(first_call == 1){
                for (i = 0; i < argc; i++) {
                        if(i){
                                fprintf(stdout,",");
                        }
                        fprintf(stdout,"%s",azColName[i]);
                }
                fprintf(stdout,"\n");
        }
        for (i = 0; i < argc; i++) {
                if(i){
                        fprintf(stdout,",");
                }
                fprintf(stdout,"%s",argv[i] ? argv[i] : "NULL");
        }
        fprintf(stdout,"\n");
/*    for (i = 0; i < argc; i++) {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }*/
        first_call = 0;
        return OK;
}

