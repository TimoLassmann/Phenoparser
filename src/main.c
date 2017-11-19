#include "phenoparser.h"


int check_if_db_exists_otherwise_create(struct parameters* param);

int query_OMIM_and_insert_results(struct parameters* param);
int search_and_insert_disease(struct parameters* param, char* search_term);

int make_table_output(struct parameters* param);

static void startElement(void *userData, const XML_Char *name, const XML_Char **atts);
static void characterDataHandler(void *userData, const XML_Char *s, int len);
static void endElement(void *userData, const XML_Char *name);
static size_t parseStreamCallback(void *contents, size_t length, size_t nmemb, void *userp);

int enter_term(struct OMIM_list* ol, char* name,char* value );

int callback(void *NotUsed, int argc, char **argv, char **azColName);


struct OMIM_list* init_omim_list(int n);
int resize_omim_list(struct OMIM_list* ol, int add);
void free_omim_list(struct OMIM_list* ol);


struct OMIM* init_omim_entry(void);
void free_omim(struct OMIM* omim);


int main (int argc, char * argv[])
{
        struct parameters* param = NULL;
        
        
        if(argc == 1){
                /* print help */
                tlog.echo_build_config();
                RUN(print_global_help(argc,argv));
                return OK;
        }
        
        if(strncmp(argv[1],"insert",6) == 0){
                /* get parameters */
                RUNP(param = get_insert_param(argc,argv));
                /* create output structure ...  */
                RUN(check_if_db_exists_otherwise_create(param));

                RUN(query_OMIM_and_insert_results(param));

                
        }else if(strncmp(argv[1],"panel",5) == 0){
                RUNP(param = get_panel_param(argc,argv));
                RUN(make_table_output(param));
                        
        }else{
                ERROR_MSG("Option %s not recognized.",argv[1]);
        }
        free_param(param);
        return EXIT_SUCCESS;
ERROR:        
        free_param(param);
        return EXIT_FAILURE;
}

int make_table_output(struct parameters* param)
{
        sqlite3 *sqlite_db = NULL;
        int rc;          

        char* buffer = NULL;
        int buffer_len = 1024;     
 
        ASSERT(param != NULL, "No parameters.");

        MMALLOC(buffer,sizeof(char) * buffer_len);
       
        snprintf(buffer,buffer_len,"SELECT patient_id AS ID, patient.DiseaseSearch AS DISEASE,  MIMgene.phenotypeMimNumber AS MIMnumber, diseaseMIM.phenotypeDescription AS DESC,  MIMgene.gene AS GENE  FROM  patient INNER JOIN diseaseMIM ON diseaseMIM.DiseaseSearch = patient.DiseaseSearch INNER JOIN MIMgene ON MIMgene.phenotypeMimNumber = diseaseMIM.phenotypeMimNumber   WHERE patient_id == \"%s\";",param->patient_id);

        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }


        rc = sqlite3_exec(sqlite_db, buffer, callback, 0, 0);
    
        if (rc != SQLITE_OK ) {
                fprintf(stderr, "Failed to select data\n");
        } 
    
        rc = sqlite3_close(sqlite_db);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }
        return OK;
ERROR:
        return FAIL;
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


int query_OMIM_and_insert_results(struct parameters* param)
{
        FILE* fptr = NULL;
        char buffer[BUFFER_LEN];

        int len;
        
        ASSERT(param != NULL, "No parameters.");
        
        RUNP(fptr = fopen(param->phenofile,"r"));


        while (fgets(buffer, sizeof(buffer), fptr)){

                len = strlen(buffer);
     
                while(isspace((int)buffer[len]) || iscntrl((int)buffer[len])){
                        buffer[len] = 0;
                        len--;
                }
                len = 0;
                while(isspace((int)buffer[len]) || iscntrl((int)buffer[len])){
                        len++;
                }
                
                
                LOG_MSG("Searching with: \"%s\"",buffer+len);
                RUN(search_and_insert_disease(param,buffer+len));
        }
        
        
        


        fclose(fptr);

        return OK;
ERROR:
        if(fptr){
                fclose(fptr);
        }
        return FAIL;
}

int search_and_insert_disease(struct parameters* param, char* search_term)
{
        struct OMIM_list* ol = NULL;
        

        sqlite3 *sqlite_db = NULL;
        int rc;          

        char buffer[BUFFER_LEN*10];
        char search_term_tmp[BUFFER_LEN];
        CURL *curl_handle;
        CURLcode res;
        XML_Parser parser;
        struct ParserStruct state;
        int i,j,c,len; 
        /* Initialize the state structure for parsing. */
        memset(&state, 0, sizeof(struct ParserStruct));
        state.ok = 1;
        state.ol = NULL; 

        RUNP(state.ol = init_omim_list(50));
        
        /* Initialize a namespace-aware parser. */
        parser = XML_ParserCreateNS(NULL, '\0');
        XML_SetUserData(parser, &state);
        XML_SetElementHandler(parser, startElement, endElement);
        XML_SetCharacterDataHandler(parser, characterDataHandler);
 
        /* Initialize a libcurl handle. */
        curl_global_init(CURL_GLOBAL_ALL ^ CURL_GLOBAL_SSL);
        curl_handle = curl_easy_init();


        
        
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
        
        snprintf(buffer,BUFFER_LEN*10,"http://api.omim.org/api/entry/search?search=%s%s%s+AND+gm_phenotype_exists:true&include=geneMap&apiKey=%s","%22",search_term_tmp,"%22", param->omimkey);
        DPRINTF2("%s",buffer);
        //curl_easy_setopt(curl_handle, CURLOPT_URL,"http://api.omim.org/api/entry/search?search=%22epidermolysis%20bullosa%22+AND+gm_phenotype_exists:true&include=geneMap&apiKey=odBmVLYxRLO11WyREKlwyw");
        curl_easy_setopt(curl_handle, CURLOPT_URL,buffer);
        curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, parseStreamCallback);
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)parser);
 
        //       printf("Depth   Characters   Closing Tag\n");
 
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
        /* insert results into tables...  */
        
        
        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if(rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
        }
        snprintf(buffer,BUFFER_LEN*10,"INSERT OR IGNORE INTO patient  VALUES ('%s','%s');",param->patient_id,search_term);
        rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
        if( rc!=SQLITE_OK ){
                LOG_MSG("try:%s",buffer);            
                ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
        }        	
        

        ol = state.ol;
        LOG_MSG("Found: %d associations.",ol->num_entries);
        c = 0;
        for (i = 0; i <= ol->num_entries;i++){
                snprintf(buffer,BUFFER_LEN*10,"INSERT OR IGNORE INTO diseaseMIM VALUES ('%s','%s','%s');", search_term,  ol->terms[i]->phenotypeMimNumber,ol->terms[i]->phenotype);
                
                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        LOG_MSG("try:%s",buffer);
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }
                for(j = 0; j < MAX_ALT_GENE_NAMES;j++){
                        if(strcmp(ol->terms[i]->geneSymbols[j],"NA")){
                                c++;
                                snprintf(buffer,BUFFER_LEN*10,"INSERT OR IGNORE  INTO MIMgene VALUES ('%s','%s');", ol->terms[i]->phenotypeMimNumber,ol->terms[i]->geneSymbols[j]);
                                rc = sqlite3_exec(sqlite_db, buffer, 0, 0, 0);
                                if( rc!=SQLITE_OK ){
                                        LOG_MSG("try:%s",buffer);
                                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                                }
                                // fprintf(stdout,"%s %s %s\n",ol->terms[i]->phenotype,ol->terms[i]->phenotypeMimNumber, ol->terms[i]->geneSymbols[j]);
                        }
                }
        }
        LOG_MSG("Found: %d gene symbols.",c);
              
        rc = sqlite3_close(sqlite_db);
        
        
        
        /* Clean up. */

        free_omim_list(state.ol);
        free(state.characters.memory);
        XML_ParserFree(parser);
        curl_easy_cleanup(curl_handle);
        curl_global_cleanup();
        LOG_MSG("Done");
        return OK;
ERROR:
        return FAIL;
}


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
                snprintf(buffer,BUFFER_LEN,"CREATE TABLE diseaseMIM(DiseaseSearch TEXT, phenotypeMimNumber INT NOT NULL,phenotypeDescription TEXT, unique (DiseaseSearch,phenotypeMimNumber));");
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

  WORKs ! 


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








  SELECT
  trackid,
  tracks.name AS Track,
  albums.title AS Album,
  artists.name AS Artist
  FROM
  tracks
  INNER JOIN albums ON albums.albumid = tracks.albumid
  INNER JOIN artists ON artists.artistid = albums.artistid;


*/


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

int enter_term(struct OMIM_list* ol, char* name,char* value )
{
        int cur_OMIM = ol->num_entries;
        int new = 0;

        ASSERT(ol != NULL,"No list");
        /*
          NAME:mimNumber	601282
          6           53   NAME:phenotype	Epidermolysis bullosa simplex with muscular dystrophy
          6            6   NAME:phenotypeMimNumber	226670
   

        */
        DPRINTF3("entering:%s %s %d\n",name,value,new);
        if(strcmp(name,"mimNumber") == 0 ) {
                DPRINTF3("found matching tag :%s %s %d\n",name,value,new);
                DPRINTF3("existing:%s \n",ol->terms[cur_OMIM]->mimNumber);

                
                if(strcmp(ol->terms[cur_OMIM]->mimNumber,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->mimNumber,"%s", value);
                }else{
                        new = 1;  
                }
         
        }
        if(strcmp(name,"phenotype") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotype,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotype,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"phenotypeMimNumber") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeMimNumber,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeMimNumber,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"phenotypeMappingKey") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeMappingKey ,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeMappingKey,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"phenotypeInheritance") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeInheritance,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeInheritance,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"sequenceID") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->sequenceID,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->sequenceID,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"chromosome") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosome,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosome,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"chromosomeSymbol") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeSymbol,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeSymbol,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"chromosomeSort") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeSort,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeSort,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"chromosomeLocationStart") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeLocationStart,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeLocationStart,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"chromosomeLocationEnd") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeLocationEnd,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeLocationEnd,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"transcript") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->transcript,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->transcript,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"cytoLocation") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->cytoLocation,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->cytoLocation,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"computedCytoLocation") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->computedCytoLocation,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->computedCytoLocation,"%s", value);
                }else{
                        new = 1;  
                }                
        }

        if(strcmp(name,"geneInheritance") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->geneInheritance,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->geneInheritance,"%s", value);
                }else{
                        new = 1;  
                }                
        }
        if(strcmp(name,"geneSymbols") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->geneSymbols[0],"NA")==0){
                        char seps[]  = " ,\t\n";

                        char* token;
                        char var[64];
                       
                        int i = 0;

                        token = strtok (value, seps);
                        while (token != NULL)
                        {
                                sscanf (token, "%s", var);
                                snprintf(ol->terms[cur_OMIM]->geneSymbols[i],128,"%s",var);

                                token = strtok (NULL, seps);
                                i++;
                                if(i == MAX_ALT_GENE_NAMES){
                                        break;
                                }
                        }


                      
                }else{
                        new = 1;  
                }                
        }
        DPRINTF3("entering:%s %s %d\n",name,value,new);
        if(new == 1){
                if(cur_OMIM + 1 == ol->num_malloced){
                        resize_omim_list(ol,10);
                }
                ol->num_entries++;
                RUN(enter_term(ol,name,value ));
                
        }

        return OK;
ERROR:
        return FAIL; 
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
        MMALLOC(omim->sequenceID,sizeof(char) * 128);// = NULL;
        MMALLOC(omim->transcript,sizeof(char) * 128);// = NULL;

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
        snprintf(omim->sequenceID,128,"NA");
        snprintf(omim->transcript,128,"NA");


        

        omim->geneSymbols = NULL;
        MMALLOC(omim->geneSymbols,sizeof(char*) * MAX_ALT_GENE_NAMES);
        for(i = 0; i < MAX_ALT_GENE_NAMES;i++){
                omim->geneSymbols[i] = NULL;
                MMALLOC(omim->geneSymbols[i] ,sizeof(char) * 128);// = NULL;
                snprintf(omim->geneSymbols[i],128,"NA");
        }
        
        return omim;
ERROR:
        free_omim(omim);
        return NULL;
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
                MFREE(omim->sequenceID);
                MFREE(omim->transcript);
                MFREE(omim);
        }
}


