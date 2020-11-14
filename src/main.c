#include "tldevel.h"

#include "tlmisc.h"

#include "phenoparser.h"


#define str(x)          # x
#define xstr(x)         str(x)

char* safe_string_paste( int n_string, ... );

int query_omim_and_insert_results(struct parameters* param);

int make_table_output(struct parameters* param);

int read_phenolyzer_seed_gene_list(struct parameters* param);
int read_phenolyzer_merge_gene_scores(struct parameters* param);
int get_term_list(struct parameters* param);

static int insert_phenoterm_into_list(struct series_list* sl, char* term);

int main (int argc, char * argv[])
{
        struct parameters* param = NULL;

        if(argc == 1){
                /* print help */
                //tlog.echo_build_config();
                RUN(print_global_help(argc,argv));
                return OK;
        }
        if(strncmp(argv[1],"insert",6) == 0){
                /* get parameters */
                RUNP(param = get_insert_param(argc,argv));
                /* create output structure ...  */
                RUN(check_if_db_exists_otherwise_create(param));
                RUN(query_omim_and_insert_results(param));

        }else if(strncmp(argv[1],"termlist",8) == 0){

                RUNP(param = get_termlist_param(argc,argv));
                RUN(check_if_db_exists_otherwise_create(param));
                RUN(get_term_list(param));

        }else if(strncmp(argv[1],"readphe",7) == 0){

                RUNP(param = get_readphe_param(argc,argv));
                /* create output structure ...  */
                RUN(check_if_db_exists_otherwise_create(param));
                RUN(read_phenolyzer_seed_gene_list(param));
                RUN(read_phenolyzer_merge_gene_scores(param));

        }else if(strncmp(argv[1],"panel",5) == 0){

                RUNP(param = get_panel_param(argc,argv));
                RUN(make_table_output(param));

        }else{
                ERROR_MSG("Option %s not recognized.",argv[1]);
        }

        LOG_MSG("Done.");
        free_param(param);
        return EXIT_SUCCESS;
ERROR:
        free_param(param);
        return EXIT_FAILURE;

}

int read_phenolyzer_seed_gene_list(struct parameters* param)
{
        FILE* f_ptr = NULL;
        char line[BUFSIZ];
        int line_num = 1;
        int r;

        sqlite3 *sqlite_db = NULL;
        int rc;
        char buffer[BUFFER_LEN*10];

        int Rank;
        char Gene[BUFFER_LEN];
        int Id;
        double Score;
        char Status[BUFFER_LEN];
        char* sql_query = NULL;
        int sql_len = 0;


        ASSERT(param != NULL ,"No param.");

        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if(rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
        }

        //.merge_gene_scores
        //.seed_gene_list

        snprintf(buffer,BUFFER_LEN*10,"%s.seed_gene_list",param->phenofile);
        if(!my_file_exists(buffer)){
                ERROR_MSG("File: %s does not exist!",buffer);
        }

        RUNP(f_ptr = fopen(buffer,"r"));

        // SEED LIST looks like this:
        //Rank	Gene	ID	Score
        // 1 FGD1	2245	1
        // 2 UBE3A	7337	0.02885
        // FINAL list (I presume expanded like this:
        //Rank	Gene	ID	Score	Status
        //1	FGD1	2245	1	SeedGene
        //2	ELMO1	9844	0.06035	Predicted
        //3	PIKFYVE	200576	0.05363	Predicted
        //4	FGD3	89846	0.05362	Predicted

        line_num = 1;
        while(fgets(line, BUFSIZ, f_ptr)){
                if(line_num> 1){
                        r = sscanf(line,"%d %"xstr(BUFSIZ)"s %d %lf %"xstr(BUFSIZ)"s", &Rank, Gene, &Id,&Score,Status);
                        if(r == 4){
                                Status[0] = 'N';
                                Status[1] = 'A';
                                Status[2] = 0;

                        }

                        /* DPRINTF3("%s%d\t%s\t%d\t%e\t%s\n",line,Rank,Gene,Id,Score,Status); */
                        if(r != 4 && r != 5 ){
                                ERROR_MSG("Problem reading (%d): %s",r,line);
                        }
                        /* CREATE TABLE phenolyzer(patient_id TEXT NOT NULL,gene TEXT NOT NULL,identifier INT,score REAL, status  unique (patient_id,gene);  */
                        if(Rank <= 1000){
                                RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE INTO phenolyzer VALUES ('%s','%s','%d','%e','%s');",param->patient_id,Gene,Id,Score,Status ));
                                rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
                                if( rc!=SQLITE_OK ){
                                        LOG_MSG("try:%s",sql_query);
                                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                                }
                        }
                        //sscanf( dtm, "%s %s %d  %d", weekday, month, &day, &year );
                }

                line_num++;
        }
        fclose(f_ptr);

        //RUN(enter_detailed_evidence_information(param));

        rc = sqlite3_close(sqlite_db);
        if( rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));
        }
        MFREE(sql_query);
        return OK;
ERROR:
        if(sql_query){
                MFREE(sql_query);
        }
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

int read_phenolyzer_merge_gene_scores(struct parameters* param)
{
        FILE* f_ptr = NULL;
        char line[BUFSIZ];
        int line_num = 1;

        int new_gene = 1;

        char buffer[BUFSIZ*10];
        char append_buffer[BUFSIZ*10 + BUFSIZ*10];
        char gene[BUFSIZ];

        /* entries we want to capture */

        struct capture{
                char database_type[BUFSIZ]; /* Type of database : OMIM., umls, ORPHANET etc */
                char database_id[BUFSIZ];   /* ID number of entry i.e. C0265210 or  117550*/
                char data_base_name[BUFSIZ];      /* name of database - king of redundant */
                char disease_description[BUFSIZ]; /* name of disease  */
                char hpo_term[BUFSIZ];            /* HPO term triggering match  */
                double score;

        } cap;
        /* To gethyperlinks into the output datatable all we need to do us
         * insert html code into the information below and add the "escape =
         * FALSE" option to the R DT datatable call */

        char* description = NULL;
        char* sql_query = NULL;
        int description_len;
        int description_alloc_len = 256;
        int sql_len = 0;
        int i;

        sqlite3 *sqlite_db = NULL;
        int rc;

        int num_genes = 0;

        ASSERT(param != NULL ,"No param.");

        //.merge_gene_scores
        //.seed_gene_list

        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if(rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
        }

        snprintf(buffer,BUFFER_LEN*10,"%s.merge_gene_scores",param->phenofile);

        if(!my_file_exists(buffer)){
                ERROR_MSG("File: %s does not exist!",buffer);
        }

        MMALLOC(description,sizeof(char) * description_alloc_len);
        description_len = 0;
        //MMALLOC(sql_query_buffer,sizeof(char) *sql_query_buffer_alloc_len );


        RUNP(f_ptr = fopen(buffer,"r"));
        line_num = 1;
        while(fgets(line, BUFSIZ, f_ptr)){
                if(line_num> 1){
                        if(new_gene){
                                sscanf(line,"%"xstr(BUFFER_LEN)"s\t",gene);
                                //fprintf(stdout,"Gene: %s\n",gene);

                                new_gene = 0;
                                num_genes++;
                        }else{
                                if(strlen(line) == 1){
                                        new_gene = 1;
                                        description[description_len] = 0;

                                        RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE INTO phenolyzerGeneData  VALUES ('%s','%s','%s');",param->patient_id,gene,description));

                                        rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
                                        if( rc!=SQLITE_OK ){
                                                LOG_MSG("try:%s",sql_query);
                                                ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                                        }

                                        //fprintf(stdout,"%s (%d)\n%s\n",gene, num_genes,description);
                                        description_len = 0;
                                }else{
                                        /* This is interesting. Instead of using
                                         * %s etc is it possible to tell scanf
                                         * function to match everything
                                         * excluding a particular character. For
                                         * example %[^:]: means match everything
                                         * until a : is encountered. The next :
                                         * means scanf w2ill skip over this
                                         * character. */
                                        sscanf(line,"%"xstr(BUFFER_LEN)"[^:]:%"xstr(BUFFER_LEN)"[^ ] %"xstr(BUFFER_LEN)"[^\t]\t%"xstr(BUFFER_LEN)"[^\t]\t%"xstr(BUFFER_LEN)"[^\t]\t%lf", cap.database_type,cap.database_id,cap.data_base_name,cap.disease_description,cap.hpo_term,&cap.score);
                                        if(strncmp(cap.database_type,"OMIM",4) == 0){


                                                snprintf(buffer,BUFFER_LEN * 10,"<a href=\"https://www.omim.org/entry/%s\">%s</a>", cap.database_id,cap.database_id);

                                                /* safe_string_paste(5,"<a href=\"https://www.omim.org/entry/", */
                                                /*                   cap.database_id, */
                                                /*                   "\">", */
                                                /*                   cap.database_id, */
                                                /*                   "</a>"); */



                                        }else if(strncmp(cap.database_type,"ORPHANET",8) == 0){
                                                snprintf(buffer,BUFFER_LEN * 10,"<a href=\"http://www.orpha.net/consor/cgi-bin/OC_Exp.php?Lng=GB&Expert=%s\">%s</a>", cap.database_id,cap.database_id);
                                                /* safe_string_paste(5,"<a href=\"http://www.orpha.net/consor/cgi-bin/OC_Exp.php?Lng=GB&Expert=", */
                                                /*                   cap.database_id, */
                                                /*                   "\">", */
                                                /*                   cap.database_id, */
                                                /*                   "</a>"); */

                                        }else if(strncmp(cap.database_type,"umls",4) == 0){
                                                snprintf(buffer,BUFFER_LEN * 10,"<a href=\"https://www.ncbi.nlm.nih.gov/medgen/%s\">%s</a>", cap.database_id,cap.database_id);
                                                //          C0175695
                                                /* safe_string_paste(5,"<a href=\"https://www.ncbi.nlm.nih.gov/medgen/", */
                                                /*                   cap.database_id, */
                                                /*                   "\">", */
                                                /*                   cap.database_id, */
                                                /*                   "</a>"); */

                                        }else{
                                                //safe_string_paste(1, cap.database_id);
                                                snprintf(buffer,BUFFER_LEN * 10,"%s",cap.database_id);
                                        }


                                        /* create line parser line input to be concatenated with other info */

                                        snprintf(append_buffer,BUFFER_LEN*10 + BUFFER_LEN*10,"%s:%s %s %s %s %0.4f<br />",cap.database_type, buffer, cap.data_base_name, cap.disease_description,cap.hpo_term, cap.score);



                                        //fprintf(stdout,"%s\n%s\n%s\n%s\n%s\n%f (score) \n%s\n \n", cap.database_type, cap.database_id,cap.data_base_name, cap.disease_description, cap.hpo_term, cap.score,buffer);
                                        //fprintf(stdout,"%s\n",append_buffer);
                                        /* change omim id to hyperlink....  */
                                        for(i = 0; i < strlen(append_buffer);i++){
                                                description[description_len] = append_buffer[i];
                                                description_len++;
                                                if(description_len == description_alloc_len){
                                                        description_alloc_len = description_alloc_len << 1;
                                                        MREALLOC(description,sizeof(char) * description_alloc_len);
                                                }
                                        }
                                        append_buffer[0] = 0;
                                        /* DPRINTF3("%d %s",strlen(line),line); */
                                }
                        }
                        //DPRINTF3("%d %s",strlen(line),line);

                }
                line_num++;
        }
        fclose(f_ptr);

        MFREE(description);
        MFREE(sql_query);
        rc = sqlite3_close(sqlite_db);
        if( rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));

        }


        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        if(description){
                MFREE(description);
        }
        if(sql_query){
                MFREE(sql_query);
        }
        return FAIL;
}






int action_insert_into_sqlite( struct OMIM_list* ol,struct parameters* param, char* search_term, struct series_list* sl)
{
        sqlite3 *sqlite_db = NULL;
        int rc;
        char* sql_query = NULL;
        int sql_len = 0;
        int i,j,c;

        /* struct string_struct* tmp = NULL; */


        ASSERT(ol != NULL,"No List");
        ASSERT(param != NULL,"No parameters");


        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if(rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
        }

        RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE INTO patient  VALUES ('%s','%s');",param->patient_id,search_term));
        rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
        if( rc!=SQLITE_OK ){
                LOG_MSG("try:%s",sql_query);
                ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
        }


        LOG_MSG("Found: %d associations.",ol->num_entries);
        c = 0;
        for (i = 0; i <= ol->num_entries;i++){
                RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE INTO diseaseMIM VALUES ('%s','%s','%s','%s');", search_term,  ol->terms[i]->phenotypeMimNumber,ol->terms[i]->phenotype,ol->terms[i]->phenotypeInheritance ));

                rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
                if( rc != SQLITE_OK ){
                        LOG_MSG("try:%s",sql_query);
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }
                for(j = 0; j < MAX_ALT_GENE_NAMES;j++){
                        if(strcmp(ol->terms[i]->geneSymbols[j],"NA")){
                                c++;
                                RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE  INTO MIMgene VALUES ('%s','%s');", ol->terms[i]->phenotypeMimNumber,ol->terms[i]->geneSymbols[j]));
                                rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
                                if( rc!=SQLITE_OK ){
                                        LOG_MSG("try:%s",sql_query);
                                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                                }
                                if(sl){
                                        insert_phenoterm_into_list(sl, ol->terms[i]->phenotypicSeriesNumber);
                                }

                                /* fprintf(stdout,"%s %s %s %s %s\n",ol->terms[i]->phenotypeInheritance,  ol->terms[i]->phenotype,ol->terms[i]->phenotypeMimNumber, ol->terms[i]->geneSymbols[j], ol->terms[i]->phenotypicSeriesNumber); */
                        }
                }
        }
        LOG_MSG("Found: %d gene symbols.",c);

        rc = sqlite3_close(sqlite_db);
        MFREE(sql_query);
        return OK;
ERROR:
        if(sql_query){
                MFREE(sql_query);
        }
        return FAIL;
}

int insert_phenoterm_into_list(struct series_list* sl, char* term)
{
        /* check if already present */
        int i,c;
        if(!strcmp(term,"NA")){
                return OK;
        }
        c = 1;
        for(i = 0; i < sl->len;i++){
                if(!strncmp(sl->list[i],term, 128)){
                        c = 0;
                        break;
                }
        }
        if(c){
                sl->list[sl->len] = NULL;
                MMALLOC(sl->list[sl->len],sizeof(char) * 128);
                snprintf(sl->list[sl->len],128,"%s", term);
                sl->len++;
                if(sl->len == sl->alloc_len){
                        sl->alloc_len = sl->alloc_len + sl->alloc_len /2;
                        MREALLOC(sl->list, sizeof(char*) * sl->alloc_len);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int get_term_list(struct parameters* param)
{
        char* sql_query = NULL;
        int sql_len = 0;
        sqlite3 *sqlite_db = NULL;
        int rc;

        char** term_list = NULL;
        FILE* f_ptr = NULL;
        char line[BUFSIZ];

        int num_terms = 0;
        int alloc_numterms = 16;
        int i,j,c;
        int pos;

        ASSERT(param != NULL ,"No param.");
        ASSERT(param->local_sqlite_database_name != NULL,"No sqlite database name given.");

        MMALLOC(term_list,sizeof(char*) * alloc_numterms);
        for(i = 0; i < alloc_numterms;i++){
                term_list[i] = NULL;
                MMALLOC(term_list[i],sizeof(char) * BUFSIZ);
        }

        for(i = 0; i < param->num_infiles;i++){
                if(!my_file_exists(param->infile[i])){
                        WARNING_MSG("File: %s does not exist!", param->infile[i]);
                }else{
                        RUNP(f_ptr = fopen(param->infile[i],"r"));

                        while(fgets(line, BUFSIZ, f_ptr)){
                                /* DPRINTF3("%s",line); */
                                pos = 0;
                                for(j = 0;j < BUFSIZ;j++){
                                        /* DPRINTF3("%c",line[j]); */
                                        if(line[j] == ';'){
                                                term_list[num_terms][pos] = 0;
                                                /* DPRINTF2("%s",term_list[num_terms]); */

                                                pos = 0;
                                                num_terms++;
                                                if(num_terms == alloc_numterms){
                                                        alloc_numterms = alloc_numterms << 1;
                                                        /* DPRINTF3("Realloc:%d",alloc_numterms); */
                                                        MREALLOC(term_list,sizeof(char*) * alloc_numterms);
                                                        for(c = num_terms;c < alloc_numterms;c++){
                                                                term_list[c] = NULL;
                                                                MMALLOC(term_list[c],sizeof(char) * BUFSIZ);
                                                        }
                                                }
                                        }else if(iscntrl(line[j])){
                                                term_list[num_terms][pos] = 0;
                                                /* DPRINTF2("%s",term_list[num_terms]); */

                                                pos = 0;
                                                num_terms++;
                                                if(num_terms == alloc_numterms){
                                                        alloc_numterms = alloc_numterms << 1;
                                                        /* DPRINTF3("Realloc:%d",alloc_numterms); */
                                                        MREALLOC(term_list,sizeof(char*) * alloc_numterms);
                                                        for(c = num_terms;c < alloc_numterms;c++){
                                                                term_list[c] = NULL;
                                                                MMALLOC(term_list[c],sizeof(char) * BUFSIZ);
                                                        }
                                                }
                                                break;
                                        }else{
                                                term_list[num_terms][pos] = line[j];
                                                pos++;
                                        }
                                }

                        }
                        fclose(f_ptr);
                }
        }

        if(param->outfile){
                LOG_MSG("Writing to file: %s.",param->outfile);
                if(my_file_exists(param->outfile)){
                        WARNING_MSG("File: %s do exist; will overwrite!",param->outfile);
                }
                RUNP(f_ptr = fopen(param->outfile,"w"));
        }else{
                f_ptr = stdout;
        }

        for(i = 0 ; i< num_terms;i++){
                fprintf(f_ptr,"%s\n",term_list[i]);


        }
        if(param->outfile){
                fclose(f_ptr);
        }

        /* Insert into mysql databas */
        /* Don't need new table - just insert into patient */

        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if(rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
        }

        for(i = 0 ; i< num_terms;i++){
                RUNP(sql_query = create_query_string(sql_query,&sql_len,"INSERT OR IGNORE INTO patient  VALUES ('%s','%s');",param->patient_id,term_list[i]));
                rc = sqlite3_exec(sqlite_db, sql_query, 0, 0, 0);
                if( rc!=SQLITE_OK ){
                        LOG_MSG("try:%s",sql_query);
                        ERROR_MSG("sqlite3_exec failed: %s\n", sqlite3_errmsg(sqlite_db));
                }

        }

        rc = sqlite3_close(sqlite_db);
        if( rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));

        }


        for(i = 0; i < alloc_numterms;i++){
                MFREE(term_list[i]);
        }
        MFREE(term_list);
        MFREE(sql_query);

        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        if(term_list){
                for(i = 0; i < alloc_numterms;i++){
                        MFREE(term_list[i]);
                }
                MFREE(term_list);
        }
        if(sql_query){
                MFREE(sql_query);
        }
        return FAIL;
}

int make_table_output(struct parameters* param)
{
        FILE* fptr = NULL;
        sqlite3 *sqlite_db = NULL;
        int rc;

        sqlite3_stmt *pStmt;
        char* buffer = NULL;
        int buffer_len = 1024;

        ASSERT(param != NULL, "No parameters.");
        ASSERT(param->outfile != NULL ,"No outfile.");

        MMALLOC(buffer,sizeof(char) * buffer_len);

        snprintf(buffer,buffer_len,"%s_omim.csv",param->outfile);
        RUNP(fptr = fopen(buffer,"w"));



        snprintf(buffer,buffer_len,"SELECT patient_id AS ID, patient.DiseaseSearch AS DISEASE,  MIMgene.phenotypeMimNumber AS MIMnumber, diseaseMIM.phenotypeDescription AS DESC, diseaseMIM.phenotypeInheritance AS DESC, MIMgene.gene AS GENE  FROM  patient INNER JOIN diseaseMIM ON diseaseMIM.DiseaseSearch = patient.DiseaseSearch INNER JOIN MIMgene ON MIMgene.phenotypeMimNumber = diseaseMIM.phenotypeMimNumber   WHERE patient_id == \"%s\";",param->patient_id);

        rc = sqlite3_open(param->local_sqlite_database_name, &sqlite_db);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_open failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }

        rc = sqlite3_prepare(sqlite_db, buffer, -1, &pStmt, 0);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_prepare failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }

        while ( sqlite3_step(pStmt) !=SQLITE_DONE) {
                int i;
                int num_cols = sqlite3_column_count(pStmt);

                for (i = 0; i < num_cols; i++)
                {
                        if(i){
                                /* fprintf(fptr,","); */
                                fprintf(fptr,"\t");
                        }
                        switch (sqlite3_column_type(pStmt, i))
                        {
                        case (SQLITE3_TEXT):
                                fprintf(fptr,"%s", sqlite3_column_text(pStmt, i));
                                break;
                        case (SQLITE_INTEGER):
                                fprintf(fptr,"%d", sqlite3_column_int(pStmt, i));
                                break;
                        case (SQLITE_FLOAT):
                                fprintf(fptr,"%g", sqlite3_column_double(pStmt, i));
                                break;
                        default:
                                break;
                        }
                }
                fprintf(fptr,"\n");
        }


        sqlite3_finalize(pStmt);
        fclose(fptr);
        /*
          rc = sqlite3_exec(sqlite_db, buffer, callback, 0, 0);

          if (rc != SQLITE_OK ) {
          fprintf(stderr, "Failed to select data\n");
          }*/

        // Print Phenolyzer outpuit...
        snprintf(buffer,buffer_len,"%s_phenolyzer.csv",param->outfile);

        RUNP(fptr = fopen(buffer,"w"));





        //snprintf(buffer,buffer_len,"SELECT * FROM  phenolyzer   WHERE patient_id == \"%s\" ORDER BY score DESC;",param->patient_id);


        snprintf(buffer,buffer_len,"SELECT phenolyzer.patient_id ,phenolyzer.gene ,phenolyzer.identifier ,phenolyzer.score , phenolyzerGeneData.gene_evidence   FROM  phenolyzer INNER JOIN phenolyzerGeneData ON (phenolyzer.gene == phenolyzerGeneData.gene AND  phenolyzer.patient_id == phenolyzerGeneData.patient_id) WHERE phenolyzer.patient_id == \"%s\" ORDER BY score DESC;",param->patient_id);


        rc = sqlite3_prepare(sqlite_db, buffer, -1, &pStmt, 0);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_prepare failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }

        while ( sqlite3_step(pStmt) !=SQLITE_DONE) {
                int i;
                int num_cols = sqlite3_column_count(pStmt);

                for (i = 0; i < num_cols; i++)
                {
                        if(i){
                                /* fprintf(fptr,","); */
                                fprintf(fptr,"\t");
                        }
                        switch (sqlite3_column_type(pStmt, i))
                        {
                        case (SQLITE3_TEXT):
                                fprintf(fptr,"%s", sqlite3_column_text(pStmt, i));
                                break;
                        case (SQLITE_INTEGER):
                                fprintf(fptr,"%d", sqlite3_column_int(pStmt, i));
                                break;
                        case (SQLITE_FLOAT):
                                fprintf(fptr,"%g", sqlite3_column_double(pStmt, i));
                                break;
                        default:
                                break;
                        }
                }
                fprintf(fptr,"\n");
        }


        sqlite3_finalize(pStmt);


        /* write out table with terms used...  */
        // Print Phenolyzer outpuit...
        snprintf(buffer,buffer_len,"%s_terms.csv",param->outfile);

        RUNP(fptr = fopen(buffer,"w"));





        snprintf(buffer,buffer_len,"SELECT * FROM  patient  WHERE patient_id == \"%s\";",param->patient_id);


        rc = sqlite3_prepare(sqlite_db, buffer, -1, &pStmt, 0);
        if( rc!=SQLITE_OK ){
                fprintf(stderr, "sqlite3_prepare failed: %s\n", sqlite3_errmsg(sqlite_db));
                sqlite3_close(sqlite_db);
                exit(1);
        }

        while ( sqlite3_step(pStmt) !=SQLITE_DONE) {
                int i;
                int num_cols = sqlite3_column_count(pStmt);

                for (i = 0; i < num_cols; i++)
                {
                        if(i){
                                /* fprintf(fptr,","); */
                                fprintf(fptr,"\t");
                        }
                        switch (sqlite3_column_type(pStmt, i))
                        {
                        case (SQLITE3_TEXT):
                                fprintf(fptr,"%s", sqlite3_column_text(pStmt, i));
                                break;
                        case (SQLITE_INTEGER):
                                fprintf(fptr,"%d", sqlite3_column_int(pStmt, i));
                                break;
                        case (SQLITE_FLOAT):
                                fprintf(fptr,"%g", sqlite3_column_double(pStmt, i));
                                break;
                        default:
                                break;
                        }
                }
                fprintf(fptr,"\n");
        }


        sqlite3_finalize(pStmt);
        rc = sqlite3_close(sqlite_db);
        if( rc!=SQLITE_OK ){
                ERROR_MSG("sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));

        }
        fclose(fptr);
        MFREE(buffer);
        return OK;
ERROR:
        if(fptr){
                fclose(fptr);
        }
        if(sqlite_db){
                rc = sqlite3_close(sqlite_db);
                if( rc!=SQLITE_OK ){
                        fprintf(stderr, "sqlite3_close failed: %s\n", sqlite3_errmsg(sqlite_db));
                }
        }
        if(buffer){
                MFREE(buffer);
        }
        return FAIL;
}



int query_omim_and_insert_results(struct parameters* param)
{
        FILE* fptr = NULL;
        /* struct rbtree_root* series = NULL; */
        /* struct string_struct* tmp = NULL; */
        char buffer[BUFFER_LEN];
        int i;
        int len;

        ASSERT(param != NULL, "No parameters.");

        struct series_list* sl = NULL;

        MMALLOC(sl, sizeof(struct series_list));
        sl->alloc_len = 256;
        sl->len = 0;
        sl->list = NULL;
        MMALLOC(sl->list, sizeof(char*) * sl->alloc_len);

        /* RUNP(series = make_string_tree()); */
        //RUN(phenotype_series_search(param,"PS106210"));



        LOG_MSG("Phenotype search done");


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
                RUN(search_and_insert_disease(param,buffer+len,  sl));
        }




        //series->print_tree(series,NULL);


        /* series->flatten_tree(series); */
        /* ASSERT((series->cur_data_nodes == series->num_entries),"fail"); */
        for(i = 0; i < sl->len ;i++){

                LOG_MSG("Searching for terms in phenotypic series: %s", sl->list[i]);
                RUN(phenotype_series_search(param, sl->list[i]));
        }
        fclose(fptr);

        if(sl){
                for(i = 0; i < sl->len ;i++){
                        MFREE(sl->list[i]);
                }
                MFREE(sl->list);
                MFREE(sl);
        }

        // this stuff is allocated by libcrypto and libssl when curl is called by not freed.
        // this is no problem but when debugging with valgrind I'l like to see 'my' mem leaks!!
        //CONF_modules_free();
        //ERR_remove_state(0);
        //ENGINE_cleanup();
        //CONF_modules_unload(1);
        //ERR_free_strings();
        //EVP_cleanup();
        //CRYPTO_cleanup_all_ex_data();
        //sk_SSL_COMP_free(SSL_COMP_get_compression_methods());

        return OK;
ERROR:
        if(fptr){
                fclose(fptr);
        }
        return FAIL;
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
        /* DPRINTF3("entering:%s %s %d\n",name,value,new); */
        if(strcmp(name,"mimNumber") == 0 ) {
                /* DPRINTF3("found matching tag :%s %s %d\n",name,value,new); */
                /* DPRINTF3("existing:%s \n",ol->terms[cur_OMIM]->mimNumber); */


                if(strcmp(ol->terms[cur_OMIM]->mimNumber,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->mimNumber,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->mimNumber));

                }else{
                        new = 1;
                }

        }
        if(strcmp(name,"phenotype") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotype,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotype,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->phenotype));
                }else{
                        new = 1;
                }
        }
        if(strcmp(name,"phenotypeMimNumber") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeMimNumber,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeMimNumber,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->phenotypeMimNumber));
                }else{
                        new = 1;
                }
        }
        if(strcmp(name,"phenotypeMappingKey") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeMappingKey ,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeMappingKey,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->phenotypeMappingKey));
                }else{
                        new = 1;
                }
        }
        if(strcmp(name,"phenotypeInheritance") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypeInheritance,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypeInheritance,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->phenotypeInheritance));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"phenotypicSeriesNumber") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->phenotypicSeriesNumber,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->phenotypicSeriesNumber,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->phenotypicSeriesNumber));
                }else{
                        new = 1;
                }
        }



        if(strcmp(name,"sequenceID") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->sequenceID,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->sequenceID,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->sequenceID));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"chromosome") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosome,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosome,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->chromosome));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"chromosomeSymbol") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeSymbol,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeSymbol,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->chromosomeSymbol));
                }else{
                        new = 1;
                }
        }
        if(strcmp(name,"chromosomeSort") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeSort,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeSort,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->chromosomeSort));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"chromosomeLocationStart") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeLocationStart,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeLocationStart,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->chromosomeLocationStart));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"chromosomeLocationEnd") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->chromosomeLocationEnd,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->chromosomeLocationEnd,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->chromosomeLocationEnd));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"transcript") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->transcript,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->transcript,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->transcript));
                }else{
                        new = 1;
                }
        }
        if(strcmp(name,"cytoLocation") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->cytoLocation,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->cytoLocation,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->cytoLocation));
                }else{
                        new = 1;
                }
        }

        if(strcmp(name,"computedCytoLocation") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->computedCytoLocation,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->computedCytoLocation,"%s", value);
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->computedCytoLocation));
                }else{
                        new = 1;
                }
        }


        if(strcmp(name,"geneInheritance") == 0 ) {
                if(strcmp(ol->terms[cur_OMIM]->geneInheritance,"NA")==0){
                        sprintf(ol->terms[cur_OMIM]->geneInheritance,"%s", value);
                        /* DPRINTF3("gene inheritance: %d.",ol->terms[cur_OMIM]->geneInheritance ); */
                        RUN(remove_comma_quote(ol->terms[cur_OMIM]->geneInheritance));
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
        /* DPRINTF3("entering:%s %s %d\n",name,value,new); */
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






char* safe_string_paste( int n_string, ... )
{
        char* tmp = NULL;
        int len = 0;
        int i;

        va_list ap;
        va_start( ap, n_string);
        for(i = 0; i < n_string; i++){
                len += strlen(va_arg(ap, char* ));

        }
        va_end(ap);

        /* plus one for terminator */
        len++;

        MMALLOC(tmp, sizeof(char) * len);
        tmp[0] = 0;

        va_start( ap, n_string);
        for(i = 0; i < n_string; i++){
                tmp = strncat(tmp, va_arg(ap, char*), BUFFER_LEN);
                //strlen += strnlen(va_arg( intArgumentPointer, int ), BUFFER_LEN);

        }
        va_end(ap);
        return tmp;
ERROR:
        return NULL;
}
