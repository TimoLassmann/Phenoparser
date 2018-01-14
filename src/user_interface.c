
#define OPT_OMIM_KEY 1
#define OPT_PATIENT_ID 2
#define OPT_SQLITE_DB 3
#define OPT_PHENOFILE 4
#define OPT_OUTFILE 5


#include "phenoparser.h"


static struct parameters* init_param(void);

static int byg_count(char* pattern, char* text);
static int print_insert_help(int argc, char * argv[]);
static int print_panel_help(int argc, char * argv[]);
static int print_termlist_help(int argc, char * argv[]);
int print_readphe_help(int argc, char * argv[]);

int print_global_help(int argc, char * argv[])
{
        const char usage[] = " <command>";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);	
        fprintf(stdout,"Commands:\n\n");
	
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"insert","Inserts patient data and populates OMIM tables." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"panel","Retrieve gene panel from local database." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"termlist","Merges terms from multiple files." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"readphe","Reads phenolyzer output." ,"[NA]"  );



        
        fprintf(stdout,"\n");
        
        return OK; 
}


int print_panel_help(int argc, char * argv[])
{
        const char usage[] = " panel [-options] ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);	
        fprintf(stdout,"Options:\n\n");
	
       
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--id","Patient ID." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--db","Local database name" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--out","Output file name PREFIX" ,"[NA]"  );
        fprintf(stdout,"\n");
        return OK;
}

int print_insert_help(int argc, char * argv[])
{
        const char usage[] = " insert [-options] ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);	
        fprintf(stdout,"Options:\n\n");
	
       
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--id","Patient ID." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--pheno","File with patient disease terms" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--key","OMIM key." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--db","Local database name" ,"[NA]"  );
        fprintf(stdout,"\n");
        return OK;
}

int print_termlist_help(int argc, char * argv[])
{
        const char usage[] = " termlist <file1> <file2> ...  ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--out","Output file name" ,"[NA]"  );
       
        return OK;
}

int print_readphe_help(int argc, char * argv[])
{
        const char usage[] = " readphe [-options] ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);	
        fprintf(stdout,"Options:\n\n");
	
       
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--id","Patient ID." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--pheno","File with patient disease terms" ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--db","Local database name" ,"[NA]"  );
        fprintf(stdout,"\n");
        fprintf(stdout,"\nNOTE:\n\tThis program will only import the top 1000 genes.\n");
        
        return OK;
}

struct parameters* get_readphe_param(int argc, char * argv[])
{
        struct parameters* param = NULL;     
        int c, help; 

        RUN(print_program_header(argv,"Loads Phenolyzer results into a database."));

        
        if(argc == 2){
                RUN(print_readphe_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }


        help = 0;
        c = 0;

        RUNP(param = init_param());

        
        while (1){
                static struct option long_options[] ={
                        {"db",required_argument,0,OPT_SQLITE_DB},
                        {"id",required_argument,0, OPT_PATIENT_ID},
                        {"pheno",required_argument,0, OPT_PHENOFILE},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_PHENOFILE:
                        param->phenofile = optarg;
                        break;
                case OPT_PATIENT_ID:
                        param->patient_id = optarg;
                        break;
                case OPT_SQLITE_DB:
                        param->local_sqlite_database_name = optarg;
                        break;
                        
                case 'h':
                        help = 1;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }

        if(help){
                RUN(print_readphe_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }
        ASSERT(param->local_sqlite_database_name != NULL,"No database.");
        ASSERT(param->patient_id != NULL,"No patient.");
        ASSERT(param->phenofile != NULL,"No phenotype file.");
        
        
        RUN(log_command_line(argc, argv));
        
        LOG_MSG("Read param.");       
        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}

struct parameters* get_termlist_param(int argc, char * argv[])
{
        struct parameters* param = NULL;
        int c, help;
        RUN(print_program_header(argv,"Collects terms from files and merges them."));
        help = 0;
        c = 0;

        if(argc == 2){
                RUN(print_termlist_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }
        
        RUN(log_command_line(argc, argv));


        c = 0;
        RUNP(param = init_param());
        while (1){
                static struct option long_options[] ={
                        {"out",required_argument,0, OPT_OUTFILE},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case OPT_OUTFILE:
                        param->outfile = optarg;
                        break;
                case 'h':
                        help = 1;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }
        if(help){
                RUN(print_termlist_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }

		
        MMALLOC(param->infile,sizeof(char*)* (argc-optind));
        optind++;               /* first command is always termlist; otherwise we not end up here... */
        c = 0;
        while (optind < argc){
                param->infile[c] =  argv[optind++];
                c++;
        }
        param->num_infiles= c;
        ASSERT(c != 0,"No input files!");
        LOG_MSG("Read param.");       
        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}
         

struct parameters* get_insert_param(int argc, char * argv[])
{
        struct parameters* param = NULL;     
        int c, help; 

        RUN(print_program_header(argv,"Retrieves data from OMIM and stores in local sqlite db."));

        
        if(argc == 2){
                RUN(print_insert_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }


        help = 0;
        c = 0;

        RUNP(param = init_param());

        
        while (1){
                static struct option long_options[] ={
                        {"db",required_argument,0,OPT_SQLITE_DB},
                        {"id",required_argument,0, OPT_PATIENT_ID},
                        {"pheno",required_argument,0, OPT_PHENOFILE},
                        {"key",required_argument,0, OPT_OMIM_KEY},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_PHENOFILE:
                        param->phenofile = optarg;
                        break;
                case OPT_OMIM_KEY:
                        param->omimkey = optarg; 
                        break;
                case OPT_PATIENT_ID:
                        param->patient_id = optarg;
                        break;
                case OPT_SQLITE_DB:
                        param->local_sqlite_database_name = optarg;
                        break;
                        
                case 'h':
                        help = 1;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }

        if(help){
                RUN(print_insert_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }
        ASSERT(param->omimkey != NULL,"No omim key.");
        ASSERT(param->local_sqlite_database_name != NULL,"No database.");
        ASSERT(param->patient_id != NULL,"No patient.");
        ASSERT(param->phenofile != NULL,"No phenotype file.");
        
        
        RUN(log_command_line(argc, argv));
        
        LOG_MSG("Read param.");       
        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}


struct parameters* get_panel_param(int argc, char * argv[])
{
        struct parameters* param = NULL;     
        int c, help; 

        RUN(print_program_header(argv,"Retrieves data from OMIM and stores in local sqlite db."));

        
        if(argc == 2){
                RUN(print_panel_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }


        help = 0;
        c = 0;
        RUNP(param = init_param());

        while (1){
                static struct option long_options[] ={
                        {"db",required_argument,0,OPT_SQLITE_DB},
                        {"id",required_argument,0, OPT_PATIENT_ID},
                        {"out",required_argument,0, OPT_OUTFILE},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_OUTFILE:
                        param->outfile = optarg;
                        break;
                case OPT_PATIENT_ID:
                        param->patient_id = optarg;
                        break;
                case OPT_SQLITE_DB:
                        param->local_sqlite_database_name = optarg;
                        break;
                        
                case 'h':
                        help = 1;
                        break;
                case '?':
                        exit(EXIT_FAILURE);
                        break;
                default:
                        WARNING_MSG("unrecognized option.");
                        break;
                }
        }

        if(help){
                RUN(print_panel_help(argc,argv));
                exit(EXIT_SUCCESS); 
        }
        ASSERT(param->local_sqlite_database_name != NULL,"No database.");
        ASSERT(param->patient_id != NULL,"No patient.");
        ASSERT(param->outfile != NULL,"No output filename.");

        
        
        RUN(log_command_line(argc, argv));
        
        LOG_MSG("Read param.");       
        return param;
ERROR:
        LOG_MSG("Something went wrong. Try using the -h option.");
        free_param(param);
        return NULL;
}


struct parameters* init_param(void)     
{
        struct parameters* param = NULL;
        
        MMALLOC(param, sizeof(struct parameters));
        param->omimkey = NULL;
        param->local_sqlite_database_name = NULL;
        param->patient_id = NULL;
        param->phenofile = NULL;
        param->infile = NULL;
        param->num_infiles = 0;
        
        return param;
ERROR:
        free_param(param);
        return NULL;
}


void free_param(struct parameters* param)
{
        //int i; 
        if(param){
                if(param->infile){
                        MFREE(param->infile);
                }
                MFREE(param);
        }
}



int byg_count(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }
	
        int m = (int) strlen(pattern);
        int n = (int) strlen(text);
        int count = 0;
	
        if(m > n){
                return -1;
        }
        int mb = (1 << (m-1));
	
        for (i= 0;i < m;i++){
                T[(int)toupper(pattern[i])] |= (1 << i);
        }
	
        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)toupper(text[i])];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}



