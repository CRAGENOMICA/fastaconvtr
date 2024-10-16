#include "main_fasta2ms2.h"
#include "read_fasta.h"
#include "tfasta.h"

#include "log.h"
#include <unistd.h>
	/* CHUNK is the size of every compressed piece of data */
	#define CHUNK                16384  /* bytes 0x4000 */
    #define MAX_FZPRINTF_MESSAGE 0x4000 /* 16384 bytes */

int bzgetc(FILE * file_handle, BGZF *z) {
	// read one character from the file
	int ret = bgzf_getc(z);
	
	if (ret == -1)
	{
		return EOF;
	}
	// // print ret as a character
	// putchar(ret);  // Print the character to stdout
	// // flush the output buffer
	// fflush(stdout);
	return ret;
}

FILE * bzopen(const char * filename, const char *opentype, BGZF **z) {

	FILE *file_handle = fopen(filename, opentype);
	if (file_handle == NULL)
	{
		return NULL;
	}
	int file_descriptor = fileno(file_handle);
	if (file_descriptor < 0)
	{
		return NULL;
	}
	// open BGZF file
	*z = bgzf_dopen(file_descriptor, opentype);
	if (*z == NULL)
	{
		return NULL;
	}
	return file_handle;
}

int bzprintf(FILE *file_handle, BGZF *z, const char *message, ...)
{
	char *buffer = NULL;
	size_t message_len = 0;
	size_t sent_chars = 0;
	int replace_args = 0;
	char buffer_block[CHUNK + 1]; /* +1 for the \x0 */
	va_list args;
	int rest = 0;

	message_len = strlen(message);
	replace_args = (message_len < MAX_FZPRINTF_MESSAGE);

	if (replace_args)
	{
		buffer = (char *)malloc(2 * MAX_FZPRINTF_MESSAGE * sizeof(char));
		va_start(args, message);
		vsnprintf(buffer, 2 * MAX_FZPRINTF_MESSAGE, message, args);
		va_end(args);
	}
	else
	{
		buffer = (char *)message;
	}

	int bytes_written = bgzf_write(z, buffer, strlen(buffer));
	if (bytes_written < 0)
	{
		fprintf(stderr, "Error writing to BGZF file\n");
		bgzf_close(z);
		fclose(file_handle);
		// return bytes_written;
	}
	return bytes_written;
}

int main(int argc, const char * argv[]) {
	//log_set_level(LOG_DEBUG);
	//log_trace("Hello %s", "world");
	//log_debug("Hello %s %d", "world",1);

	// struct to hold all command line arguments
	fastaconvtr_args_t args;
	// set all values to 0
	memset(&args, 0, sizeof(fastaconvtr_args_t));
	// set the number of arguments
	args.argc = argc;

	int x;

	// int ploidy,outgroup;
	// long int slide,window;
	int arg = 0;
	char *f;
	char msformat[10];
    int printtfasta;
	
	int nsam;
	long int lenR,lenT,lenS;
	double lenP;
	long int *vector_pos;
	double *vector_sizepos;
	char *matrix_pol;
	long int nmissing;
	float summatrix_sizepos;
	int *mis_pos; /*vector with the number of Ns in each position*/
	float *fnut;
	// int Physical_length;
    // int include_unknown;
	//int refasta;
	//int tfasta;
    long int *wgenes;
    long int nwindows;
    long int *masked_wgenes;
    long int masked_nwindows;
	float CpG=0.0;
	float GCs=0.0;
	float svratio;
	int *CpGp,*Ap,*Cp,*Gp,*Tp,*GCp; /*CpG,A,C,G,T, GC per position 1/0*/
	
	long int *Pp=0; /*number of position (for effect sizes)*//*not used*/
	long int nV=0; /*number of variants at file ov variant weights (effect sizes)*//*not used*/
	float *wV=0;/*weight at variant (effect sizes)*//*not yet functional although we can recover*/
    float *wP=0;/*weight for each position*/
    float *wPV=0;/*weight for the variant at each position*/
	
	int *svp;/* transitions (1) or transversions (2) per polymorphic position*/
	float **pwmatrix_miss;/*pairwise matrix of differences per position*/
	double **sum_sam;
	double **nsites1_pop;
	double **nsites2_pop;
	double **nsites3_pop;
	double **nsites1_pop_outg;
	double **nsites2_pop_outg;
	double **nsites3_pop_outg;

    char *chr_name;
    // char chr_name_all[ MSP_MAX_NAME];

	// char file_in[ MSP_MAX_FILENAME];
	// char file_out[MSP_MAX_FILENAME];
	// char file_GFF[MSP_MAX_FILENAME];
	// char file_effsz[MSP_MAX_FILENAME];
	// char file_Wcoord[MSP_MAX_FILENAME];
    // char file_wps[MSP_MAX_FILENAME];
    // char file_masked[MSP_MAX_FILENAME];
    
	char file_log[MSP_MAX_FILENAME];
	
	// FILE *file_input	= 	0;
//	FILE *file_es   	= 	0;
	FILE *file_output	=	stdout;
	FILE *file_wcoor    =   0;
    FILE *file_ws   	= 	0;
    FILE *file_msk   	= 	0;

    // FILE *file_logerr   =   stdout; // TODO :: Better use stderr
	FILE *error_log_file = 0; // file_logerr


    // SGZip file_input_gz;
	//  SGZip file_es_gz;
	// SGZip file_output_gz;
    BGZF* file_output_gz;

    BGZF *file_wcoor_gz; // reading
    BGZF *file_ws_gz; // reading
    BGZF *file_msk_gz; // reading
    // SGZip file_logerr_gz;

    // TODO :: handle index where ?
	// struct SGZIndex file_output_gz_index;          /* This is the index for the output gz file. */
    // init_gzindex_structure(&file_output_gz_index); /* IMPORTANT TO INITIALIZE!*/

	/* GFF variables */
	// int 	args.gfffiles			= 0;
	args.gfffiles			= 0;
	// char 	subset_positions[ MSP_MAX_GFF_WORDLEN ];		
	// char 	code_name[ MSP_GENETIC_CODETYPE_LEN ];
	char 	genetic_code[ MSP_GENCODE_COMBINATIONS ];
	// char 	criteria_transcript[ MSP_GFF_CRITERIA_MSG_LEN ];
	
	/*sort samples*/
	// int 
	args.int_total_nsam_order=0; /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
	// int *
	args.sort_nsam=0; /*vector with the position in the file*/
	int sort_index;/*used in the loop*/
	
	/* variables defining more data*/
	// int 
	args.npops = 0;
	/* Number of samples for each population, each element a population */
	// int *	
	args.vint_perpop_nsam = NULL;
	/* Population index */
	int pop_index = 0;
	/* Sum of all nsam of all populations */
	int int_total_nsam = 0;
    
    unsigned long i/*,j,k*/;
    unsigned long nscaffolds;
    char **chr_name_array;
    char **chr_length_array;

    memset( args.chr_name_all, 0, MSP_MAX_NAME);

	memset( args.file_in,  0, MSP_MAX_FILENAME);
	memset( args.file_out, 0, MSP_MAX_FILENAME);
	
	memset( args.file_GFF, 0, MSP_MAX_FILENAME);
	memset( args.file_effsz,  0, MSP_MAX_FILENAME);
    memset( args.file_Wcoord, 0, MSP_MAX_FILENAME);
    memset( args.file_wps, 0, MSP_MAX_FILENAME);
    memset( args.file_masked, 0, MSP_MAX_FILENAME);

    memset( file_log, 0, MSP_MAX_FILENAME);
	

	/*defaults*/
    printtfasta = 1;
	args.format[0] = 't'; /*default tfasta*/
	args.ploidy = 1;
	args.outgroup = 0;
	args.slide = 0;
	args.window = 0;
	args.refasta = 0;
	args.tfasta = 0;
    nwindows = 0;
	wgenes = 0;
    masked_nwindows = 0;
    masked_wgenes = 0;
	args.Physical_length = 1;
    args.include_unknown = 0;
	strcpy( args.criteria_transcript,"long\0");
	strcpy( args.subset_positions,"\0");
    strcpy( args.input_format,"fasta\0");
    strcpy( msformat,"tfasta\0");
    
	#ifdef DEBUG
	log_set_level(LOG_TRACE);
	#else
	log_set_level(LOG_INFO);
	#endif
	// set program name as fastaconvtr
  	const char *program_name = "fastaconvtr";
	log_start(program_name, argc, argv);


    if((args.vint_perpop_nsam = (int *) calloc( (unsigned long)1, sizeof(int) )) == 0) {
        printf("Error allocating memory");
        exit(1);
    }
    /**/
	
	if( argc > 1 ) 
	{
		arg = 1;
		while(arg < argc)
		{
			if( argv[arg][0] != '-' ) 
			{
                if(argv[arg][0] == '>')
                    break;
                // printf(" argument should be -%s ?\n", argv[arg]);
				log_error(" argument should be -%s ?\n", argv[arg]);
                usage();
                exit(1);
			}
			
			switch (argv[arg][1])
			{
                case 'i' : /* input file*/
                    arg++;
                    strcpy( args.file_in, argv[arg] );
                    break;
                case 'o' : /* output file */
                    arg++;
                    strcpy( args.file_out, argv[arg] );

					
                    // if( (file_output = bzopen( args.file_out, "w", &file_output_gz)) == 0) {
                    //     // fprintf(stdout,"\n It is not possible to write in the output file %s\n", file_out);
					// 	log_error("It is not possible to write in the output file %s\n", args.file_out);
                    //     exit(1);
                    // }
                    // file_output_gz.index = &file_output_gz_index;
                    /* Here, after openning the GZ file, the index is assigned to its output gz file. */
                    
                    strcpy(file_log, args.file_out);
                    strcat(file_log,".log");
                    // if( (file_logerr = fzopen( file_log, "w", &file_logerr_gz)) == 0) {
					if ((error_log_file = fopen(file_log, "w")) == 0) {

                        // fprintf(stdout,"\n It is not possible to write the log file %s.", file_log);
						log_error("It is not possible to write the log file %s.", file_log);
                        exit(1);
                    }

                    /*printf("\nOpen log file...");*/
                    // fzprintf(file_logerr,&file_logerr_gz,"\nOpen log file...");
					// TODO :: set log level from command line
      				log_add_fp(error_log_file, LOG_DEBUG);
					log_debug("Open log file...");
                    break;
                case 'F' : /* F fa or tfa */
                    arg++;
                    args.input_format[0] = argv[arg][0];
                    if(args.input_format[0] == 'f') strcpy(args.input_format,"fasta\0");
                    if(args.input_format[0] == 't') strcpy(args.input_format,"tfasta\0");
                    
                    if(args.input_format[0] != 'f' && args.input_format[0] != 't') {
                        //printf("\n Error in -r argument: only the values 'f' (fasta) and 't' (tfasta) are allowed");
						log_error("Error in -F argument: only the values 'f' (fasta) and 't' (tfasta) are allowed");
                        usage();
                        exit(1);
                    }
                    break;
				case 'f' : /*output format: f fasta, t tfasta, m ms format  */
					arg++;
					args.format[0] = argv[arg][0];
					if(args.format[0] == 'm') strcpy(msformat,"ms\0");
					if(args.format[0] == 'x') strcpy(msformat,"x\0");
					if(args.format[0] == 't') strcpy(msformat,"tfasta\0");
                    if(args.format[0] == 'f') {strcpy(msformat,"fasta\0");args.refasta=1;}
					if(args.format[0] == '0') strcpy(msformat,"null\0");
						
					if(args.format[0] != 'm' && args.format[0] != 'x' && args.format[0] != 't' && args.format[0] != 'f' && args.format[0] != '0') {
						// printf("\n Error in -f argument: only the values 'm' (ms), 'f' (fasta), t (tfasta) or 0(nothing) are allowed");
						log_error("Error in -f argument: only the values 'm' (ms), 'f' (fasta), t (tfasta) or 0(nothing) are allowed");
						usage();
						exit(1);
					}
					break;
                case 'p' : /* p Ploidy, 1: haploid, 2: diploid 4: diploid. lowercase is half N*/
                    arg++;
                    args.ploidy = (int)atoi(argv[arg]);
                    if(args.ploidy != 1 && args.ploidy != 2 && args.ploidy != 4) {
                        // printf("\n Error in -p argument: only the values 1 or 2 or 4 (if Ns are counted as lowercase a=AN, c=CN, g=GN, t=TN) are allowed.");
                        log_error("Error in -p argument: only the values 1 or 2 or 4 (if Ns are counted as lowercase a=AN, c=CN, g=GN, t=TN) are allowed.");
						usage();
                        exit(1);
                    }
                    break;
                case 'T' : /* print DNA sequence*/
                    arg++;
                    printtfasta = (int)atoi(argv[arg]);
                    if(printtfasta != 0 && printtfasta != 1) {
                        // printf("\n Error in -T argument: only the values 0 or 1 are allowed.");
						
                        usage();
                        exit(1);
                    }
                    break;
				case 'g': /* g GFF file name, AND more words
						   2nd : synonymous, nonsynonymous, silent or whatever
						   3rd : selected genetic code name or "Other"
						   next 64th's : in case of 'Other', provide 64 IUPAC letters of each codon. 
						   * Check order. 
						   */
					arg++;
					strcpy( args.file_GFF, argv[arg]);					
					arg++;
					strcpy( args.subset_positions, argv[arg] );
					
					args.gfffiles = 1;	
					
					/* Go if known coding option - */
					if( ( strcmp(args.subset_positions,"synonymous")==0 || 
						 strcmp(args.subset_positions,"nonsynonymous")==0 || 
						 strcmp(args.subset_positions,"silent")==0)) 
					{
						arg++;
						strcpy( args.code_name,argv[arg] ); 
						
						if( strcmp(args.code_name, "Nuclear_Universal") == 0 ) 
						{
							genetic_code[0] = 'F';
							genetic_code[1] = 'F';
							genetic_code[2] = 'L';
							genetic_code[3] = 'L';
							genetic_code[4] = 'S';
							genetic_code[5] = 'S';
							genetic_code[6] = 'S';
							genetic_code[7] = 'S';
							genetic_code[8] = 'Y';
							genetic_code[9] = 'Y';
							genetic_code[10] = '*';
							genetic_code[11] = '*';
							genetic_code[12] = 'C';
							genetic_code[13] = 'C';
							genetic_code[14] = '*';
							genetic_code[15] = 'W';
							genetic_code[16] = 'L';
							genetic_code[17] = 'L';
							genetic_code[18] = 'L';
							genetic_code[19] = 'L';
							genetic_code[20] = 'P';
							genetic_code[21] = 'P';
							genetic_code[22] = 'P';
							genetic_code[23] = 'P';
							genetic_code[24] = 'H';
							genetic_code[25] = 'H';
							genetic_code[26] = 'Q';
							genetic_code[27] = 'Q';
							genetic_code[28] = 'R';
							genetic_code[29] = 'R';
							genetic_code[30] = 'R';
							genetic_code[31] = 'R';
							genetic_code[32] = 'I';
							genetic_code[33] = 'I';
							genetic_code[34] = 'I';
							genetic_code[35] = 'M';
							genetic_code[36] = 'T';
							genetic_code[37] = 'T';
							genetic_code[38] = 'T';
							genetic_code[39] = 'T';
							genetic_code[40] = 'N';
							genetic_code[41] = 'N';
							genetic_code[42] = 'K';
							genetic_code[43] = 'K';
							genetic_code[44] = 'S';
							genetic_code[45] = 'S';
							genetic_code[46] = 'R';
							genetic_code[47] = 'R';
							genetic_code[48] = 'V';
							genetic_code[49] = 'V';
							genetic_code[50] = 'V';
							genetic_code[51] = 'V';
							genetic_code[52] = 'A';
							genetic_code[53] = 'A';
							genetic_code[54] = 'A';
							genetic_code[55] = 'A';
							genetic_code[56] = 'D';
							genetic_code[57] = 'D';
							genetic_code[58] = 'E';
							genetic_code[59] = 'E';
							genetic_code[60] = 'G';
							genetic_code[61] = 'G';
							genetic_code[62] = 'G';
							genetic_code[63] = 'G';
						}
						else if(strcmp(args.code_name,"mtDNA_Drosophila")==0) 
						{
							genetic_code[0] = 'F';
							genetic_code[1] = 'F';
							genetic_code[2] = 'L';
							genetic_code[3] = 'L';
							genetic_code[4] = 'S';
							genetic_code[5] = 'S';
							genetic_code[6] = 'S';
							genetic_code[7] = 'S';
							genetic_code[8] = 'Y';
							genetic_code[9] = 'Y';
							genetic_code[10] = '*';
							genetic_code[11] = '*';
							genetic_code[12] = 'C';
							genetic_code[13] = 'C';
							genetic_code[14] = 'W';
							genetic_code[15] = 'W';
							genetic_code[16] = 'L';
							genetic_code[17] = 'L';
							genetic_code[18] = 'L';
							genetic_code[19] = 'L';
							genetic_code[20] = 'P';
							genetic_code[21] = 'P';
							genetic_code[22] = 'P';
							genetic_code[23] = 'P';
							genetic_code[24] = 'H';
							genetic_code[25] = 'H';
							genetic_code[26] = 'Q';
							genetic_code[27] = 'Q';
							genetic_code[28] = 'R';
							genetic_code[29] = 'R';
							genetic_code[30] = 'R';
							genetic_code[31] = 'R';
							genetic_code[32] = 'I';
							genetic_code[33] = 'I';
							genetic_code[34] = 'M';
							genetic_code[35] = 'M';
							genetic_code[36] = 'T';
							genetic_code[37] = 'T';
							genetic_code[38] = 'T';
							genetic_code[39] = 'T';
							genetic_code[40] = 'N';
							genetic_code[41] = 'N';
							genetic_code[42] = 'K';
							genetic_code[43] = 'K';
							genetic_code[44] = 'S';
							genetic_code[45] = 'S';
							genetic_code[46] = 'S';
							genetic_code[47] = 'S';
							genetic_code[48] = 'V';
							genetic_code[49] = 'V';
							genetic_code[50] = 'V';
							genetic_code[51] = 'V';
							genetic_code[52] = 'A';
							genetic_code[53] = 'A';
							genetic_code[54] = 'A';
							genetic_code[55] = 'A';
							genetic_code[56] = 'D';
							genetic_code[57] = 'D';
							genetic_code[58] = 'E';
							genetic_code[59] = 'E';
							genetic_code[60] = 'G';
							genetic_code[61] = 'G';
							genetic_code[62] = 'G';
							genetic_code[63] = 'G';
						}
						else if( strcmp(args.code_name,"mtDNA_Mammals") == 0 ) 
						{
							genetic_code[0] = 'F';
							genetic_code[1] = 'F';
							genetic_code[2] = 'L';
							genetic_code[3] = 'L';
							genetic_code[4] = 'S';
							genetic_code[5] = 'S';
							genetic_code[6] = 'S';
							genetic_code[7] = 'S';
							genetic_code[8] = 'Y';
							genetic_code[9] = 'Y';
							genetic_code[10] = '*';
							genetic_code[11] = '*';
							genetic_code[12] = 'C';
							genetic_code[13] = 'C';
							genetic_code[14] = 'W';
							genetic_code[15] = 'W';
							genetic_code[16] = 'L';
							genetic_code[17] = 'L';
							genetic_code[18] = 'L';
							genetic_code[19] = 'L';
							genetic_code[20] = 'P';
							genetic_code[21] = 'P';
							genetic_code[22] = 'P';
							genetic_code[23] = 'P';
							genetic_code[24] = 'H';
							genetic_code[25] = 'H';
							genetic_code[26] = 'Q';
							genetic_code[27] = 'Q';
							genetic_code[28] = 'R';
							genetic_code[29] = 'R';
							genetic_code[30] = 'R';
							genetic_code[31] = 'R';
							genetic_code[32] = 'I';
							genetic_code[33] = 'I';
							genetic_code[34] = 'M';
							genetic_code[35] = 'M';
							genetic_code[36] = 'T';
							genetic_code[37] = 'T';
							genetic_code[38] = 'T';
							genetic_code[39] = 'T';
							genetic_code[40] = 'N';
							genetic_code[41] = 'N';
							genetic_code[42] = 'K';
							genetic_code[43] = 'K';
							genetic_code[44] = 'S';
							genetic_code[45] = 'S';
							genetic_code[46] = '*';
							genetic_code[47] = '*';
							genetic_code[48] = 'V';
							genetic_code[49] = 'V';
							genetic_code[50] = 'V';
							genetic_code[51] = 'V';
							genetic_code[52] = 'A';
							genetic_code[53] = 'A';
							genetic_code[54] = 'A';
							genetic_code[55] = 'A';
							genetic_code[56] = 'D';
							genetic_code[57] = 'D';
							genetic_code[58] = 'E';
							genetic_code[59] = 'E';
							genetic_code[60] = 'G';
							genetic_code[61] = 'G';
							genetic_code[62] = 'G';
							genetic_code[63] = 'G';
						}
						else if( strcmp(args.code_name,"Other") == 0 ) {
							for(x=0;x<64;x++) {
								arg++;
								if(argv[arg][0] == '-') {
									// printf("\n Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
									log_error("Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
									usage();
									exit(1);
								}
								genetic_code[x] = atoi(argv[arg]);
							}
						}
						else {
							// printf(" %s: Unknown code, sorry", code_name);
							log_error(" %s: Unknown code, sorry", args.code_name);
							exit(1);
						}	
					}
					break;
				case 'c' : /* c Criteria used for analyzing the transcripts */
					/* Basically, max or min */
					arg++;
					strcpy(args.criteria_transcript,argv[arg]);
					if(strcmp( args.criteria_transcript, "max")!=0  && 
					   strcmp( args.criteria_transcript, "min")!=0  && 
					   strcmp( args.criteria_transcript, "first")!=0   && 
					   strcmp( args.criteria_transcript, "long")!=0  ) 
					{
						// printf("\n Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
						log_error("Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
						usage();
						exit(1);
					}
					break;
				case 'G' : /* outgroup present or not */
					arg++;
					args.outgroup = (int)atoi(argv[arg]);
					if(args.outgroup != 0 && args.outgroup != 1) {
						// printf("\n Error in -o argument: only the values 0 or 1 are allowed.");
						log_error("Error in -o argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;
				case 'w' : /* window size (ms format)*/
					arg++;
					args.window = (long int)atol(argv[arg]);
					break;
				case 's' : /* slide size (ms format)*/
					arg++;
					args.slide = (long int)atol(argv[arg]);
					break;					
				case 'P' : /* physical length or effective length (only valid positions) */
					arg++;
					args.Physical_length = (int)atoi(argv[arg]);
					if(args.Physical_length != 0 && args.Physical_length != 1) {
						// printf("\n Error in -P argument: only the values 0 or 1 are allowed.");
						log_error("Error in -P argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;
				/*case 'r' :
					arg++;
					refasta = (int)atoi(argv[arg]);
					if(refasta != 0 && refasta != 1) {
						printf("\n Error in -r argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;*/
				case 'O' : /* O the order of each individual in the original data, Warning!, followed by
							O numbers indicating the order (0 is the first) in case samples are not consecutive. Only for fasta data!
							*/
					arg++;
					args.int_total_nsam_order = atoi(argv[arg]);
					if((args.sort_nsam = (int *) calloc( (unsigned long)args.int_total_nsam_order, sizeof(int) )) == 0) {
						// printf("Error allocating memory");
						log_error("Error allocating memory for sort_nsam");
						exit(1);
					}
					for( sort_index = 0; sort_index < args.int_total_nsam_order; sort_index++ ) 
					{
						arg++;
						args.sort_nsam[ sort_index ] = atoi( argv[arg] );
					}					
					break;
				case 'N' : /* N number of populations, Warning!, followed by
							N numbers indicating the sample size of each population
							*/
                    free(args.vint_perpop_nsam);

                    arg++;
					args.npops = atoi(argv[arg]);
					if((args.vint_perpop_nsam = (int *) calloc( (unsigned long)args.npops, sizeof(int) )) == 0) {
						//printf("Error allocating memory");
						log_error("Error allocating memory for vint_perpop_nsam");
						exit(1);
					}
					int_total_nsam = 0;
					for( pop_index = 0; pop_index < args.npops; pop_index++ ) 
					{
						arg++;
						args.vint_perpop_nsam[ pop_index ] = atoi( argv[arg] );
						int_total_nsam += args.vint_perpop_nsam[ pop_index ];
					}					
					break;
				/*case 't' : 
                 arg++;
					tfasta = (int)atoi(argv[arg]);
					if(tfasta != 0 && tfasta != 1) {
						printf("\n Error in -t argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;*/
                case 'u' : /* -u missing data considered  */
                    arg++;
                    args.include_unknown = (int)atoi(argv[arg]);
                    if(args.include_unknown != 0 && args.include_unknown != 1) {
                        //printf("\n Error in -u argument: only the values 0 or 1 are allowed.");
						log_error("Error in -u argument: only the values 0 or 1 are allowed.");
                        usage();
                        exit(1);
                    }
                    break;
				case 'W' : /* file with the coordinates of each window [init end](overwrite options -w and -s)*/
					arg++;
					strcpy(args.file_Wcoord, argv[arg] );					
					break;
                case 'E' : /*file with the weight for each position */
                    arg++;
                    strcpy(args.file_wps, argv[arg]);
                    break;
				case 'e' : /*file with the effect size of each variant */ /*NOT YET INCLUDED IN THE PROGRAM!!!!!*/
					arg++;
					strcpy(args.file_effsz, argv[arg]);
					break;
                case 'm' : /* file with the coordinates of each window [init end] to be masked by Ns*/
                    arg++;
                    strcpy(args.file_masked, argv[arg] );
                    break;
                case 'n' : /* name of the scaffold to analyze*/
                    arg++;
                    strcpy( args.chr_name_all, argv[arg] );
                    break;
				case 'h' : /* h HELP */
					usage();
					exit(0);
					break;
			}
			arg++;
		}
	}
	else {
		// Print Usage and return an error
		usage();
		exit(1);
	}
    
    /*default*/
    /*
    #TO FORBIDE
    # [-g -c] + [-E]
    # [-P] without [-g -c]
    # [-P] without [-E]
    # [-f tfasta] + [-w -s]
    # [-f fasta]  + [-w -s]
    # [-W] + [-w -s]
     */
    if(args.input_format[0] != 'f' && args.ploidy != 1) {
        // fzprintf(file_logerr,&file_logerr_gz,"\n the option -p 2 is only available with fasta IUPAC input format");
        // printf("\nError: The option -p 2 is only available with fasta IUPAC input format\n");
		log_error("The option -p 2 is only available with fasta IUPAC input format");
        exit(1);
    }
	if(args.file_Wcoord[0]!=0 && (args.slide > 0 && args.window>0)) {
        // fzprintf(file_logerr,&file_logerr_gz,"\n the option -W (coordinates file) is incompatible with definitions of -w (slide) and -s (slide)");
        // printf("\nError: The option -W (coordinates file) is incompatible with definitions of -w (slide) and -s (slide)\n");
		log_error("The option -W (coordinates file) is incompatible with definitions of -w (slide) and -s (slide)");
        exit(1);
	}/*
	if(file_Wcoord[0]!=0 && (format[0] == 'f')) {
        fzprintf(file_output,&file_output_gz,"\n the option -W (coordinates file) is only incompatible with output -f m");
        exit(1);
	}
	if((slide > 0 && window>0) && (format[0] == 't' || format[0] == 'f')) {
        fzprintf(file_output,&file_output_gz,"\n the options -w (slide) and -s (slide) is only compatible with the output -f m");
        exit(1);
	}*/
    if(args.file_wps[0]!=0 && args.file_GFF[0]!=0) {
        // fzprintf(file_logerr,&file_logerr_gz,"\n the option -g (gff file) is incompatible with option -E (weighting file)\n");
		log_error("The option -g (gff file) is incompatible with option -E (weighting file)");
        exit(1);
    }
    if((args.format[0]=='t') && (/**/args.file_Wcoord[0]!=0 || /**/(args.slide > 0 && args.window>0) || args.Physical_length == 0)) {
        //fzprintf(file_logerr,&file_logerr_gz,"\nError: The options -W or -w -s or -P are not effective using the output format -f tfasta\n");
        //printf("\nError: The options -W or -w -s or -P are not effective using the output format -f tfasta\n");
		log_error("The options -W or -w -s or -P are not effective using the output format -f tfasta");
		exit(0);
    }
    if((args.format[0]=='f') && ((args.slide > 0 && args.window>0) || args.Physical_length == 0)) {
        // fzprintf(file_logerr,&file_logerr_gz,"\nError: The options -w -s or -P are not effective using the output format -f fasta\n");
        //printf("\nError: The options -w -s or -P are not effective using the output format -f fasta\n");
        log_error("The options -w -s or -P are not effective using the output format -f fasta");
		exit(0);
    }
    if(args.outgroup==1 && args.npops==0) {
        //fzprintf(file_logerr,&file_logerr_gz,"\nError:  the option -G (outgroup) needs to define option -N: the population samples of at least two pops\n");
        //printf("\nError:  the option -G (outgroup) needs to define option -N: the population samples of at least two pops\n");
        log_error("The option -G (outgroup) needs to define option -N: the population samples of at least two pops");
		exit(1);
    }	
    if( args.format[0] == 'x' && args.npops==0) {
        //fzprintf(file_logerr,&file_logerr_gz,"\n the option -f x require the option -N ");
		log_error("The option -f x require the option -N ");
        exit(1);
    }
    if(args.include_unknown == 1 && args.format[0]=='m') {
        args.include_unknown = 0;
        //fzprintf(file_logerr,&file_logerr_gz,"Warning: The option -u 1 is only allowed in case GFF file is defined and the output format is not ms\n");
        //printf("Warning: The option -u 1 is only allowed in case GFF file is defined and the output format is not ms\n");
		log_warn("The option -u 1 is only allowed in case GFF file is defined and the output format is not ms");
    }
    if(strcmp(args.chr_name_all,"") == 0 &&
       ((args.input_format[0] == 'f' && (args.file_GFF[0] != '\0' || args.file_effsz[0] != '\0' || args.file_Wcoord[0] != '\0' || args.file_wps[0] != '\0' || args.file_masked[0] != '\0')) ||
        (args.input_format[0] == 't' || args.format[0] == 't'))) {
        //fzprintf(file_logerr,&file_logerr_gz,"\nError: the name of the scaffold (option -n) must be defined\n");
        //printf("\nError: the name of the scaffold (option -n) must be defined\n");
        log_error("The name of the scaffold (option -n) must be defined");
		exit(1);
    }
    if(args.vint_perpop_nsam[0]==0) {
        args.npops = 1;
    }
    
    /*Define arrays and vectors*/
    
    if((fnut = (float *)calloc((unsigned long)4,sizeof(float))) == NULL) {
        //fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. get_obsdata.5 \n");
        log_error("Error: memory not reallocated. fnut");
		exit(1);
    }
    /* Definition of a File Stream Buffer, for buffered IO */
    if( (f = (char *)malloc((unsigned long)BUFSIZ*10)) == NULL ) {
        //fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. get_obsdata.4 \n");
        log_error("Error: memory not reallocated for buffered IO f");
		exit(1);
    }

    
    
	
	/* Opening files */
	// if( args.file_in[0] == '\0' ) {
	// 	file_input = stdin;
	// }
	// else {
	// 	if( (file_input = fzopen( args.file_in, "r", &file_input_gz)) == 0) {
	// 		// fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the input file %s\n", file_in);
	// 		log_error("It is not possible to open the input file %s\n", args.file_in);
	// 		exit(1);
	// 	}
	// }
    // setbuf(file_input,f);

   
    if(read_index_file(args.chr_name_all,&nscaffolds,&chr_name_array,&chr_length_array)) {
        // printf("\nError reading the scaffold names file %s\n",chr_name_all);
		log_error("Error reading the scaffold names file %s\n",args.chr_name_all);

        exit(1);
    }
    /*separate all values of the list chr_name_all in chr_name_array: */
    /* Only do the list if input and output is tfa*//*
    nscaffolds = 1;
    if(format[0] == 't' && input_format[0] == 't' ) {
        chr_name_array = (char **)calloc(nscaffolds,sizeof(char *));
        chr_name_array[0] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
        j=0;
        while(chr_name_all[j] != '\0') {
            k=0;
            while(chr_name_all[j] != ',' && chr_name_all[j] != '\0' && j < MSP_MAX_NAME) {
                chr_name_array[nscaffolds-1][k] = chr_name_all[j];
                j++; k++;
            }
            if(chr_name_all[j] == ',') {
                nscaffolds += 1;
                chr_name_array = (char **)realloc(chr_name_array,nscaffolds*sizeof(char *));
                chr_name_array[nscaffolds-1] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
                j++;
            }
       }
    }
    else {
        chr_name_array = (char **)calloc(1,sizeof(char *));
        chr_name_array[0] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
        strcpy(chr_name_array[0],chr_name_all);
    }
    */
     // data structure to hold tfasta weights file and index
	wtfasta_file *wtfasta;
	wtfasta = NULL;
	/*open the file for weigth for positions, if included*/
	if (args.file_wps[0] == '\0')
	{
		file_ws = 0;
	}
	else
	{
		if ((file_ws = bzopen(args.file_wps, "r", &file_ws_gz)) == 0)
		{
			// fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the weighting file %s\n", file_wps);
			log_error("It is not possible to open the weighting file %s\n", args.file_wps);
			exit(1);
		}

		log_info("Using tfasta weights file %s", args.file_wps);
		log_info("Validating tfasta weights file %s", args.file_wps);
		// validate weights file exists and is readable
		if (access(args.file_wps, R_OK) != 0)
		{
			log_error("Error: the tfasta weights file %s does not exist or is not readable.", args.file_wps);
			exit(1);
		}

		if ((wtfasta = (wtfasta_file *)malloc(sizeof(wtfasta_file))) == NULL)
		{
			// fprintf(file_logerr, "\nError: memory not reallocated. mstatspop.c.00 \n");
			log_fatal("Error: can not allocate memory for wtfasta_file structure.");
			exit(1);
		}

		// validate weights file is in TFAv2.0 format
		// memset the structure to 0
		memset(wtfasta, 0, sizeof(wtfasta_file));
		int ret_status = init_wtfasta_file(wtfasta, args.file_wps);
		if (ret_status == TFA_ERROR)
		{
			log_error("Can not initialize tfasta weights file %s", args.file_wps);
			exit(1);
		}
		if (ret_status == TFA_INVALID_FORMAT)
		{
			log_error("The tfasta weights file %s is not in wTFAv2.0 format", args.file_wps);
			exit(1);
		}
		if (ret_status == TFA_INVALID_INDEX)
		{
			log_error("The index file for %s is invalid", args.file_wps);
			exit(1);
		}
		log_info("tfasta weights file %s is valid", args.file_wps);
	}

	tfasta_file *tfasta;
	tfasta = NULL;
	//if (args.format[0] == 't')
	// if input format is tfasta
	if (args.input_format[0] == 't')
	{
		

		// read tfasta file
		// validate input tfasta file and index file

		log_info("Validating tfasta file %s", args.file_in);

		// validate tfasta file exists and is readable and in TFAv2.0 format
		if (access(args.file_in, R_OK) != 0)
		{
			// fprintf(file_logerr, "\nError: the file %s does not exist or is not readable.\n", file_in);
			log_error("Error: the file %s does not exist or is not readable.", args.file_in);
			exit(1);
		}

		if ((tfasta = (tfasta_file *)malloc(sizeof(tfasta_file))) == NULL)
		{
			
			log_fatal("Error: can not allocate memory for tfasta_file structure.");
			exit(1);
		}

		// validate weights file is in TFAv2.0 format
		// memset the structure to 0 
		memset(tfasta, 0, sizeof(tfasta_file));

		
		int ret_status = init_tfasta_file(tfasta, args.file_in);
		if (ret_status == TFA_ERROR)
		{
			log_error("Can not initialize tfasta file %s", args.file_in);
			exit(1);
		}
		if (ret_status == TFA_INVALID_FORMAT)
		{
			log_error("The file %s is not in TFAv2.0 format", args.file_in);
			exit(1);
		}
		if (ret_status == TFA_INVALID_INDEX)
		{
			log_error("The index file for %s is invalid", args.file_in);
			exit(1);
		}
	}

	// prepare the output file
	if (args.format[0] == 't') // tfasta is compressed
	{
		file_output = bzopen(args.file_out, "wb", &file_output_gz); // fopen(args.file_out, "wb");
		if (file_output == NULL)
		{
			log_error("It is not possible to write in the output file %s\n", args.file_out);
			exit(1);
		}
		
	}
	else {
		// use wu for uncompressed output
		file_output = bzopen(args.file_out, "wu", &file_output_gz); // fopen(args.file_out, "wb");

		if (file_output == NULL)
		{
			log_error("It is not possible to write in the output file %s\n", args.file_out);
			exit(1);
		}
	}

	/*do a loop using each value of chr_names_all into chr_name.*/
    for(i=0;i<nscaffolds;i++) {
        /*open the file for effect sizes, if included*/
        /*
        if( file_effsz[0] == '\0') {
            file_es = 0;
            wV = 0;
        }
        else {
            if( (file_es = fzopen( file_effsz, "r", &file_es_gz)) == 0) {
                fzprintf(file_output,&file_output_gz,"\n It is not possible to open the effect sizes file %s\n", file_effsz);
                exit(1);
            }
        }
        */
        /* Opening coordinates file */
        if( args.file_Wcoord[0] == '\0' ) {
            file_wcoor = 0;
            nwindows = 0;
        }
        else {
			
            if( (file_wcoor = bzopen( args.file_Wcoord, "r", &file_wcoor_gz)) == 0) {
                // fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the coordinates file %s\n", file_Wcoord);
				log_error("It is not possible to open the coordinates file %s", args.file_Wcoord);
                exit(1);
            }
        }
        /* Opening mask coordinates file */
        if( args.file_masked[0] == '\0' ) {
            file_msk = 0;
            masked_nwindows = 0;
        }
        else {
            if( (file_msk = bzopen( args.file_masked, "r", &file_msk_gz)) == 0) {
                // fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the masked coordinates file %s\n", file_masked);
				log_error("It is not possible to open the masked coordinates file %s", args.file_masked);
                exit(1);
            }
        }
 
        chr_name = chr_name_array[i];
        
        /*read the file for weigth for positions, if included*/
        if( args.file_wps[0] != '\0') {
			long int w_size = 0;
            if(read_weights_positions_file(
				wtfasta,
				&wP,&wPV,&wV,chr_name,i,&w_size) == 0) {
                //fzprintf(file_logerr,&file_logerr_gz,"Error processing weighting file %s\n", file_wps);
				log_error("Error processing weighting file %s", args.file_wps);
                exit(1);
            }
        }
        /*read the file for effect sizes, if included*/
        /*
        if( file_effsz[0] != '\0') {
            if(read_weights_file(file_es,&file_es_gz,file_output,&file_output_gz, file_logerr, &file_logerr_gz,&wV,&Pp,&nV,chr_name) == 0) {
                fzprintf(file_output,&file_output_gz,"Error processing effect sizes file %s\n", file_effsz);
                exit(1);
            }
        }
        fzclose(file_es, &file_es_gz);
         */
        /* read coordinates file */
        if( args.file_Wcoord[0] != '\0' ) {
            if(read_coordinates(
				file_wcoor,file_wcoor_gz,
				file_output,file_output_gz, 
				//file_logerr, &file_logerr_gz,
				&wgenes,&nwindows,chr_name) == 0) {
                //fzprintf(file_logerr,&file_logerr_gz,"Error processing coordinates file %s\n", file_Wcoord);
				log_error("Error processing coordinates file %s", args.file_Wcoord);
                exit(1);
            }
            args.window = -1;
            args.slide = -1;
            // fzclose(file_wcoor, &file_wcoor_gz);
			bgzf_close(file_wcoor_gz);
        }
        /* read mask coordinates file */
        if( args.file_masked[0] != '\0' ) {
            if(read_coordinates(
				file_msk,file_msk_gz,
				file_output,
				file_output_gz,
				//file_logerr, &file_logerr_gz,
				&masked_wgenes,
				&masked_nwindows,
				chr_name) == 0) {
                // fzprintf(file_logerr,&file_logerr_gz,"Error processing masked coordinates file %s\n", file_masked);
				log_error("Error processing masked coordinates file %s", args.file_masked);
                exit(1);
            }
            // fzclose(file_msk, &file_msk_gz);
			bgzf_close(file_msk_gz);
        }

		if (args.format[0] == 't')
		{

			char *v2_header = "##fileformat=TFAv2.0\n";
			bzprintf(file_output, file_output_gz, v2_header);
			args.tfasta = 1;
		}

		/*print all the argv. Header*/
// TODO :: This make it hard to make an automated test
// TODO :: can be disabled for debuging and testing and redirect that to the log file instead
#ifndef DEBUG
		bzprintf(file_output, &file_output_gz, "#fastaconvtr ");
		for (x = 1; x < arg; x++)
		{
			bzprintf(file_output, &file_output_gz, "%s ", argv[x]);
		}
		bzprintf(file_output, &file_output_gz, "\n");
#else
		bzprintf(file_output, file_output_gz, "#fastaconvtr ");
		for (x = 1; x < arg; x++)
		{
			bzprintf(file_output, file_output_gz, "%s ", argv[x]);
		}
		bzprintf(file_output, file_output_gz, "\n");
#endif
		if (read_fasta(
				tfasta,
				// file_input,
				// 		&file_input_gz,
				file_output,
				file_output_gz,
				// file_logerr,
				// &file_logerr_gz,
				// input_format,
				&nsam,
				&lenR,
				&lenT,
				&lenP,
				&lenS,
				&vector_pos,
				&matrix_pol,
				// ploidy,
				// gfffiles,
				// file_GFF,
				// subset_positions,
				genetic_code,
				// criteria_transcript,
				// format,
				// outgroup,
				&vector_sizepos,
				&svratio,
				&summatrix_sizepos,
				&nmissing,
				&mis_pos,
				fnut,
				&CpG,
				&GCs,
				wV,
				&svp,
				&pwmatrix_miss,
				/*file_es,&file_es_gz,*/
				// file_in,
				// file_out,
				// refasta,
				// tfasta,
				Pp,
				&CpGp,
				&Ap,
				&Cp,
				&Gp,
				&Tp,
				&GCp,
				// &sort_nsam,
				// &int_total_nsam_order,
				// vint_perpop_nsam,
				// npops,
				&sum_sam,
				&nsites1_pop,
				&nsites2_pop,
				&nsites3_pop,
				&nsites1_pop_outg,
				&nsites2_pop_outg,
				&nsites3_pop_outg,
				wP,
				wPV,
				file_ws,
				file_ws_gz,
				wgenes,
				nwindows,
				// include_unknown,
				masked_wgenes,
				masked_nwindows,
				chr_name,
				i,
				nscaffolds,
				&args) == 0)
		{
			// fzprintf(file_logerr,&file_logerr_gz,"Error processing input data.\n");
			log_error("Error processing input data.");
			exit(1);
		}

		if(args.format[0] == 'm') {
            if(args.slide == 0 && args.window == 0) {
                args.slide = lenR;
                args.window = lenR;
            }
			if (write_msfile(file_output, 
							file_output_gz,
							 // file_logerr,&file_logerr_gz,
							 nsam, 
							 lenR, 
							 lenT, 
							 lenP, 
							 lenS, 
							 vector_pos, 
							 vector_sizepos, 
							 matrix_pol,
							 // args.slide, 
							 // args.window, 
							 svratio, 
							 summatrix_sizepos, 
							 nmissing, 
							 mis_pos, 
							 // args.format,
							 fnut, 
							 // args.Physical_length, 
							 CpG, 
							 GCs, 
							 wV, 
							 nV, 
							 svp, 
							 pwmatrix_miss, 
							 // args.tfasta,
							 Pp, 
							 CpGp, 
							 Ap, 
							 Cp, 
							 Gp, 
							 Tp, 
							 GCp, 
							 wgenes, 
							 nwindows, 
							 // args.vint_perpop_nsam, 
							 // args.npops, 
							 sum_sam,
							 nsites1_pop, 
							 nsites2_pop, 
							 nsites3_pop,
							 nsites1_pop_outg, 
							 nsites2_pop_outg, 
							 nsites3_pop_outg,
							 // args.outgroup,
							 &args) == 0)
			{
				// fzprintf(file_logerr,&file_logerr_gz,"Error printing %s ms file.\n",msformat);
				log_error("Error printing %s ms file.", args.format);
				exit(1);
			}
			/*
            if(format[0] == 'e' || format[0] == 'x') free(mis_pos);
            free(fnut);
            free(sort_nsam);
            */
        }
        /*free all arrays*/
        if(file_wcoor) free(wgenes);
        if(file_msk) free(masked_wgenes);
        if(file_ws/* || file_es*/) free(wP);
        if(file_ws/* || file_es*/) free(wPV);
        if(file_ws) free(wV);
//        if(file_es) free(Pp);
        free(vector_sizepos);
        if(args.format[0] == 'm') {
            free(matrix_pol);
            free(vector_pos);
            free(mis_pos);
            free(CpGp);
            free(Tp);
            free(Ap);
            free(Cp);
            free(Gp);
            free(GCp);
            free(svp);
            for(x=0;x<lenR;x++) {
                free(nsites1_pop[x]);
                free(nsites2_pop[x]);
                free(nsites3_pop[x]);
                free(nsites1_pop_outg[x]);
                free(nsites2_pop_outg[x]);
                free(nsites3_pop_outg[x]);
                free(sum_sam[x]);
                free(pwmatrix_miss[x]);
            }
            free(nsites1_pop);
            free(nsites2_pop);
            free(nsites3_pop);
            free(nsites1_pop_outg);
            free(nsites2_pop_outg);
            free(nsites3_pop_outg);
            free(sum_sam);
            free(pwmatrix_miss);
        }
    }
    
	/*!added. Here, we ensure that all files are closed before exiting */
	// fzclose(file_output, &file_output_gz);
	// Close the BGZF stream and the original file stream
    bgzf_close(file_output_gz);
	// create an index file for the output file
	if (args.tfasta) {
		log_info("Creating index file for the output file %s", args.file_out);
		BGZF *output_fp = bgzf_open(args.file_out, "r");
		tbx_t *tbx;
		

		tbx = tbx_index(output_fp, 0, &tfasta_conf);
		bgzf_close(output_fp);
		if (!tbx)
			 {
				log_error("Failed to create index file for the output file %s", args.file_out);
				exit(1);

			 }
		int ret = hts_idx_save_as(tbx->idx, args.file_out, NULL, HTS_FMT_TBI);
		tbx_destroy(tbx);

	}
    fclose(file_output);



	// fzclose(file_input, &file_input_gz);
    
	if(tfasta) {
		close_tfasta_file(tfasta);
		free(tfasta);
	}


	if(file_ws) 
		//fzclose(file_ws, &file_ws_gz);
		bgzf_close(file_ws_gz);

    //fzprintf(file_logerr,&file_logerr_gz,"\nProgram Ended\n");
	log_info("Program Ended");
    // fzclose(file_logerr, &file_logerr_gz);

    for(i=0;i<nscaffolds;i++) {
        free(chr_name_array[i]);
    }
    free(chr_name_array);
    free(args.vint_perpop_nsam);
    free(fnut);
    free(args.sort_nsam);
    // free(f);

	if (error_log_file)
    	fclose(error_log_file);
    /*
     * Test the just created GZ and INDEX files:
     * -----------------------------------------
    {
		FILE *h = 0;
		SGZip gz;

		h = fzopen("../Examples/output.tfa.gz", "r", &gz);

		if (h != NULL) {

			struct SGZIndex idx;
			load_index_from_file(gz.index_file_name, &idx);

		  long int row_num = 0;
			while (fzseek(h, &gz, &idx, NULL, &row_num, false) == GZ_OK) {

				char ch = ' ';
				while((!fzeof(h, &gz)) && (ch != '\n') && (ch != '\x0')) {
					ch = fzgetc(h, &gz);
					printf("%c", ch);
				}

				row_num++;
			}

			unload_all_index_positions(&idx);

			fzclose(h, &gz);
		}
		exit(0);
    }
    */


    return 0;
}
void usage(void) 
{
	printf(FASTA2MS2);
	printf("\nFlags:\n");
    printf("      -F [input format file: f (fasta), t (tfasta)]\n");/*fasta or tfasta formats are only available*/
    printf("      -i [path and name of the input file (text or gz indexed)]\n");
    printf("      -f [output format file: t (tfasta), f (fasta), m (ms), 0(nothing)]\n");
    printf("      -o [path and name of the output file (include extension .gz except ms files)]\n");
    printf("      -n [name of the file containing the name(s) of scaffold(s) and their length (separated by a tab), one per line (ex. fai file)]\n");
    printf("   OPTIONAL PARAMETERS:\n");
    printf("      -h [help and exit]\n");
    printf("      -P [define window lengths in 'physical' positions (1) or in 'effective' positions (0)]. DEFAULT: 1\n");
    printf("      -O [#_nsam] [Reorder samples: number order of first sample, number 0 is the first sample] [second sample] ...etc.\n");
    printf("      -W [for ms and fasta outputs, file with the coordinates of each window: (one header plus nlines with init end]\n");/*new!*/
    printf("      -N [#_pops] [#samples_pop1] ... [#samples_popN] (necessary in case to indicate the outgroup population)\n");
    printf("      -G [outgroup included (1) or not (0), last population (1/0)]. DEFAULT: 0\n");
    printf("      -u [Missing counted (1) or not (0) in weights given GFF annotation]. DEFAULT: 0\n");
    printf("      -m [masking regions: file indicating the start and the end of regions to be masked by Ns]\n");
    printf("     Outputing ms format:\n");
    printf("      -w [window size]. DEFAULT: Total_length\n");
    printf("      -s [slide size]. DEFAULT: Total_length\n");
    printf("     Inputing fasta format:\n");
    printf("      -p [if fasta input,\n");
    printf("             haplotype: 1 (single sequence)\n");
    printf("             genotype:  2 or 4 (two diploid mixed sequences in IUPAC format. If p=4 lowercase will be considered as one haplotype missing!). DEFAULT: 1\n");
    printf("     Annotation file and weight options:\n");
    printf("      -g [GFF_file]\n");
    printf("         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (whatever annotated)]\n");
    printf("         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]\n");
    printf("         [if 'Other', introduce the single letter code for the 64 triplets in the order UUU UUC UUA UUG ... etc.]\n");
    printf("      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT: long\n");
    printf("      -E [instead -g & -c, input file with weights for positions: include three columns with a header, first the physical positions (1...end), second the weight for positions and third a boolean weight for the variant (eg. syn variant but nsyn position)]\n");
    printf("      -T [in case define -g and -c and output is TFasta, option to print (1) or not (0) the DNA sequence]. DEFAULT: 1 (print)\n");
    /*printf("      -r [rewrite the fasta file for selected region (not valid for silent/syn/nsyn) (1/0)]\n");*/
    /*printf("      -t [rewrite the fasta transposed file including the weight of each position and variant, if available) (1/0)]\n");*//*new!*/
    /*printf("\     -e [input file with effect sizes for variants: include two columns with a header, first the physical positions and second the weight]\n");*/
	printf("\n");
}
