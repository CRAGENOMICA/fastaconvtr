/*
 *  common.h
 *  fasta2ms2
 *
 *  Created by Sebastian E. Ramos Onsins on 27/11/2012.
 *
 */

#ifndef COMMON_
#define COMMON_

#ifdef	__cplusplus
extern "C" {
	#endif
	
	#define _CRT_SECURE_NO_DEPRECATE
	
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include <math.h>
    #include <htslib/bgzf.h>

	#define FULL_VERSION "v." VERSION_NUMBER " (" BUILD_NUMBER ")"
	#define FASTA2MS2 "#fastaconvtr" //FULL_VERSION " Sebastian E. Ramos-Onsins.\n"
    #define FASTAMS2VERSION "version " VERSION_NUMBER " (" BUILD_NUMBER ")\n"
    #define FASTAMS2TITLE "Conversion of fasta/ms/tfasta files to tfasta/ms/fasta \nand " \
                          "Construction of WEIGHT files FROM GTF files to be used in mstatspop.\n"
    #define FASTAMS2AUTHORS "Sebastian E. Ramos-Onsins and Ahmed Hafez.\n"

    #define MSP_MAX_FILENAME			(unsigned long) 4096 /**< @brief Maximum Filename Length allowed */
	#define MSP_MAX_GFF_WORDLEN         (unsigned long) 20
	#define MSP_GENETIC_CODETYPE_LEN	(unsigned long) 50	/* e.g. "Nuclear universal" */
	#define MSP_GENCODE_COMBINATIONS    (unsigned long) 64 /* 4^3 */
	#define MSP_GFF_CRITERIA_MSG_LEN    (unsigned long) 20 /* e.g. "MIN" */
    #define MSP_MAX_FILELINE_LEN		(unsigned long) 4096
    #define MSP_MAX_NAME                (unsigned long) 4096
	#define CHI_INTERVAL 10

	#ifndef NULL
		#define NULL	0
	#endif

	typedef struct
	{
		char file_in[ MSP_MAX_FILENAME];
		char file_out[MSP_MAX_FILENAME];

		char format[1];  /*output format: f fasta, t tfasta, m ms format  */
    	char input_format[10]; /*input format F fa or tfa */

		int refasta;
		int tfasta;
		int ploidy,outgroup;
		long int slide,window;

		/* GFF variables */
	    char file_GFF[MSP_MAX_FILENAME];
		int 	gfffiles;//			= 0;
		char 	subset_positions[ MSP_MAX_GFF_WORDLEN ];		
		char 	code_name[ MSP_GENETIC_CODETYPE_LEN ];
		//char 	genetic_code[ MSP_GENCODE_COMBINATIONS ];
		char 	criteria_transcript[ MSP_GFF_CRITERIA_MSG_LEN ];

		char chr_name_all[ MSP_MAX_NAME];
		char file_effsz[MSP_MAX_FILENAME];
		char file_Wcoord[MSP_MAX_FILENAME];
		char file_wps[MSP_MAX_FILENAME];
		char file_masked[MSP_MAX_FILENAME];


		int Physical_length;
    	int include_unknown;

		/*sort samples*/
		int int_total_nsam_order;//=0; /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
		int *sort_nsam ; // =0; /*vector with the position in the file*/

		/* variables defining more data*/
		int npops ;
		/* Number of samples for each population, each element a population */
		int *	vint_perpop_nsam ;
		int printtfasta;
		int argc;
        char file_weights_char[ MSP_MAX_FILENAME ];
        char file_mask_char[ MSP_MAX_FILENAME ];
	} fastaconvtr_args_t;


	int bzprintf(FILE *file_handle, BGZF *z, const char *message, ...);
	FILE * bzopen(const char *filename, const char *opentype, BGZF **z) ;
	int bzgetc(FILE * file_handle, BGZF *z) ;
	
	#ifdef	__cplusplus
	}
#endif

#endif /* COMMON */
