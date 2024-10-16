/*
 *  main_fasta2ms2.h
 *  xcode_project
 *
 *  Created by Sebastian E. Ramos Onsins on 27/11/2012.
 *
 */

#ifndef MAIN_FASTA2MS_
#define MAIN_FASTA2MS_

#ifdef	__cplusplus
extern "C" {
	#endif
	
	#include "common.h"
	
    #include "tfasta.h"
	void usage(void);

	int read_fasta(
    tfasta_file *tfasta,
    // FILE *file_input, 
    // SGZip *file_input_gz, 
    FILE *file_output, 
    BGZF *file_output_gz,
    //FILE *file_logerr, 
    //SGZip *file_logerr_gz, 
    // char *input_format,
    int *nsam,
	long int *lenR,
    long int *length_al_real, 
    double *length_al, 
    long int *length_seg,
	long int **matrix_pos, 
    char **matrix_pol, 
    // int ploidy, 
	// int gfffiles, 
    // char *name_fileinputgff,  # args.file_GFF
    // char *subset_positions,  # args.subset_positions
	char *genetic_code, 
    //char *criteria_transcript, # args.criteria_transcript
    // char *format, 
    // int outgroup_presence, # args.outgroup
	double **matrix_sizepos, 
    float *svratio,
    float *summatrix_sizepos,
    long int *nmissing,
    int **mis_pos, 
    float *fnut,
	float *CpG, 
    float *GCs, 
    float *wV, 
    int **svp, 
    float ***pwmatrix_miss,
    /*FILE *file_es, SGZip *file_es_gz,*/ 
    // char *file_in, 
    // char *file_out,
    // int refasta,
    // int tfasta,
	long int *Pp,
    int **CpGp,
    int **Ap,
    int **Cp,
    int **Gp,
    int **Tp,
    int **GCp,
    //int **sort_nsam,
    //int *int_total_nsam_order,
    //int *nsamuser,
    //int npops, 
    double ***sum_sam,
	double ***nsites1_pop,
    double ***nsites2_pop,
    double ***nsites3_pop,
    double ***nsites1_pop_outg,
    double ***nsites2_pop_outg,
    double ***nsites3_pop_outg,
	float *wP,
    float *wPV, 
    FILE *file_ws, 
    BGZF *file_ws_gz,
    long int *wgenes, 
    long int nwindows, 
    // int include_unknown,
    long int *masked_wgenes, 
    long int masked_nwindows,
    char *chr_name,
    unsigned long first,
    unsigned long nscaffolds,
    fastaconvtr_args_t *args,
	int printtfasta);
	
int write_msfile(
	FILE *file_output,
	BGZF *file_output_gz,
	int nsam, 
	long int lenR, 
	long int lenT, 
	double lenP, 
	long int lenS, 
	long int *vector_pos, 
	double *vector_sizepos,
	char *matrix_pol, 
	// long int slide, //  args.slide, 
	// long int window, //  args.window,
	float svratio, 
	float summatrix_sizepos, 
	long int nmissing, 
	int *mis_pos, 
	// char *format, // args.format,
	float *fnut,
	// int Physical_length, // args.Physical_length,  
	float CpG, 
	float GCs, 
	float *wV, 
	long int nV, 
	int *svp, 
	float **pwmatrix_miss, 
	// int tfasta, // args.tfasta
	long int *Pp, 
	int *CpGp, 
	int *Ap, 
	int *Cp, 
	int *Gp, 
	int *Tp, 
	int *GCp, 
	long int *wgenes, 
	long int nwindows, 
	// int *nsamuser,   args.vint_perpop_nsam
	// int npops,  // args.npops
	double **sum_sam,
	double **nsites1_pop, 
	double **nsites2_pop, 
	double **nsites3_pop, 
	double **nsites1_pop_outg, 
	double **nsites2_pop_outg, 
	double **nsites3_pop_outg,
	// int outgroup, // args.outgroup
	fastaconvtr_args_t *args);
	
	/*int read_weights_file(FILE *file_es, SGZip *file_es_gz, FILE *file_output, SGZip *file_output_gz, FILE *file_logerr, SGZip *file_logerr_gz, float **wV, long int **Pp, long int *nV, char *chr_name);*/
	// int read_weights_positions_file(
    //     FILE *file_ws, 
    //     SGZip *file_ws_gz,
    //     FILE *file_output, SGZip *file_output_gz, 
    //     //FILE *file_logerr, SGZip *file_logerr_gz, 
    //     float **wP, float **wPV, float **wV,char *chr_name,unsigned long first);
    // int read_coordinates(
    //     FILE *, SGZip *, 
    //     FILE *, SGZip *, 
    //     //FILE *file_logerr, SGZip *file_logerr_gz,
    //     long int **, long int *,char *chr_name);
    int read_index_file(char *chr_name_all, unsigned long *nscaffolds,char ***chr_name_array,char ***chr_length_array);
	
#ifdef	__cplusplus
}
#endif

#endif /* MAIN_FASTA2MS_ */
