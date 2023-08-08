#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <mex.h>
#include <time.h>

#include "complex.hpp"

#define MAXSPIKES 1000000
#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

/**
 *
 * Mex function for `sim_efferent_model`
 *
 * Mex function that parses MATLAB inputs, appropriately converts them into a form 
 * suitable for passing to the model C code, calls the model on the converted inputs, and 
 * then appropriately converts the outputs and passed them back to MATLAB.
 *
 * Note that there are no safety checks built into this function --- all checking should
 * happen in MATLAB before passing to this function.
 *
 * @param nlhs Number of return values (i.e., [n]umber [l]eft [h]and [s]ide values)
 * @param plhs mxArray of pointers to output variables (i.e., [p]ointers to [l]eft [h]and
 * [s]ide values). Obtain a pointer wiht `mexGetPr`.
 * @param nrhs Number of return values (i.e., [n]umber [l]eft [h]and [s]ide values)
 * @param prhs mxArray of pointers to input variables (i.e., [p]ointers to [r]eft [h]and
 * [s]ide values). Obtain a pointer wiht `mexGetPr`.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* 
	 * Declare signature for `model_efferent_wrapper`, our entry point to the model.
	 * This function is defined in `model.c`, and is a simplified version of the full 
	 * `model` function with a more stable signature than the full `model` function (making 
	 * it easier to write the Mex wrapper).
	 */
	void model_efferent_wrapper(
		double *,   // px
		double **,  // randNums_hsr
		double **,  // randNums_lsr
        double *,   // cf
		int,        // n_chan
		double,     // tdres
		int,        // totalstim
		double,     // cohc
        double,     // cihc
		int,        // species
		int,        // powerlaw_mode
		double,     // ic_tau_e
 		double,     // ic_tau_i
        double,     // ic_delay
		double,     // ic_amp
		double,     // ic_inh
		double,     // moc_cutoff
        double,     // moc_beta_wdr
		double,     // moc_offset_wdr
		double,     // moc_beta_ic
        double,     // moc_offset_ic
		double,     // moc_weight_wdr
        double,     // moc_weight_ic
		double,     // moc_width_wdr
		double **,  // ihcout
        double **,  // anrateout_hsr
		double **,  // anrateout_lsr
		double **,  // icout
        double **   // gain
	);
							
	/* todo remove */
 	int indexcf, indextime;
	int fc, i, lp, l1, l2;
        mwSize ihcoutsize[2], anrateout_hsrsize[2], anrateout_lsrsize[2], icoutsize[2], gainoutsize[2];
	 mwSize index;	
   
	/* Check for proper number of arguments */
	if (nrhs != 23) 
	{
		mexErrMsgTxt("model requires 23 input arguments.");
	}; 

	if (nlhs != 5)  
	{
		mexErrMsgTxt("model requires 5 output argument.");
	};

	/* Get sizes of inputs */
	int totalstim = mxGetN(prhs[0]);

	/* De-reference (and, where needed, cast to int) scalar input values */
	int n_chan = (int) *mxGetPr(prhs[4]);
    double tdres = *mxGetPr(prhs[5]);
	double cohc	= *mxGetPr(prhs[6]);
	double cihc	= *mxGetPr(prhs[7]);
	int species = (int) *mxGetPr(prhs[8]);
	double ic_tau_e = *mxGetPr(prhs[9]);
	double ic_tau_i = *mxGetPr(prhs[10]);
	double ic_delay = *mxGetPr(prhs[11]);
	double ic_amp = *mxGetPr(prhs[12]);
	double ic_inh = *mxGetPr(prhs[13]);
	double moc_cutoff = *mxGetPr(prhs[14]);
	double moc_beta_wdr = *mxGetPr(prhs[15]);
	double moc_offset_wdr = *mxGetPr(prhs[16]);
	double moc_beta_ic = *mxGetPr(prhs[17]);
	double moc_offset_ic = *mxGetPr(prhs[18]);
	double moc_weight_wdr = *mxGetPr(prhs[19]);
	double moc_weight_ic = *mxGetPr(prhs[20]);
	double moc_width_wdr = *mxGetPr(prhs[21]);
	int powerlaw_mode = (int) *mxGetPr(prhs[22]);

	/* 
	 * Handle pressure vector by copying data from MATLAB mxArray pointer to dynamically
	 * allocated array in C
	 */
	double *px_mex = mxGetPr(prhs[0]);
	double *px = (double*) calloc(totalstim, sizeof(double));
	for (int i = 0; i < totalstim; i++) {
		px[i] = px_mex[i];
	}

	/* 
	 * Handle CF vector by copying data from MATLAB mxArray pointer to dynamically
	 * allocated array in C
	 */
	double *cf_mex = mxGetPr(prhs[3]);
	double *cf = (double*) calloc(n_chan, sizeof(double));
	for (int i = 0; i < n_chan; i++) {
		cf[i] = cf_mex[i];
	}

	/* 
	 * Handle fGn matrices by copying data from MATLAB mxArray pointer to dynamically
	 * allocated matrix in C. Note that because MATLAB arrays are in column-major order 
	 * while C expects a row-major order, we must reorder the elements in memory before 
	 * passing to C.
	 */
	double *randNums_hsrtmp	= mxGetPr(prhs[1]);
	double *randNums_lsrtmp = mxGetPr(prhs[2]);
	double *randNums_hsr[n_chan];
	double *randNums_lsr[n_chan];
	for (int i = 0; i < n_chan; i++) {
        randNums_hsr[i] = (double*) calloc(totalstim, sizeof(double));
        randNums_lsr[i] = (double*) calloc(totalstim, sizeof(double));
    }
	
	for (int i = 0; i < (totalstim * n_chan); i++) {
		/* Determine Cartesian indices for linear index i given C-style row-major order */
		indexcf = (int) fmod(index, n_chan);
		indextime = (int) (index/n_chan);

		/* Store i-th element in corresponding location in randNums matrices */
		randNums_hsr[indexcf][indextime] = randNums_hsrtmp[i];
		randNums_lsr[indexcf][indextime] = randNums_lsrtmp[i];
	}
		
	/* Create an array for the return argument */
	double *ihcout[n_chan];
	double *anrateout_hsr[n_chan];
	double *anrateout_lsr[n_chan];
	double *icout[n_chan];
	double *gain[n_chan];
    
    ihcoutsize[0] = n_chan;
	ihcoutsize[1] = totalstim;

	plhs[0] = mxCreateNumericArray(2, ihcoutsize, mxDOUBLE_CLASS, mxREAL);    
    plhs[1] = mxCreateNumericArray(2, ihcoutsize, mxDOUBLE_CLASS, mxREAL); 
    plhs[2] = mxCreateNumericArray(2, ihcoutsize, mxDOUBLE_CLASS, mxREAL);  
	plhs[3] = mxCreateNumericArray(2, ihcoutsize, mxDOUBLE_CLASS, mxREAL);    
    plhs[4] = mxCreateNumericArray(2, ihcoutsize, mxDOUBLE_CLASS, mxREAL); 
     	
    
	/* Assign pointers to the outputs */
	for (int i = 0; i < n_chan; i++) {
    	ihcout[i] = (double*) calloc(totalstim, sizeof(double));
		anrateout_hsr[i] = (double*) calloc(totalstim, sizeof(double));
		anrateout_lsr[i] = (double*) calloc(totalstim, sizeof(double));
		icout[i] = (double*) calloc(totalstim, sizeof(double));
		gain[i] = (double*) calloc(totalstim, sizeof(double));
	}
	
	double *ihcouttmp;
	double *anrateout_hsrtmp;
	double *anrateout_lsrtmp;
	double *icouttmp;
	double *gaintmp;
	
	ihcouttmp	  = mxGetPr(plhs[0]);   
	anrateout_hsrtmp= mxGetPr(plhs[1]);  
	anrateout_lsrtmp	= mxGetPr(plhs[2]);  		
	icouttmp	= mxGetPr(plhs[3]);  
	gaintmp= mxGetPr(plhs[4]);  
	
	//mexPrintf("ANmodel: Zilany, Bruce, Ibrahim, and Carney : Auditory Nerve Model\n");
	/* run the model */
	model_efferent_wrapper(
		px, 
		randNums_hsr,
		randNums_lsr,
		cf,
		n_chan,
		tdres,
		totalstim,
		cohc,
		cihc,
		species,
		powerlaw_mode,
		ic_tau_e,
		ic_tau_i,
		ic_delay,
		ic_amp,
		ic_inh,
		moc_cutoff,
		moc_beta_wdr,
		moc_offset_wdr,
		moc_beta_ic,
		moc_offset_ic,
		moc_weight_wdr,
		moc_weight_ic,
		moc_width_wdr,
		ihcout,
		anrateout_hsr,
		anrateout_lsr,
		icout,
		gain
	); 
   	for (index=0; index< (n_chan * totalstim); index++){
	indexcf = (int) (fmod(index,n_chan));
	indextime = (int) (index/n_chan);	
	*ihcouttmp++=ihcout[indexcf][indextime];
	*anrateout_hsrtmp++=anrateout_hsr[indexcf][indextime];
	*anrateout_lsrtmp++=anrateout_lsr[indexcf][indextime];
	*icouttmp++=icout[indexcf][indextime];
	*gaintmp++=gain[indexcf][indextime];
	}

	for (int i = 0; i < n_chan; i++) {
		free(randNums_hsr[i]);
		free(randNums_lsr[i]);
		free(ihcout[i]);
		free(anrateout_hsr[i]);
		free(anrateout_lsr[i]);
		free(icout[i]);
		free(gain[i]);
	}
	free(px);
	free(cf);
	// free(randNums_hsrarray);
	// free(randNums_lsrarray);
}