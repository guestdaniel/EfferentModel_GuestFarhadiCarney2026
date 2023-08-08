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
 * Mex function that parses MATLAB inputs and appropriately converts them into formats/types
 * suitable for passing to a simplified wrapper around the full efferent model,
 * `sim_efferent_model`, calls the efferent model, and then handles converting outputs into
 * a format for return to MATLAB.
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
	/* Declare signature for `model_efferent_wrapper`, our entry point to the model */
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
							
	/* Declare pointers we need for handling input arrays (px, cf, randNums) */
	double *randNums_hsrarray, *randNums_lsrarray;
 	int indexcf, indextime;
	int fc, i, lp, l1, l2, pxbins,cfbins;
        mwSize ihcoutsize[2], anrateout_hsrsize[2], anrateout_lsrsize[2], icoutsize[2], gainoutsize[2];
	 mwSize total_num_of_elements_hsr, total_num_of_elements_lsr, index;	
   
	/* Check for proper number of arguments */
	if (nrhs != 23) 
	{
		mexErrMsgTxt("model requires 23 input arguments.");
	}; 

	if (nlhs != 5)  
	{
		mexErrMsgTxt("model requires 5 output argument.");
	};

	/* Get sizes of inputs as necessary */
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

	/* Handle pressure vector (px) */
	double *px_mex = mxGetPr(prhs[0]);
	double *px = (double*) mxCalloc(totalstim, sizeof(double));
	for (int i = 0; i < totalstim; i++) {
		px[i] = px_mex[i];
	}

	/* Handle CF vector (cf) */
	double *cf_mex = mxGetPr(prhs[3]);
	double *cf = (double*) mxCalloc(n_chan, sizeof(double));
	for (int i = 0; i < n_chan; i++) {
		cf[i] = cf_mex[i];
	}
	
	/* Assign pointers to the inputs */
	double *randNums_hsrtmp		= mxGetPr(prhs[1]);
	double *randNums_lsrtmp     = mxGetPr(prhs[2]);

	/* Check with individual input arguments */
	pxbins = mxGetN(prhs[0]);
	
	// px = (double*)mxCalloc(totalstim,sizeof(double)); 
	// for (lp=0; lp<pxbins; lp++)
	// 	px[lp] = pxtmp[lp];
	// double *px = mxGetPr(prhs[0]);
			
	cfbins = mxGetN(prhs[3]);
	
	total_num_of_elements_hsr = mxGetNumberOfElements(prhs[1]);
	total_num_of_elements_lsr = mxGetNumberOfElements(prhs[2]);
	
	 double *randNums_hsr[n_chan];
	 double *randNums_lsr[n_chan];
	 
	randNums_hsrarray = (double*)mxCalloc(cfbins*totalstim,sizeof(double)); 
    randNums_lsrarray = (double*)mxCalloc(cfbins*totalstim,sizeof(double)); 
	
	for (int i = 0; i < n_chan; i++) {
        randNums_hsr[i] = (double*) calloc(totalstim, sizeof(double));
        randNums_lsr[i] = (double*) calloc(totalstim, sizeof(double));
    }

	// randNums_hsrarray column major order 
	for (index=0; index<total_num_of_elements_hsr; index++){

		randNums_hsrarray[index]=*randNums_hsrtmp++;
		randNums_lsrarray[index]=*randNums_lsrtmp++;
		//mexPrintf("%g\n",randNums_hsrarray[index]);
	}
	
	for (index=0; index<total_num_of_elements_hsr; index++){
			indexcf = (int) (fmod(index,cfbins));
	indextime = (int) (index/cfbins);
		
		
	randNums_hsr[indexcf][indextime]=randNums_hsrarray[index];
	randNums_lsr[indexcf][indextime]=randNums_lsrarray[index];
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
   	for (index=0; index<total_num_of_elements_hsr; index++){
	indexcf = (int) (fmod(index,cfbins));
	indextime = (int) (index/cfbins);	
	*ihcouttmp++=ihcout[indexcf][indextime];
	*anrateout_hsrtmp++=anrateout_hsr[indexcf][indextime];
	*anrateout_lsrtmp++=anrateout_lsr[indexcf][indextime];
	*icouttmp++=icout[indexcf][indextime];
	*gaintmp++=gain[indexcf][indextime];
	}

	mxFree(px);
	mxFree(cf);
	for (int i = 0; i < n_chan; i++) {
		free(randNums_hsr[i]);
		free(randNums_lsr[i]);
		free(ihcout[i]);
		free(anrateout_hsr[i]);
		free(anrateout_lsr[i]);
		free(icout[i]);
		free(gain[i]);
	}
	mxFree(randNums_hsrarray);
	mxFree(randNums_lsrarray);
}