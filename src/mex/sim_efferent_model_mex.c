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

/* This function is the MEX "wrapper", to pass the input and output variables between the .dll or .mexglx file and Matlab */						
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Declare variables for this function */
	double *px, *cf, tdres, cohc, cihc, ic_tau_e, ic_tau_i, ic_delay, ic_amp, 
	        ic_inh, moc_cutoff, moc_beta_wdr, moc_offset_wdr, moc_beta_ic,
	        moc_offset_ic, moc_weight_wdr, moc_weight_ic, moc_len_integ;
	double *randNums_hsrarray;
	double *randNums_lsrarray;
	int indexcf, indextime;
	int fc, i, lp, l1, l2, n_chan, totalstim, pxbins,species,cfbins;
        mwSize ihcoutsize[2], anrateout_hsrsize[2], anrateout_lsrsize[2], icoutsize[2], gainoutsize[2];
	 mwSize total_num_of_elements_hsr, total_num_of_elements_lsr, index;	
	double *randNums_hsrtmp, *randNums_lsrtmp, *pxtmp, *cftmp, *tdrestmp, *nchantmp, *cohctmp, *cihctmp,  *speciestmp, *ic_tau_etmp, *ic_tau_itmp;
    double *ic_delaytmp, *ic_amptmp, *ic_inhtmp, *moc_cutofftmp, *moc_beta_wdrtmp, *moc_offset_wdrtmp, *moc_beta_ictmp, *moc_offset_ictmp, *moc_weight_wdrtmp, *moc_weight_ictmp, *moc_len_integtmp;      
    
	/* Declare the signature for the C function we call to actually 
	   implement the model --- TODO: move to separate file */
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
							
	/* Check for proper number of arguments */
	if (nrhs != 23) 
	{
		mexErrMsgTxt("model requires 23 input arguments.");
	}; 

	if (nlhs != 5)  
	{
		mexErrMsgTxt("model requires 5 output argument.");
	};
	
	/* Assign pointers to the inputs */
	pxtmp				= mxGetPr(prhs[0]);
	randNums_hsrtmp		= mxGetPr(prhs[1]);
	randNums_lsrtmp     = mxGetPr(prhs[2]);
	cftmp  				= mxGetPr(prhs[3]);
	nchantmp     		= mxGetPr(prhs[4]);
    tdrestmp            = mxGetPr(prhs[5]);
	cohctmp      		= mxGetPr(prhs[6]);
	cihctmp    			= mxGetPr(prhs[7]);
	speciestmp	        = mxGetPr(prhs[8]);
	ic_tau_etmp    		= mxGetPr(prhs[9]);
	ic_tau_itmp 	    = mxGetPr(prhs[10]);
	ic_delaytmp         = mxGetPr(prhs[11]);
	ic_amptmp    		= mxGetPr(prhs[12]);
	ic_inhtmp		    = mxGetPr(prhs[13]);
	moc_cutofftmp		= mxGetPr(prhs[14]);
	moc_beta_wdrtmp     = mxGetPr(prhs[15]);
	moc_offset_wdrtmp   = mxGetPr(prhs[16]);
	moc_beta_ictmp		= mxGetPr(prhs[17]);
	moc_offset_ictmp    = mxGetPr(prhs[18]);
	moc_weight_wdrtmp   = mxGetPr(prhs[19]);
	moc_weight_ictmp    = mxGetPr(prhs[20]);
	moc_len_integtmp    = mxGetPr(prhs[21]);
	int powerlaw_mode = mxGetPr(prhs[22])[0];

	/* Check with individual input arguments */
	pxbins = mxGetN(prhs[0]);
	totalstim = pxbins;  
	
	px = (double*)mxCalloc(totalstim,sizeof(double)); 
	for (lp=0; lp<pxbins; lp++)
		px[lp] = pxtmp[lp];
			
	if (pxbins==1)
		mexErrMsgTxt("px must be a row vector\n");
	
	n_chan = (int)nchantmp[0];
	
	if (nchantmp[0]!=n_chan)
		mexErrMsgTxt("n_chan must an integer.\n");
	if (n_chan<1)
		mexErrMsgTxt("n_chan must be greater than 0.\n");
	
	cfbins = mxGetN(prhs[3]);
	
	cf = (double*)mxCalloc(cfbins,sizeof(double)); 
	for (fc=0; fc <n_chan; fc++)
		cf[fc] = cftmp[fc];
	
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
		
	cihc = cihctmp[0];
	cohc = cohctmp[0];
	
	species = (int) speciestmp[0];
	if (speciestmp[0]!=species)
		mexErrMsgTxt("species must an integer.\n");
	if (species<1 || species>3)
		mexErrMsgTxt("Species must be 1 for cat, or 2 or 3 for human.\n");
	
	ic_tau_e = ic_tau_etmp[0];
	ic_tau_i = ic_tau_itmp[0];
	ic_delay = ic_delaytmp[0];
	ic_amp   = ic_amptmp[0];
	ic_inh   = ic_inhtmp[0];
	moc_cutoff   = moc_cutofftmp[0];
	moc_beta_wdr = moc_beta_wdrtmp[0];
	moc_offset_wdr = moc_offset_wdrtmp[0];
	moc_beta_ic = moc_beta_ictmp[0];
	moc_offset_ic = moc_offset_ictmp[0];
	moc_weight_wdr = moc_weight_wdrtmp[0];
	moc_weight_ic = moc_weight_ictmp[0];
	moc_len_integ = moc_len_integtmp[0];
	
    tdres = tdrestmp[0];	

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
		moc_len_integ,
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