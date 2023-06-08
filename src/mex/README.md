This is the Mex wrapper code for the new auditory subcortical model with efferent gain
control. 

# Installation and usage
1. Point MATLAB to this folder, and then run `compile_mex.m` in MATLAB. This will compile
    all of the necessary `.c` files into `.obj` files and build a `.mex*` file in the
    folder, if it is successful! Things may/may not work depending on what versions of
    MATLAB and Mex you have installed, so let me know if there are any issues here and we
    can troubleshoot! 

2. Call `sim_efferent_model` from MATLAB while this folder is on your path to run the model... Right now,
    this is still a work-in-progress function and won't behave exactly like old model
    functions, since a lot of the features and behavior are still in flux... key points are:
- The first two arguments are the (1) row-vector sound waveform and (2) the row-vector of
    CFs you want.
- (Matrix-valued) outputs are [ihc, hsr, lsr, ic, gain] , each in the shape (n_chan,
    n_sample) 
- There is currently no padding of inputs/outputs to give the model time to "settle in" or
    to avoid clipping off the end of responses that train beyond the duration of the stimulus. For the time being, I would recommend zero-padding your stimulus with some 10-20 ms of silence if you notice any funny business, but eventually we will automate this problem away (so don't worry much about it presently)
- Many more model parameters are exposed than before (e.g., IC parameters, MOC parameters), so to avoid having a big long function
    call with 10+ positional arguments to remember the order of, we use the new(ish) MATLAB
    arguments syntax (https://www.mathworks.com/help/matlab/ref/arguments.html), so that any
    arguments other than x and cf are passed as key-value pairs like `moc_cutoff=0.2` ... 
- Every parameter has a default value if you don't explicitly override it, and these default
    values are visible in the code for `sim_efferent_model` (lines 77-96). Right now, the
    defaults are set up to include fractional Gaussian noise and some moderate gain control
    from both the WDR and IC pathways, and the IC pathway is configured with some sensible
    default values (~1-2 ms time constants)

3. Two example simulations are available in `demo.m` 
# Notes

# Changelog