This is the Mex wrapper code for the new auditory subcortical model with efferent gain
control. 

# Installation and usage
1. Make sure you have a C compiler installed and set up with MATLAB. Windows users can
   install "MATLAB Support for MinGW-w64 C/C++ Compiler" from the Add-on Explorer and it
   should be correctly auto-configured.

2. Point MATLAB to this folder, and then run `compile_mex.m` in MATLAB. This will compile
    all of the necessary `.c` files into some `.obj` files and a `.mex*` file in the folder,
    if it is successful! Things may/may not work depending on what versions of MATLAB and
    Mex you have installed, so let me know if there are any issues here and we can
    troubleshoot! 

3. Call `sim_efferent_model` from MATLAB while this folder is on your path to run the
    model... Right now, this is still a work-in-progress function and won't behave exactly
    like old model functions, since a lot of the features and behavior are still in flux...
    key points are:
- The first two arguments are the (1) row-vector sound waveform and (2) the row-vector of
    CFs you want.
- (Matrix-valued) outputs are [ihc, hsr, lsr, ic, gain] , each in the shape (n_chan,
    n_sample) 
- There is currently no padding of inputs/outputs to give the model time to "settle in" or
    to avoid clipping off the end of responses that train beyond the duration of the
    stimulus. For the time being, I would recommend zero-padding your stimulus with some
    10-20 ms of silence if you notice any funny business, but eventually we will automate
    this problem away (so don't worry much about it presently)
- Many more model parameters are exposed than before (e.g., IC parameters, MOC parameters),
     so to simplify function calls, we use the new(ish) MATLAB arguments syntax
    (https://www.mathworks.com/help/matlab/ref/arguments.html). Any arguments other than the
    sound waveform `x` and and the CF vector are passed as key-value pairs like
    `moc_cutoff=0.2` ... 
- Every such model parameter has a default value if you don't explicitly override it, and
    these default values are visible in the code for `sim_efferent_model` (lines 77-96).
    Right now, the defaults are set up to include fractional Gaussian noise and some
    moderate gain control from both the WDR and IC pathways, and the IC pathway is
    configured with some sensible default values (~1-2 ms time constants)

4. Two example simulations are available in `demo.m` 

# Notes

# Changelog
- 7/10/2023: Corrected bug in `sim_efferent_model.m` whereby ffGN was synthesized with
  incorrect parameter values for the `sigma` parameter
- 7/6/2023: Modified `ffGn.m` and `sim_efferent_model.m` to better handle the `noiseType`
  argument. Now, global RNG state should be unaffected by selecting the "frozen" (i.e.,
  `noiseType==0`) ffGn.
- 6/29/2023: New version of model code and `sim_efferent_model.m` that adjusts default
  parameter values for MOC lowpass filter
- 6/21/2023: New version of `sim_efferent_model_mex.c` that fixes bugs related to freeing of
  dynamically allocated memory
- 6/8/2023: New version of `sim_efferent_model_mex.c` that fixes a bug wherein passed values
  for `moc_width_wdr` were ignored and replaced with zeros