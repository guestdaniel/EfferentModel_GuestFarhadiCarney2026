# Introduction
This is the Mex wrapper for the ??? (202x) auditory-nerve and midbrain model.

# Installation and usage
1.  Make sure you have a C compiler installed and set up with MATLAB. 
	Windows users can install "MATLAB Support for MinGW-w64 C/C++ Compiler" 
	from the Add-on Explorer.

2.  Change your working directory to the folder containing these files.
    In the command window, `cd("C:\path\to\this\folder")`.

3.  Compile the Mex wrapper running `compile_mex.m` in MATLAB 
	(either by using the "Run" button in the editor window, or by calling 
	`compile_mex` in the Command Window). This will compile all of the 
	necessary `.c` files into a Mex function that MATLAB can use to run
    the model.

4.  Call `sim_efferent_model` from MATLAB while this folder is on your path 
	to run the model. The pattern is:
		`[~, hsr, ~, ~, ~] = sim_efferent_model(x, cf, param1=val1, ..., paramN=valN)`
	- The first two arguments are the (1) row-vector sound waveform and 
	  (2) the row-vector of CFs.
	- (Matrix-valued) outputs are [ihc, hsr, lsr, ic, gain], each in the 
	  shape (n_chan, n_sample) 
	- There is currently no padding of inputs/outputs to give the model 
	  time to "settle in" or to simulate points in time beyond the final 
	  sample of the stimulus. We recommend that you zero-pad your stimulus
	  with some time at the beginning (~10 ms) and at the end (~50 ms) to
      ensure your responses are not affected.
	- All model parameters are exposed as key-value parameter combinations,
	  as in `moc_cutoff=0.2`, passed to `sim_efferent_model` after `x` and
	  `cf` are passed.
	- All model parameters have default values if they're not explicitly
      overridden. To examine default parameters, look in the code of 
      `sim_efferent_model` for the `arguments` block. They key defaults
	  are described briefly here:
		>> Fractional Gaussian noise is included with 2014 parameters
        >> Power-law adaptation is implemented with a new approximation
		   scheme described in Guest and Carney (20xx)
		>> "Rabbit" is the default species (i.e., `species=1`)
		>> The efferent system includes both WDR and IC-driven gain 
           control. See Farhadi et al. (2023) for more information on
           the basic architecture, and Guest et al. (202x) for information 
           on the inclusion of "cross-channel" gain control
		>> The IC model is a band-enhanced neuron with a BMF near 60-100 Hz

4.  Numerous detailed example simulations and plots are available in another
    m file, `demo.m`. 

# Changelog
Changes to the MATLAB/Mex model code are documented here, while changes
to the model code itself are documented in separately in the main model
code in `model_changelog.txt`.

- 11/13/2023, DRG
  Updated model code to use vector-valued COHC/CIHC instead of 
  scalar-valued, so now different channels can have different COHC/CIHC 
  values in multi-channel simulations. Updated this wrapper file to 
  reflect these changes. By default, a vector of ones with the same size
  as the CF vector is passed for both COHC and CIHC. 

- 11/6/2023, DRG
  Updated this document and re-released code with updated model code that
  allows for simulating power-law adaptation with a set of parallel 
  exponentially adapating processes. This approximation is much faster than
  true PLA but retains most of its key properties. (see model changelog for 
  more information). The wrapper was modified to set this scheme as the 
  default implementation (i.e., `powerlaw_mode=2`). Also, added the model 
  changelog to the release folder (previously, it was confusingly not
  included with these files).

- 10/9/2023, DRG
  Updated some documentation in the main wrapper files, converted README
  to .txt and updated, modified Mex compilation file to (hopefully)
  correctly handle different platforms.

- 08/11/2023, DRG
  Cleaned up Mex files, updated documentation in this file, updated this 
  file to include validator functions applied to inputs to ensure that
  only inputs within a sensible range can be provided.

- 07/27/2023, DRG
  Added the fast power-law adaptation approximation based on a parallel
  set of exponential adaptation processes. Enabled by passing
  `powerlaw_mode=2`.

- 07/10/2023, DRG
  Corrected synthesis of fractional Gaussian noise, which was previously
  using inappropriate values for noise variance from the 2018 model code,
  instead of correct values from the 2014 model code.

- 06/29/2023, DRG
  Adjusted MOC lowpass filter cutoff value to 0.64 Hz, which should more
  closely match the value used in the old single-channel efferent model.

- 06/21/2023, DRG
  Fixed a memory leak in the model due to incorrect freeing of
  dynamically allocated memory in the Mex file.

- 06/08/2023, DRG
  New version of Mex file that fixes bug related to incorrectly type
  casting `moc_width_wdr`