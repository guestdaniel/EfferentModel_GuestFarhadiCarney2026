%% compile_mex.m
% Author: Daniel R. Guest
% Date: 7/26/2023
% 
% Notes:
% - Unclear whether we need to include `model.h` (header file) in
%   compilation steps

% Delete existing intermediate files and compilex Mex file
delete *.obj;
delete *.mex*;

% Compile individual C files
mex -c complex.c sfie.c adaptation.c model.c;  

% Compile Mex wrapper
mex sim_efferent_model_mex.c complex.obj sfie.obj adaptation.obj model.obj;