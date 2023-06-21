% Delete existing intermediate files and compilex Mex file
delete *.obj;
delete *.mex*;

% Compile individual C files
mex -c complex.c sfie.c adaptation.c model.c;

% Compile Mex wrapper
mex sim_efferent_model_mex.c complex.obj sfie.obj adaptation.obj model.obj;