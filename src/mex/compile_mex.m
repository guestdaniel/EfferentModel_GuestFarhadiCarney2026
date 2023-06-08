clear all;
mex -c complex.c sfie.c adaptation.c model.c;
mex sim_efferent_model_mex.c complex.obj sfie.obj adaptation.obj model.obj;