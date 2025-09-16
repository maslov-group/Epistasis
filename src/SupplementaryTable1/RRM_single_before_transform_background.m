%read dataset
data=readtable("RRM2_sample_80_650M_esm2.csv");

%independent x = mut21, dependent y = log(expt_double/expt_mut1)
expt_double_fitness = data.expMutDouble;
expt_mut_1 = data.expMut1;
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;

%independent x = mut12, dependent y = log(expt_double/expt_mut2)
expt_mut_2 = data.expMut2;
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;

R = corrcoef(cat(1,data.mut21,data.mut12), cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2));

R(1,2)

rho = corr(cat(1,data.mut21,data.mut12), cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2),'Type','Spearman');

rho


