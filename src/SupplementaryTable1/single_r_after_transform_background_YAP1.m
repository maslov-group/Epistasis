%read YAP1 dataset
data=readtable("20_esm2_650_unique_single_mutations_YAP1.csv");
data1=readtable("80_esm2_650_unique_single_mutations_YAP1.csv");
ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = mut21, dependent y = log(expt_double/expt_mut1)
expt_double_fitness = data.log_fitness;
expt_mut_1 = data.Mut1_expt_fitness;
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;

%independent x = mut12, dependent y = log(expt_double/expt_mut2)
expt_mut_2 = data.Mut2_expt_fitness;
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;

%fitted on 20%
f = fit(cat(1,data.mut21,data.mut12),cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2),ft,fo)

%correlation computation on 80%
fitted_function_single_mut_background = ft(f.b,f.c,cat(1,data1.mut21,data1.mut12));
expt_double_fitness = data1.log_fitness;
expt_mut_1 = data1.Mut1_expt_fitness;
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;
expt_mut_2 = data1.Mut2_expt_fitness;
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;
R = corrcoef(fitted_function_single_mut_background, cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2));
R(1,2)