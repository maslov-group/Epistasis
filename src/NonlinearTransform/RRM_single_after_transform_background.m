data=readtable("RRM2_sample_20_650M_esm2.csv");
ft = fittype('-0.5.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = mut21, dependent y = log(expt_double/expt_mut1)
%expt_double_fitness = log(data.DoubleMutantFitness);
%expt_mut_1 = log(data.Mut1Fitness);
expt_double_fitness = data.expMutDouble;
expt_mut_1 = data.expMut1;
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;

%independent x = mut12, dependent y = log(expt_double/expt_mut2)
%expt_mut_2 = log(data.Mut2Fitness);
expt_mut_2 = data.expMut2;
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;

f = fit(cat(1,data.mut21,data.mut12),cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2),ft,fo)
fitted_function_single_mut_background = ft(f.b,f.c,cat(1,data.mut21,data.mut12));

R = corrcoef(fitted_function_single_mut_background, cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2));
R(1,2)