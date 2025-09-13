%read TEM-1 dataset
data=readtable("20_llm_tem1_esm2_650.csv");
data1=readtable("80_llm_tem1_esm2_650.csv");
ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = mut21, dependent y = log(expt_double/expt_mut1)
expt_double_fitness = log(data.DoubleMutantFitness);
expt_mut_1 = log(data.Mut1Fitness);
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;

%independent x = mut12, dependent y = log(expt_double/expt_mut2)
expt_mut_2 = log(data.Mut2Fitness);
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;

%fitted on 20%
f = fit(cat(1,data.mut21,data.mut12),cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2),ft,fo)

%correlation computation on 80%
fitted_function_single_mut_background = ft(f.b,f.c,cat(1,data1.mut21,data1.mut12));
expt_double_fitness = log(data1.DoubleMutantFitness);
expt_mut_1 = log(data1.Mut1Fitness);
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;
expt_mut_2 = log(data1.Mut2Fitness);
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;
R = corrcoef(fitted_function_single_mut_background, cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2));
R(1,2)