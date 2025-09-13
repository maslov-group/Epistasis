%read TEM-1 dataset
data=readtable("80_llm_tem1_esm2_650.csv");

%independent x = mut21, dependent y = log(expt_double/expt_mut1)
expt_double_fitness = log(data.DoubleMutantFitness);
expt_mut_1 = log(data.Mut1Fitness);
log_exptdouble_exptmut1 = expt_double_fitness - expt_mut_1;

%independent x = mut12, dependent y = log(expt_double/expt_mut2)
expt_mut_2 = log(data.Mut2Fitness);
log_exptdouble_exptmut2 = expt_double_fitness - expt_mut_2;

R = corrcoef(cat(1,data.mut21,data.mut12), cat(1,log_exptdouble_exptmut1,log_exptdouble_exptmut2));

R(1,2)