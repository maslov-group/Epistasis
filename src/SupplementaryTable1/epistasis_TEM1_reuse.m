%read dataset and unique single mutations
data=readtable("80_llm_tem1_esm2_650.csv");
data1=readtable("20_esm2_650_unique_single_mutations_TEM1.csv");

data.log_fitness = log(data.DoubleMutantFitness);
data.Mut1_expt_fitness = log(data.Mut1Fitness);
data.Mut2_expt_fitness = log(data.Mut2Fitness);
data.expt_epistasis = data.log_fitness - (data.Mut1_expt_fitness + data.Mut2_expt_fitness);

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);
%independent x = llm_single_mut, dependent y = expt_single_mut
f = fit(data1.llm_single_mut,data1.expt_single_mut,ft,fo);

f.b
f.c 

mut1_prime = ft(f.b,f.c,data.mut1); %predicted mut1_prime based on the transform
mut2_prime = ft(f.b,f.c,data.mut2); %predicted mut2_prime based on the transform

log_exp_double_exp_mut2_ = data.log_fitness - data.Mut2_expt_fitness;
log_exp_double_exp_mut1_ = data.log_fitness - data.Mut1_expt_fitness;

mut12_prime = ft(f.b,f.c,data.mut12);     %predicted mut12_prime based on the transform
mut21_prime = ft(f.b,f.c,data.mut21);     %predicted mut21_prime based on the transform

%calculate LLM predicted epistasis
total_epistasis = 0.5*(mut1_prime+mut21_prime+mut2_prime+mut12_prime)-(mut1_prime+mut2_prime);

[R,p] = corrcoef(total_epistasis, data.expt_epistasis);
R(1,2)
p(1,2)
