%read dataset and unique single mutations
data=readtable("RRM2_sample_80_650M_esm2.csv");
data1=readtable("RRM2_sample_20_650M_esm2_unique_muts.csv");

data.expt_epistasis = data.expMutDouble - (data.expMut1 + data.expMut2);

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);
%independent x = llm_mut, dependent y = exp_mut
f = fit(data1.llm_mut,data1.exp_mut,ft,fo);

f.b
f.c

mut1_prime = ft(f.b,f.c,data.mut1); %predicted mut1_prime based on the transform
mut2_prime = ft(f.b,f.c,data.mut2); %predicted mut2_prime based on the transform

log_exp_double_exp_mut2_ = data.expMutDouble - data.expMut2;
log_exp_double_exp_mut1_ = data.expMutDouble - data.expMut1;

mut12_prime = ft(f.b,f.c,data.mut12);     %predicted mut12_prime based on the transform
mut21_prime = ft(f.b,f.c,data.mut21);     %predicted mut21_prime based on the transform

%calculate LLM predicted epistasis
total_epistasis = 0.5*(mut1_prime+mut21_prime+mut2_prime+mut12_prime)-(mut1_prime+mut2_prime);

[R,p] = corrcoef(total_epistasis, data.expt_epistasis);
R(1,2)
p(1,2)