%fit using 20%, correlation on 80%
%fix a=1, keep same b and c from single mutation fitting
%epistasis = 0.5*(mut1_prime+mut21_prime+mut2_prime+mut12_prime)-(mut1_prime+mut2_prime);

data=readtable("RRM2_sample_20_650M_esm2.csv");
data1=readtable("RRM2_sample_20_650M_esm2_unique_muts.csv");
data_test = readtable("RRM2_sample_80_650M_esm2.csv");

data.expt_epistasis = data.expMutDouble - (data.expMut1 + data.expMut2);
%data = data1(data1.expt_epistasis>0,:);  %positive values only

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b', 'c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

f = fit(data1.llm_mut,data1.exp_mut,ft,fo);
%f = fit(cat(1,data.mut1,data.mut2), cat(1,data.Mut1_expt_fitness,data.Mut2_expt_fitness),ft,fo);
f.b
f.c

mut1_prime = ft(f.b,f.c,data_test.mut1); %predicted mut1_prime based on the transform
mut2_prime = ft(f.b,f.c,data_test.mut2); %predicted mut2_prime based on the transform

log_exp_double_exp_mut2_ = data.expMutDouble - data.expMut2;
log_exp_double_exp_mut1_ = data.expMutDouble - data.expMut1;

f = fit(cat(1,data.mut21,data.mut12), cat(1,log_exp_double_exp_mut1_,log_exp_double_exp_mut2_),ft,fo);
f.b
f.c

mut12_prime = ft(f.b,f.c,data_test.mut12);     %predicted mut12_prime based on the transform
mut21_prime = ft(f.b,f.c,data_test.mut21);     %predicted mut21_prime based on the transform

total_epistasis = 0.5*(mut1_prime+mut21_prime+mut2_prime+mut12_prime)-(mut1_prime+mut2_prime);

data_test.expt_epistasis = data_test.expMutDouble - (data_test.expMut1 + data_test.expMut2);
data_test.llm_epistasis = total_epistasis;
writetable(data_test, 'epi_RRM2_esm2_650M_20_80.csv');


[R,p] = corrcoef(total_epistasis, data_test.expt_epistasis);
R(1,2)
p(1,2)