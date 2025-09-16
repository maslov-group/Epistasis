%fit using 20%, correlation on 80%
%read dataset and unique single mutations
data=readtable("RRM2_sample_20_650M_esm2_unique_muts.csv");
data1=readtable("RRM2_sample_80_650M_esm2_unique_muts.csv");

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = llm_mut, dependent y = exp_mut
f = fit(data.llm_mut,data.exp_mut,ft,fo)
fitted_function_single_mut = ft(f.b,f.c,data1.llm_mut);
R = corrcoef(fitted_function_single_mut, data1.exp_mut);

R(1,2)