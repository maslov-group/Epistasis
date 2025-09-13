%read dataset
data=readtable("20_esm2_650_unique_single_mutations_TEM1.csv");
data1 = readtable("80_esm2_650_unique_single_mutations_TEM1.csv");
ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = llm_single_mut, dependent y = expt_single_mut
%fit on 20%
f = fit(data.llm_single_mut,data.expt_single_mut,ft,fo)

%correlation on 80%
fitted_function_single_mut = ft(f.b,f.c,data1.llm_single_mut);
R = corrcoef(fitted_function_single_mut, data1.expt_single_mut);

R(1,2)