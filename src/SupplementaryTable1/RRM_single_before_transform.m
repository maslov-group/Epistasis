data=readtable("RRM2_sample_20_650M_esm2_unique_muts.csv");
ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%independent x = llm_single_mut, dependent y = expt_single_mut
R = corrcoef(data.llm_mut, data.exp_mut);

R(1,2)

rho = corr(data.llm_mut, data.exp_mut,'Type','Spearman');

rho