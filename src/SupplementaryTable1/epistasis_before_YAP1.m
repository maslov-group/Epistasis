%read YAP-1 dataset (80%)
data=readtable("80_llm_yap1_esm2_650.csv");

data.expt_epistasis = data.log_fitness - (data.Mut1_expt_fitness + data.Mut2_expt_fitness);

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%LLM predicted epistasis
total_epistasis = 0.5*(data.mut1+data.mut21+data.mut2+data.mut12) - (data.mut1+data.mut2);

%experimental epistasis
data.total_epistasis = total_epistasis;

[R,p] = corrcoef(data.total_epistasis, data.expt_epistasis);
R(1,2)
p(1,2)