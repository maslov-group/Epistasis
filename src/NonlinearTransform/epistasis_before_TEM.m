%read TEM-1 dataset
data=readtable("llm_tem1_esm2_650.csv");

data.log_fitness = log(data.DoubleMutantFitness);
data.Mut1_expt_fitness = log(data.Mut1Fitness);
data.Mut2_expt_fitness = log(data.Mut2Fitness);

data.expt_epistasis = data.log_fitness - (data.Mut1_expt_fitness + data.Mut2_expt_fitness);

ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b','c'});
fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

%calculate LLM predicted epistasis
total_epistasis = 0.5*(data.mut1+data.mut21+data.mut2+data.mut12) - (data.mut1+data.mut2);

%calculate experimental epistasis
expt_epistasis = data.log_fitness - (data.Mut1_expt_fitness + data.Mut2_expt_fitness);

[R,p] = corrcoef(total_epistasis, expt_epistasis);
R(1,2)
p(1,2)

[rho,pval] = corr(total_epistasis, expt_epistasis,'Type','Spearman');
rho
pval