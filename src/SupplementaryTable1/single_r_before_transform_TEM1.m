%read dataset
data=readtable("80_esm2_650_unique_single_mutations_TEM1.csv");

%independent x = llm_single_mut, dependent y = expt_single_mut
R = corrcoef(data.llm_single_mut, data.expt_single_mut);

R(1,2)