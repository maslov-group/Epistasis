%fix a=1, keep same b and c from single mutation fitting
%fit on 20% data, correlation on 80% data

approved_proteins = readtable("approved_proteins.csv");

files = dir(fullfile("outputs", "20_esm2_650M_*.csv"));

n = height(approved_proteins);

file_names = strings(n,1);
b1_vals = zeros(n,1);
c1_vals = zeros(n,1);
b2_vals = zeros(n,1);
c2_vals = zeros(n,1);
R_vals  = zeros(n,1);
p_vals  = zeros(n,1);
R_before_vals  = zeros(n,1);
p_before_vals  = zeros(n,1);

for i = 1:n
    base = string(approved_proteins.file_name(i));
    try
        
        file1 = fullfile("outputs", "20_esm2_650M_" + base);
        file2 = fullfile("outputs", "20_esm2_650M_unique_single_mutations_" + base);
        file3 = fullfile("outputs", "80_esm2_650M_" + base);
        
        data = readtable(file1);
        data1 = readtable(file2);
        data_test = readtable(file3);
    
        %calculate epistasis before non-linear fit        
        ft = fittype('-1.*log(1+exp(-b.*(x+c)))','dependent',{'y'},'independent',{'x'},'coefficients',{'b', 'c'});
        fo = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0]);

        total_epistasis_before = 0.5*(data_test.mut1+data_test.mut21+data_test.mut2+data_test.mut12) - (data_test.mut1+data_test.mut2);
        data_test.total_epistasis_before = total_epistasis_before;
        data_test.expt_epistasis = data_test.DoubleMutantFitness - (data_test.Mut1Fitness + data_test.Mut2Fitness);
        
        f = fit(data_test.total_epistasis_before, data_test.expt_epistasis, ft,fo);

        [R_before,p_before] = corrcoef(data_test.total_epistasis_before, data_test.expt_epistasis);

        f = fit(data1.llm_single_mut, data1.expt_single_mut, ft, fo);
        b1 = f.b;
        c1 = f.c;
    
        mut1_prime = ft(f.b,f.c,data_test.mut1); %predicted mut1_prime based on the transform
        mut2_prime = ft(f.b,f.c,data_test.mut2); %predicted mut2_prime based on the transform
    
        log_exp_double_exp_mut2_ = data.DoubleMutantFitness - data.Mut2Fitness;
        log_exp_double_exp_mut1_ = data.DoubleMutantFitness - data.Mut1Fitness;
    
        f = fit(cat(1,data.mut21,data.mut12), cat(1,log_exp_double_exp_mut1_,log_exp_double_exp_mut2_),ft,fo);
        b2 = f.b;
        c2 = f.c;
        
        mut12_prime = ft(f.b,f.c,data_test.mut12);     %predicted mut12_prime based on the transform
        mut21_prime = ft(f.b,f.c,data_test.mut21);     %predicted mut21_prime based on the transform
    
        %calculate LLM predicted epistasis after non-linear fit
        total_epistasis = 0.5*(mut1_prime+mut21_prime+mut2_prime+mut12_prime)-(mut1_prime+mut2_prime);
    
        data_test.expt_epistasis = data_test.DoubleMutantFitness - (data_test.Mut1Fitness + data_test.Mut2Fitness);
        [R,p] = corrcoef(total_epistasis, data_test.expt_epistasis);
        
        file_names(i) = base;
        b1_vals(i) = b1;
        c1_vals(i) = c1;
        b2_vals(i) = b2;
        c2_vals(i) = c2;
        R_vals(i)  = R(1,2);
        p_vals(i)  = p(1,2);
        R_before_vals(i)  = R_before(1,2);
        p_before_vals(i)  = p_before(1,2);
    catch ME
        base
        disp(ME.message)
        continue;
    end
    
end

results_table = table(file_names, b1_vals, c1_vals, b2_vals, c2_vals, R_before_vals, p_before_vals, R_vals, p_vals);

results_table.Properties.VariableNames = ["file","b1","c1","b2","c2","R_before","p_before","R","p"];

writetable(results_table, fullfile("outputs","epistasis_results_25.csv"));