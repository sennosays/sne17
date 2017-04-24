using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase

@everywhere include("sne_init.jl");
@everywhere get_T(1e52,1.0);
n_exp = 10000;

my_E_cr = 1e46; 
my_sn_frac = 1.0; 


my_ts = SharedArray(Float64,n_exp,len_sne+1);
#tic()
@time @sync @parallel for n in 1:n_exp
        t_opt_T = get_T(my_E_cr,my_sn_frac);
        my_ts[n,1:end-1] = t_opt_T[1];
        my_ts[n,end] = t_opt_T[2];
end
#toc();
writedlm(string("results/T_",n_exp,"_",my_E_cr,"_",my_sn_frac,".dat"),my_ts);
