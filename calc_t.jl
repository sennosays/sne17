using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase
@everywhere include("sne_init.jl");
@everywhere get_T();
n_exp =10

my_ts = SharedArray(Float64,n_exp,len_sne+1);
tic()
@sync @parallel for n in 1:n_exp
       # println("hello")
        t_opt_T = get_T();
        my_ts[n,1:end-1] = t_opt_T[1];
        my_ts[n,end] = t_opt_T[2];
       # println("goodbye");
end
toc();
#rmprocs(workers())
writedlm(string("../results/T_",n_exp,".dat"),my_ts);
