using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase

@everywhere include("sne_init.jl");
@everywhere get_T(1e52,1.0);
n_exp = 1000;



my_ts = SharedArray(Float64,n_exp,len_sne+1);
#tic()
@time @sync @parallel for n in 1:n_exp
	my_ts[n,1:end-1], my_ts[n,end] = get_T(1e46,0.0);
end
#toc();
writedlm(string("results/T_",n_exp,"_46.0_0.0_new.dat"),my_ts);
