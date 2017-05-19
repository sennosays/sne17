using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase

@everywhere include("sne_init_stationary.jl");
@everywhere get_T(1e52,1.0);
n_exp = 1000;
N_E_cr = 10;
N_sn_frac = 10;

my_E_cr_ex = linspace(46,52,N_E_cr);
my_E_cr = 10.^my_E_cr_ex;

my_sn_frac = linspace(0.0,1.0,N_sn_frac);


my_ts = SharedArray(Float64,n_exp,len_sne+1);
@time for ii in 1:N_E_cr
	for jj in 1:N_sn_frac
		@sync @parallel for n in 1:n_exp
		        my_ts[n,1:end-1], my_ts[n,end] = get_T(my_E_cr[ii],my_sn_frac[jj]);
		end
		@sync writedlm(string("results/T_",n_exp,"_",round(log10(my_E_cr[ii]),1),"_",round(my_sn_frac[jj],2),"_stationary.dat"),my_ts);
	end
end
