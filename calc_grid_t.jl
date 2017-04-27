using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase

@everywhere include("sne_init.jl");
@everywhere get_T(1e52,1.0);
n_exp = 10000;
N_E_cr = 10; 
N_sn_frac = 10; 

my_E_cr_ex = linspace(46,52,N_E_cr); 
my_E_cr = 10.^my_E_cr_ex; 

my_sn_frac = linspace(0.1,1.0,N_sn_frac); 


my_ts = SharedArray(Float64,n_exp,len_sne+1);
@time for ii in 1:N_E_cr 
	for jj in 1:N_sn_frac 
		@sync @parallel for n in 1:n_exp
			println("E_cr: ",my_E_cr[ii]," f_sn: ",my_sn_frac[jj]); 
        		t_opt_T = get_T(my_E_cr[ii],my_sn_frac[jj]);
        		my_ts[n,1:end-1] = t_opt_T[1];
        		my_ts[n,end] = t_opt_T[2];
		end
		@sync writedlm(string("results/T_",n_exp,"_",my_E_cr[ii],"_",my_sn_frac[jj],".dat"),my_ts);
	end
end

