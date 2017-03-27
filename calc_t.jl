using DataFrames
using HypothesisTests
using Distributions
using KernelDensity
using Optim
using StatsBase
@everywhere include("sne_init_4.jl");
n_exp = 2;


@everywhere bootstrap_T(1)

tic();
my_ts_pmap = pmap(bootstrap_T,n_exp.*ones(Int,nworkers()));
toc();
my_ts = my_ts_pmap[1];

for i in 2:nworkers()
        my_ts = vcat(my_ts,my_ts_pmap[i]);
end

writecsv(string("../results/T_",n_exp*nworkers()),my_ts);
