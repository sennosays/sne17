include("sne_init.jl");

function precondprepbox_sn!(P, x, l, u, mu, t_sn::Array{sn,1})
        @inbounds @simd for i = 1:length(x)
        xi = x[i]
        li = l[i]
        ui = u[i]
        cs = t_sn[i].coefs;
        #P.diag[i] = 1/(mu*(1/(xi-li)^2 + 1/(ui-xi)^2) + 1) # +1 like identity far from edges
        P.diag[i] = mu*sum((cs./(cs.*xi+1.0)).^2)
    end
end

srand(13)
nu_sample = calc_sample_nus();
my_nu[:] = [nu(nu_sample[j,:]...) for j in 1:len_nu];

map(rm_coefs!,my_sn)
map(rm_associated_nus!,my_sn)

for i in 1:len_sne
  associated_nus = find_associated_nus(my_sn[i],my_nu);
  if length(associated_nus) > 0
    add_nu!(my_sn[i],my_nu[associated_nus])
  end
end

#map(x-> add_nu!(x,my_nu[find_associated_nus(x,my_nu)]),my_sn);

map(calc_coefs!,my_sn);

lower = [maximum(-1./my_sn[i].coefs) for i in 1:len_sne];
upper = Inf.*ones(Float64,len_sne);

println("GradientDescent w/ prep")
opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = GradientDescent);

@time opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = GradientDescent);

println(opt_T_prep.minimum,"\n")

println("GradientDescent no prep")
opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = GradientDescent);

@time opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = GradientDescent);

println(opt_T_no_prep.minimum,"\n")

println("ConjugateGradient w/ prep")
opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = ConjugateGradient);

@time opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = ConjugateGradient);

println(opt_T_prep.minimum,"\n")

println("ConjugateGradient no prep")
opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = ConjugateGradient);

@time opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = ConjugateGradient);

println(opt_T_no_prep.minimum,"\n")

println("LBFGS w/ prep")

opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = LBFGS);

@time opt_T_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),precondprep = (P, x, l, u, mu) -> precondprepbox_sn!(P, x, l, u, mu, my_sn),
    optimizer = LBFGS);

println(opt_T_prep.minimum,'\n');

println("LBFGS no prep")
opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = LBFGS);

@time opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
    zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = LBFGS);

println(opt_T_no_prep.minimum,'\n');

println("NelderMead")
opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = NelderMead);

@time opt_T_no_prep = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),
zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = NelderMead);

println(opt_T_no_prep.minimum);
