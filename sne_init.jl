using HypothesisTests
using Distributions
using KernelDensity
using Optim
using Interpolations
using StatsBase

include("cosmology.jl");

immutable nu
    mjd::Float64
    log_eng::Float64
    ang_err::Float64
    ra::Float64
    dec::Float64
end

type sn
  max_date::Float64
  ra::Float64
  dec::Float64
  z::Float64
  zenith_bin::Int
  nb::Float64
  associated_nus::Array{nu,1}
  coefs::Array{Float64,1}

  function sn(max_date::Real,ra::Float64,dec::Float64,z::Float64)
      new_sn = new(max_date,ra,dec,z);
      zenith_bin = -1;
      new_sn.associated_nus=nu[];
      new_sn.coefs =Float64[];

      return new_sn
  end
end

add_coefs!(t_sn::sn, coefs::Array{Float64,1}) = t_sn.coefs=coefs;
add_nb_and_zenith_bin_idx!(t_sn::sn, nb::Float64, kk::Int) = (t_sn.nb=nb; t_sn.zenith_bin = kk;)
add_nu!(t_sn::sn, nu::Array{nu,1}) = append!(t_sn.associated_nus,nu);
add_nu!(t_sn::sn, nu::nu) = append!(t_sn.associated_nus,nu);

rm_associated_nus!(t_sn) = t_sn.associated_nus=nu[];
rm_coefs!(t_sn) = t_sn.coefs=Float64[];

function get_nu_sample()
    nu_sample= Array(nu,len_nu);

    t_nu_sample = calc_sample_nus();

    nu_sample[:] = [nu(t_nu_sample[i,:]...) for i in 1:len_nu];


    return nu_sample
end

function calc_sample_nus(N_nus=len_nu)
    mjd = sample(nu_data[:,1],replace=true,N_nus);
    eng = sample(nu_data[:,2],replace=true,N_nus);
    ang_err = sample(nu_data[:,3],replace=true,N_nus);
    ra = rand(Uniform(0,2pi),N_nus);
    dec = sample(nu_data[:,5],replace=true,N_nus);

    return hcat(mjd,eng,ang_err,ra,dec)
end

function find_rand_N_nus(t_sn::sn)
    if t_sn.z > 0.0
        t_z = t_sn.z;
    else
        rand_z = -1.0;
        while rand_z <= 0.0
            sn_idx = round(Int,(len_sne-1)*rand()+1.0);
            rand_z = my_sn[sn_idx].z;
        end
        t_z = rand_z;
    end
    @assert(t_z > 0.0)
    t_D_L = calc_D_L(t_z);
    num_flux = erg_2_GeV.*nu_flux_coef/t_D_L^2/(1e5)^2;
    return rand(Poisson(unnormed_num_nus[t_sn.zenith_bin]*num_flux));

end

function calc_sig_sample_nus(t_sn::sn,t_N_nus::Int)
    if t_N_nus < 1
        return Array(Float64,0,5)
    else
        ang_err = sample(nu_data[:,3],replace=true,t_N_nus);
        mjd = (-0.5-rand(Poisson(13.0))+t_sn.max_date).*ones(Float64,t_N_nus);

        eng = eng_cdf_fn[t_sn.zenith_bin][rand(t_N_nus)];

        kappa = 1./ang_err.^2;

        @assert(minimum(kappa) .> 10.0)
        sample_mu = (log(rand(t_N_nus))+kappa)./kappa;
        #println(kappa)
        theta_prime_nu = acos(sample_mu);
        phi_prime_nu = 2pi*rand(t_N_nus);

        theta_sn = 0.5*pi-t_sn.dec;
        phi_sn = t_sn.ra;


        tan_phi_nu = sin(phi_prime_nu-phi_sn).*sin(theta_prime_nu);
        tan_phi_nu ./= cos(phi_prime_nu).*cos(phi_sn).*cos(theta_sn).*sin(theta_prime_nu) +
            cos(theta_sn).*sin(phi_prime_nu).*sin(phi_sn).*sin(theta_prime_nu)-
            cos(theta_prime_nu).*sin(theta_sn);

        sin_dec_nu = cos(theta_prime_nu).*cos(theta_sn)+cos(phi_prime_nu-phi_sn).*sin(theta_prime_nu).*sin(theta_sn);

        ra = atan(tan_phi_nu);
        dec = asin(sin_dec_nu);

        return hcat(mjd,eng,ang_err,ra,dec)
    end
end

function get_sample_sig_nu!(t_sn::Array{sn,1},t_nu_sample::Array{Float64,2})
    t_len_sne = length(t_sn);
    num_nus = Array(Int,t_len_sne);

    num_nus[:] = map(find_rand_N_nus,t_sn);

    nu_sample_idx = 1;
    for i in 1:t_len_sne
        if num_nus[i] > 0
            t_nu_sample[nu_sample_idx:nu_sample_idx+num_nus[i]-1,:] = calc_sig_sample_nus(t_sn[i],num_nus[i]);
            nu_sample_idx += num_nus[i];
        end
    end
    return nu_sample_idx;
end

function find_associated_nus(t_sn::sn, t_nus::Array{nu,1})
    n_hi = 19;
    n_lo = 4;

    in_time_window = [t_sn.max_date - n_hi <= t_nus[j].mjd <= t_sn.max_date - n_lo for j in 1:len_nu]

    in_ang_window = [(sin(t_sn.dec)*sin(t_nus[j].dec) + cos(t_sn.dec)*cos(t_nus[j].dec)*
        cos(t_sn.ra-t_nus[j].ra)) > cos(5.0*t_nus[j].ang_err) for j in 1:len_nu]

    return in_time_window.*in_ang_window;
end



## Function that reads in the pre-calculated values of the energy PDFs for
## Signal and background neutrinos. Signals neutrinos are assumed to come from an E^-2 flux.
## The PDFs vary as a function of neutrino energy proxy, and zenith angle.

function get_energy_pdf!(t_zenith::Array{Float64,1}, t_proxy::Array{Float64,1})

    #AA = readcsv("../data/sig_num_nus_pdf");
    #BB = readcsv("../data/atm_num_nus_pdf");

    AA = readcsv("w_energy_dep/data/sig_num_nus_pdf");
    BB = readcsv("w_energy_dep/data/atm_num_nus_pdf");



    t_zenith[:] = AA[:,1];
    t_proxy[:] = AA[1,2:end];

    ## Array of energy PDF functions. Each element corresponds to a zenith bin.
    ## This function is the result of the interpolation of the pre-calculated PDF data points

    wrapper_log_sig_energy_pdf = Array(Function,length(t_zenith)-1);
    wrapper_log_atm_energy_pdf = Array(Function,length(t_zenith)-1);

    ## read in the pre-calculated values of the energy PDFs

    t_sig_energy_pdf = AA[2:end,2:end];
    t_atm_energy_pdf = BB[2:end,2:end];

    ## create a wrapper function with the max and min values of the proxy energy already set

    t_make_interp_fn = x -> make_interp_fn(x,log10(proxy[1]),log10(proxy[end]),log10(proxy[2]/proxy[1]));

    for j in 2:len_zenith
        wrapper_log_sig_energy_pdf[j-1] = t_make_interp_fn(log10(AA[j,2:end]));
        wrapper_log_atm_energy_pdf[j-1] = t_make_interp_fn(log10(BB[j,2:end]));
    end

    return wrapper_log_sig_energy_pdf, wrapper_log_atm_energy_pdf;
end

## function that allows a wraper function to be created for a constant x_min, x_max, and Dx
## which only depends on an array of data points fed in. Creates an interpolation object that
## becomes a function that can be referenced at any point x_min < x < x_max

function make_interp_fn(zz::Array{Float64,1}, uu_min::Float64, uu_max::Float64, Du::Float64)

    t_interp = interpolate(zz,BSpline(Quadratic(Natural())),OnGrid());

    return x -> eval_interp(x,t_interp,uu_min,uu_max,Du)
end


## Function that allows for easy evaluation of an interpolated object.
## Once given information about the interpoluation, can be called simply by wrapping the function


function eval_interp(xx::Float64, yy::AbstractInterpolation, xx_min::Float64, xx_max::Float64, Dx::Float64)

    ii = (xx - xx_min)/Dx + 1.0;
    return 10^yy[ii];
end

## Create a neutrino declination PDF based the true background neutrinos.
## their indecies have already been calculated and must be input.
## the PDF is a function only the of the neutrino Dec.

function get_dec_pdf(t_nu_dec::Array{Float64,1})

    nu_dec_kde = kde(t_nu_dec);
    nu_dec_interp = InterpKDE(nu_dec_kde);

    return xx ->  pdf(nu_dec_interp,xx)./cos(xx);
end

function S_dir(t_sn::sn,t_nu::nu)
    mu = sin(t_sn.dec)*sin(t_nu.dec) + cos(t_sn.dec)*cos(t_nu.dec)*cos(t_sn.ra-t_nu.ra);
    kappa = 1/t_nu.ang_err^2;

    if kappa > 10.0
        result = kappa/(2*pi)*exp(kappa*(mu-1.0))
    elseif 0.0 < kappa <= 10.0
        result = kappa/(4*pi*sinh(kappa))*exp(kappa*mu);
    else
        error("strange kappa");
    end

    @assert(!isnan(result))

    return result;
    #return exp(-0.5*Delta_Psi^2/t_nu.ang_err^2)/(2*pi*t_nu.ang_err^2);
end


## is it ok to use a Probability Mass Function for the time signal liklihood if my
## Background is expressed as a Probability Density Function? I think it is if you take the
## view that the PMF is a PDF with step values at each integer divided by the spacing between each
## each integer (i.e., 1). In theory you could do this integration and it should come to 1.

function S_time(t_sn::sn, t_nu::nu)

    return pdf(Poisson(13.0),t_sn.max_date-t_nu.mjd-0.5)
end

function calc_T(nns::Array{Float64,1},t_sn::Array{sn,1})
     inner_sum = [sum(log(nns[i].*t_sn[i].coefs+ 1.0)) for i in 1:len_sne]

    return sum( nns .- inner_sum);
end

function grad_T!(x::Vector, storage::Vector, t_sn::Array{sn,1})
    storage[:] = -[sum(t_sn[i].coefs./(t_sn[i].coefs.*x[i]+1.0)) - 1.0 for i in 1:length(x)];
end

function hess_T!(x::Vector, storage::Matrix, t_sn::Array{sn,1})
    fill!(storage,0.0);
    for i in 1:length(x)
        storage[i,i] = sum((t_sn[i].coefs./(t_sn[i].coefs.*x[i]+1.0)).^2)
    end
end

function calc_coefs!(t_sn::sn)

    if length(t_sn.associated_nus) > 0
	    temp_s = map(x-> S_dir(t_sn,x),t_sn.associated_nus);
    	temp_s .*= map(x-> S_time(t_sn,x),t_sn.associated_nus);
    	temp_s .*= map(x-> 10^wrapper_log_sig_energy_pdf[t_sn.zenith_bin](x),[t_sn.associated_nus[j].log_eng for j in 1:length(t_sn.associated_nus)]);

    	temp_b = (1/(2*pi))*B_nu_dec([t_sn.associated_nus[j].dec for j in 1:length(t_sn.associated_nus)]);
    	temp_b ./= 16.0;
    	temp_b .*= map(x-> 10^wrapper_log_atm_energy_pdf[t_sn.zenith_bin](x),[t_sn.associated_nus[j].log_eng for j in 1:length(t_sn.associated_nus)]);

    	add_coefs!(t_sn,temp_s./temp_b./t_sn.nb);
    end
end


## Calculate the background rate in each zenith bin for the time window of interest
## in this case since n_hi = 19 and n_lo = 4 the time window is 16 days (although this needs to be checked)

function get_nb(t_nu_dec::Array{Float64},zz; minimum_nu_date = 55694.0, maximum_nu_date = 56062.0)
  cos_zenith_binned_bkg = fit(Histogram,cos(0.5*pi + t_nu_dec),zz);
  println(cos_zenith_binned_bkg.weights)
  return 16.0.*cos_zenith_binned_bkg.weights./(maximum_nu_date - minimum_nu_date);
end


function assign_nb!(t_sn::Array{sn,1}, nnb::Array{Float64,1}, ks::Array{Int,1})
    map((x,y,z)->add_nb_and_zenith_bin_idx!(x,y,z),t_sn,nnb,ks);
    nothing;
end

function get_T()
    srand(13)
    nu_sample = Array(Float64,len_nu,5);
    #nu_sample_idx = get_sample_sig_nu!(my_sn,nu_sample);
    #if no sample nus set nu_sample_idx to 1
    nu_sample_idx = 1;
    nu_sample[nu_sample_idx:end,:] = calc_sample_nus(len_nu-nu_sample_idx+1)
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

    opt_T = optimize(OnceDifferentiable(ns-> calc_T(ns,my_sn), (x,s) -> grad_T!(x,s,my_sn)),zeros(Float64,len_sne),lower,upper,Fminbox(),optimizer = ConjugateGradient);
    #opt_T = optimize(ns-> calc_T(ns,my_sn),(x,s) -> grad_T!(x,s,my_sn),(xx,ss) -> hess_T!(xx,ss,my_sn),zeros(Float64,len_sne),Newton())

    return opt_T.minimizer, opt_T.minimum
end

function bootstrap_T(NN::Int)
        my_ts = Array(Float64,NN,len_sne+1);
        for N in 1:NN
                t_opt_T = get_T();
                my_ts[N,1:end-1] = t_opt_T[1];
		            my_ts[N,end] = t_opt_T[2];
        end
        return my_ts;
end

#nu_data = readdlm("../data/modified_upgoing_nu");
#sne_data = readdlm("../data/sne_data");
srand(13);
nu_data = readdlm("w_energy_dep/data/modified_upgoing_nu");
sne_data = readdlm("w_energy_dep/data/sne_data");

len_sne = 27;
len_nu = 69227;




my_sn = Array(sn,len_sne);
my_nu = Array(nu,len_nu);

my_sn[:] = [sn(sne_data[i,:]...) for i in 1:len_sne];

len_zenith = 12;
len_proxy = 51;

zenith = Array(Float64,len_zenith);
proxy = Array(Float64,len_proxy);

#bkg_nus = convert(Array{Int,1},readcsv("../data/bkg_nus_indx")[:]);
bkg_nus = convert(Array{Int,1},readcsv("w_energy_dep/data/bkg_nus_indx")[:]);

B_nu_dec = get_dec_pdf(nu_data[bkg_nus,5]);
wrapper_log_sig_energy_pdf, wrapper_log_atm_energy_pdf = get_energy_pdf!(zenith,proxy);

@assert(len_zenith == length(zenith));
@assert(len_proxy == length(proxy));

sig_eng_cdf_fn = Array(AbstractInterpolation,len_zenith);
#sig_eng_cdf = readdlm("../data/sig_eng_cdf");
sig_eng_cdf = readdlm("w_energy_dep/data/sig_eng_cdf");

log_eng = sig_eng_cdf[:,1];

sig_eng_cdf_fn[:] = [interpolate((sort(sig_eng_cdf[:,i]),),log_eng, Gridded(Linear())) for i in 1:len_zenith] ;


#assign_nb!(my_sn,readdlm("../data/nb_data")[1:len_sne,1],convert(Array{Int,1},readdlm("../data/ks")[1:len_sne,1]))
assign_nb!(my_sn,readdlm("w_energy_dep/data/nb_data")[1:len_sne,1],convert(Array{Int,1},readdlm("w_energy_dep/data/ks")[1:len_sne,1]))


sig_cdf = readdlm("w_energy_dep/data/sig_eng_cdf");
unnormed_num_nus = readdlm("w_energy_dep/data/unnormalized_number_of_neutrinos");
cdf_log_E = sig_cdf[:,1];
eng_cdf_fn = Array(AbstractInterpolation, len_zenith-1);
eng_cdf_fn[:] = [interpolate((sort(sig_cdf[:,j]),),cdf_log_E, Gridded(Linear())) for j in 2:len_zenith];

E_cr = 1e52;
C = 18;

nu_flux_coef = (3/8)*E_cr./(4*pi*C);
