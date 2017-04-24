function dt_star_dz(this_z)

  this_opz = 1.0 + this_z;

  return 1.0/(H_0_s*this_opz*sqrt(Omega_m*this_opz^3 + Omega_L));

end

calc_D_L(this_z) = (1.0+this_z)*(calc_D_prop(this_z));

calc_D_prop(this_z) = c*quadgk(a-> (1.0+a)*dt_star_dz(a),0.0,this_z)[1];
