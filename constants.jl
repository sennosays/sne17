#Astrophysical Constants
const c = 29979245800.0;               #speed of light [cm/s]
const M_sun = 1.9891e33;             #Solar mass [g]
const k_T_0 = 8.6173325e-5 * 2.72548 #CMB blackbody temperature at z=0 [eV]

#Cosmological Constants
const H_0 = 67.74                    #Hubble constant (Plank 2015) [km s^-1 Mpc^-1]
const H_0_s = 2.18e-18               #Hubble constant (Plank 2015) [s^-1]
const Omega_m = 0.3089               #Mass fraction (Plank 2015) [1]
const Omega_L = 0.6911               #Vacuum energy fraction (Plank 2015) [1]

#Particle Constants
const charge = 4.80321e-10;    #elementary charge [esu]

const h_cgs = 6.626075540e-27  #Planck's constant in cgs units [erg s]
const hc = 1.24e-4;            #Plank's constant times speed of light [eV cm]
const hbar = 6.5821192815e-16  #Hbar [eV s]

const m_p = 1.6726e-24;        #mass of proton [g]
const rm_p = 0.938257          #rest mass energy of proton [GeV]

const m_e = 9.10938291e-28     #mass of electorn [g]
const rm_e = 5.11e-4           #rest mass energy of electron [GeV]

const rm_mu = 0.105658371535   #rest mass energy of muon [GeV]
const tau_mu = 2.2e-6          #lifetime of the muon [s]

const rm_pi_0 = 0.13497666     #rest mass energy of neutral pion [GeV]
const rm_pi_c = 0.1395701835   #rest mass energy of charged pion [GeV]

const tau_pi_c = 2.6e-8        #lifetime of a charged pion [s]
const tau_pi_0 = 8.4e-17       #lifetime of a neutrion pion [s]

const sigma_T = 6.6524e-25   #Thompson cross section [cm^2]


#Unit Conversions
const pc_2_cm = 3.08568e18;       #number of centimeters in a parsec
const Mpc_2_cm = 1e6*pc_2_cm;
const Gpc_2_cm = 1e9*pc_2_cm;     #number of centimeters in a Gigaparsec
const yrs_2_s = 3.156e7;          #number of seconds in a year
const mb_2_sqr_cm = 1e-27;        #number of squared centimeters in a millibarn
const erg_2_GeV = 624.15;         #number of GeV in an erg
const erg_2_eV = 6.24150934e11;  #number of eV in an erg

const a_rdc   = 7.5646e-15
const G_cgs   = 6.6725985e-8
const k_B_eV  = 8.617332427e-5
