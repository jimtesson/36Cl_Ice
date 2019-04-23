function [Param_cosmo_out] = clrock(Data_in,Param_site_in,Const_cosmo,Sf)
% This function computes the variables that will be used to calculate 
%   36Cl concentrations. All those variables are passed to the structure 
%   'Param_cosmo'.
%   The first part computes variables for the Schimmelfenning 2009
%   model, based on exponential attenuation length approximation. 
%   The second part computes variables used in the LSD numerical resolution.
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation           
%           Param_site : Site specific parameters.
%           Data : vector containing the data provided by the user.
%                  (data must respect the CREp format)
%           
% Output : 
%            Param_cosmo : variables to calculate 36Cl concentration 
%            concentration 
%
% version 01/08/2018, written by TESSON J.
[~,N_samples] = size(Data_in);

Param_cosmo_out = cell(N_samples);

%% LSD calculation (Sato/Heisinger)
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later

load('pmag_consts.mat'); % constants for LSD calc.
consts = pmag_consts;

%   Trajectory-traced dipolar estimate for these purposes
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

RcEst = (dd(1)*cos(d2r(Param_site_in{1}.lat)) + ...
       dd(2)*(cos(d2r(Param_site_in{1}.lat))).^2 + ...
       dd(3)*(cos(d2r(Param_site_in{1}.lat))).^3 + ...
       dd(4)*(cos(d2r(Param_site_in{1}.lat))).^4 + ...
       dd(5)*(cos(d2r(Param_site_in{1}.lat))).^5 + ...
       dd(6)*(cos(d2r(Param_site_in{1}.lat))).^6);
   
Param_cosmo.maxdepth=20000; % maximum depth for calculation (g/cm^2) (=2500 in Chronus)
depthvector=[0:1:50 60:10:Param_cosmo.maxdepth];  % depth vector
Pressure=1013.25*exp(-0.03417/0.0065*(log(288.15)-log(288.15-0.0065*Param_site_in{1}.alt)));

%Muons are time-independent as coded.
flux_muon=muonfluxsato_36(depthvector,Pressure,RcEst,consts.SPhiInf,consts,'yes');

for i=1:N_samples % loop over samples
    
    Data = Data_in{i};
    Param_site = Param_site_in{i};
    % Assign variables
    data_target = Data.target; % target chemistry
    uncert_target = Data.uncert_target; % uncertainties on target chemistry
    data_bulk = Data.bulk; % bulk chemistry

    n = size(data_target,2) ;
    if n ~= 66, error('Sample file (target) must have 66 columns'), end
    n = size(data_bulk,2) ;
    if n ~= 66, error('Sample file (bulk) must have 66 columns'), end

    chimie_targ = data_target(1:n-4) ; 
    chimie_bulk = data_bulk(1:n-4) ;

    e = Param_site.mass_depth; % mass depth (g.cm-2)
    rho_rock = Param_site.rho_rock; % rock density ((g.cm-3)

    % Muons coefficients :
    Lambda_e = Const_cosmo.Lambda_f_e ; % g.cm-2
    Lambda_mu = Const_cosmo.Lambda_mu ; % g.cm-2

    % scaling factor for nucleonic production as a function of elevation, latitude (and temporal variations)
    S_EL_f = Sf{i}.SF_St_sp;

    % scaling factor for muonic production as a function of elevation, latitude (and temporal variations)
    S_EL_mu = Sf{i}.SF_St_mu;


% CHEMICAL ELEMENTS
%
% from 1 to 10  : As Ba Be Bi Cd Ce Co Cr Cs Cu
% from 11 to 20 : Dy Er Eu Ga Gd Ge Hf Ho In La
% from 21 to 30 : Lu Mo Nb Nd Ni Pb Pr Rb Sb Sm
% from 31 to 40 : Sn Sr Ta Tb Th Tm U  V  W  Y
% from 41 to 50 : Yb Zn Zr SiO2(Si) Al2O3(Al) Fe2O3(Fe) MnO(Mn) MgO(Mg) CaO(Ca) Na2O(Na)
% from 51 to 61 : K2O(K) TiO2(Ti) P2O5(P) B Li H2Otot(H) Stot(S) CO2tot(C) O_rock O_water CltotalAMS
% 62 : [Ca] in ppm from ICP

% A_k = atomic mass of element k (g.mol-1)
A_k = [74.9 137.327 9.012182 209.0 112.4 140.1 58.9332 51.9961 132.90545 63.5] ;
A_k = [A_k 162.5 167.3 152.0 69.7 157.25 72.6 178.5 164.9 114.8 138.9] ;
A_k = [A_k 175.0 95.94 92.9 144.2 58.6934 207.2 140.9 85.4678 121.8 150.36] ;
A_k = [A_k 118.7 87.62 180.9 158.9 232.0377 168.9 238.02891 50.9 183.8 88.9] ;
A_k = [A_k 173.0 65.4 91.224 28.085 26.981538 55.845 54.93804 24.305 40.078 22.98977] ;
A_k = [A_k 39.0983 47.867 30.973761 10.811 6.941 1.008 32.065 12.01 15.999 15.999 35.453] ;

% Conversion of oxyde percents into percents of the oxyded element in the
% target
% (Elements are given directly in ppm)
ppm_targ = chimie_targ ;
ppm_targ(44) = chimie_targ(44)*A_k(44)/(A_k(44) + 2*A_k(59)) ; % Si in percent
ppm_targ(45) = chimie_targ(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ; % Al in percent
ppm_targ(46) = chimie_targ(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ; % Fe in percent
ppm_targ(47) = chimie_targ(47)*A_k(47)/(A_k(47) + A_k(59)) ; % Mn in percent
ppm_targ(48) = chimie_targ(48)*A_k(48)/(A_k(48) + A_k(59)) ; % Mg in percent
ppm_targ(49) = chimie_targ(49)*A_k(49)/(A_k(49) + A_k(59)) ; % Ca in percent
ppm_targ(50) = chimie_targ(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ; % Na in percent
ppm_targ(51) = chimie_targ(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ; % K in percent
ppm_targ(52) = chimie_targ(52)*A_k(52)/(A_k(52) + 2*A_k(59)) ; % Ti in percent
ppm_targ(53) = chimie_targ(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) ; % P in percent
ppm_targ(56) = chimie_targ(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) ; % H water in percent
O_water = chimie_targ(56)*A_k(59)/(2*A_k(56) + A_k(59)) ; % O_water in percent
ppm_targ(58) = chimie_targ(58)*A_k(58)/(A_k(58) + 2*A_k(59)) ; % C in percent

ppm_targ(59) = sum([chimie_targ(44:53) chimie_targ(58)]) - sum([ppm_targ(44:53) ppm_targ(58)]) ; % O rock in percent
ppm_targ(60) = O_water ;
ppm_targ(44:53) = ppm_targ(44:53)*1e+4 ; % in ppm
ppm_targ(56) = ppm_targ(56)*1e+4 ; % in ppm
ppm_targ(58) = ppm_targ(58)*1e+4 ; % in ppm
ppm_targ(59) = ppm_targ(59)*1e+4 ; % in ppm
ppm_targ(60) = ppm_targ(60)*1e+4 ; % in ppm

ppm_CaO_uncer = uncert_target(1)*A_k(49)/(A_k(49) + A_k(59))*1e+4 ; % uncertainties on CaO (ppm)
ppm_K2O_uncer = uncert_target(2)*2*A_k(51)/(2*A_k(51) + A_k(59))*1e+4 ; % uncertainties on K2O (ppm)
ppm_TiO2_uncer = uncert_target(3)*A_k(52)/(A_k(52) + 2*A_k(59))*1e+4 ; % uncertainties on TiO2 (ppm)
ppm_Fe2O3_uncer = uncert_target(4)*2*A_k(46)/(2*A_k(46) + 3*A_k(59))*1e+4 ; % uncertainties on  Fe2O3 (ppm)

% Conversion of oxyde percents into percents of the oxyded element in the
% target (Elements are given directly in ppm)
ppm_bulk = chimie_bulk ;
ppm_bulk(44) = chimie_bulk(44)*A_k(44)/(A_k(44) + 2*A_k(59)) ; % Si in percent
ppm_bulk(45) = chimie_bulk(45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ; % Al in percent
ppm_bulk(46) = chimie_bulk(46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ; % Fe in percent
ppm_bulk(47) = chimie_bulk(47)*A_k(47)/(A_k(47) + A_k(59)) ; % Mn in percent
ppm_bulk(48) = chimie_bulk(48)*A_k(48)/(A_k(48) + A_k(59)) ; % Mg in percent
ppm_bulk(49) = chimie_bulk(49)*A_k(49)/(A_k(49) + A_k(59)) ; % Ca in percent
ppm_bulk(50) = chimie_bulk(50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ; % Na in percent
ppm_bulk(51) = chimie_bulk(51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ; % K in percent
ppm_bulk(52) = chimie_bulk(52)*A_k(52)/(A_k(52) + 2*A_k(59)) ; % Ti in percent
ppm_bulk(53) = chimie_bulk(53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) ; % P in percent
ppm_bulk(56) = chimie_bulk(56)*2*A_k(56)/(2*A_k(56) + A_k(59)) ; % H water in percent
O_water = chimie_bulk(56)*A_k(59)/(2*A_k(56) + A_k(59)) ; % O_water in percent (H20 - H)
ppm_bulk(58) = chimie_bulk(58)*A_k(58)/(A_k(58) + 2*A_k(59)) ; % C in percent

ppm_bulk(59) = sum([chimie_bulk(44:53) chimie_bulk(58)]) - sum([ppm_bulk(44:53) ppm_bulk(58)]) ; % O rock in percent
ppm_bulk(60) = O_water ;
ppm_bulk(44:53) = ppm_bulk(44:53)*1e+4 ; % in ppm
ppm_bulk(56) = ppm_bulk(56)*1e+4 ; % in ppm
ppm_bulk(58) = ppm_bulk(58)*1e+4 ; % in ppm
ppm_bulk(59) = ppm_bulk(59)*1e+4 ; % in ppm
ppm_bulk(60) = ppm_bulk(60)*1e+4 ; % in ppm

% Num_k = Atomic number of element k
Num_k = [33 56 4 83 48 58 27 24 55 29] ;
Num_k = [Num_k 66 68 63 31 64 32 72 67 49 57] ;
Num_k = [Num_k 71 42 41 60 28 82 59 37 51 62] ;
Num_k = [Num_k 50 38 73 65 90 69 92 23 74 39] ;
Num_k = [Num_k 70 30 40 14 13 26 25 12 20 11] ;
Num_k = [Num_k 19 22 15 5 3 1 16 6 8 8 17] ;

% Xi_k = average log-decrement of energy loss per collision for element k
Xi_k = [0 0 0 0 0 0 0 0.038 0 0] ;
Xi_k = [Xi_k 0 0 0 0 0.013 0 0 0 0 0] ;
Xi_k = [Xi_k 0 0 0 0 0 0 0 0 0 0.013] ;
Xi_k = [Xi_k 0 0 0 0 0 0.008594 0.008379 0 0 0] ;
Xi_k = [Xi_k 0 0 0 0.07 0.072 0.035 0.036 0.08 0.049 0.084] ;
Xi_k = [Xi_k 0.05 0.041 0.06321 0.174 0.264 1 0 0.15776 0.12 0.12 0.055] ;

% sigma_sc_k = neutron scattering x-section of element k (barns)

sigma_sc_k = [0 0 0 0 0 0 0 3.38 0 0] ;
sigma_sc_k = [sigma_sc_k 0 0 0 0 172 0 0 0 0 0] ;
sigma_sc_k = [sigma_sc_k 0 0 0 0 0 0 0 0 0 38] ;
sigma_sc_k = [sigma_sc_k 0 0 0 0 13.55 0 9.08 0 0 0] ;
sigma_sc_k = [sigma_sc_k 0 0 0 2.04 1.41 11.35 2.06 3.414 2.93 3.038] ;
sigma_sc_k = [sigma_sc_k 2.04 4.09 3.134 4.27 0.95 20.5 0 4.74 3.76 3.76 15.8] ;

% sigma_th_k = thermal neutron absorbtion x-section of element k (barns)
sigma_th_k = [0 0 0 0 0 0 0 3.1 0 0] ;
sigma_th_k = [sigma_th_k 0 0 0 0 41560 0 0 0 0 0] ;
sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 9640] ;
sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 0] ;
sigma_th_k = [sigma_th_k 0 0 0 0.17 0.23 2.56 13.3 0.063 0.43 0.53] ;
sigma_th_k = [sigma_th_k 2.15 6.1 0.2 767 70.5 0.33 0 0.0034 0.0002 0 33.5] ;
    
% I_a_k = dilute resonance integral for absorption of epithermal neutrons by element k (barns)
I_a_k = [0 0 0 390 0 0 0 1.6 0 0] ;
I_a_k = [I_a_k 0 0 0 0 390 0 0 0 0 0] ;
I_a_k = [I_a_k 0 0 0 0 0 0 0 0 0 1400] ;
I_a_k = [I_a_k 0 0 0 0 83.3 0 277 0 0 0] ;
I_a_k = [I_a_k 0 0 0 0.082 0.17 1.36 13.4 0.038 0.233 0.311] ;
I_a_k = [I_a_k 1 3.1 0.079 343 0 0 0 0.0018 0.000269 0.000269 13.83] ;

% f_d_k = proportion of muons stopped in element k that are captured by the nucleus
f_d_k = [0 0 0 0 0 0 0 0 0 0] ;
f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
f_d_k = [f_d_k 0 0 0 0.671 0.582 0.906 0 0.538 0.864 0.432] ;
f_d_k = [f_d_k 0.83 0 0 0 0 0 0 0.09 0.223 0 0] ; 

% Y_n = average neutron yield per captured muon
Y_n = [0 0 0 0 0 0 0 0 0 0] ;
Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
Y_n = [Y_n 0 0 0 0.86 1.26 1.125 0 0.6 0.75 1] ;
Y_n = [Y_n 1.25 0 0 0 0 0 0 0.76 0.8 0 0] ;

% S_i = mass stopping power (MeV/(g.cm-2))
S_i = [0 0 0.000529 0 0 0 0 0 0 0] ;
S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
S_i = [S_i 0 0 0 0.000454 0.000444 0.000351 0 0.000461 0.000428 0.000456] ;
S_i = [S_i 0.000414 0.000375 0.000433 0.000527 0.000548 0 0.000439 0.000561 0.000527 0.000527 0] ;

% Y_U_n = neutron yield (n/an/g/ppm de U)
Y_U_n = [0 0 265 0 0 0 0 0 0 0] ;
Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
Y_U_n = [Y_U_n 0 0 0 0.69 5.1 0.19 0 5.8 0 14.5] ;
Y_U_n = [Y_U_n 0.45 0 0 62.3 21.1 0 0 0.45 0.23 0.23 0] ;

% Y_TH_n = neutron yield (n/an/g/ppm de Th)
Y_Th_n = [0 0 91.2 0 0 0 0 0 0 0] ;
Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
Y_Th_n = [Y_Th_n 0 0 0 0.335 2.6 0.205 0 2.6 0 6.8] ;
Y_Th_n = [Y_Th_n 0.305 0 0 19.2 9.6 0 0 0.18 0.079 0.079 0] ;

Avogadro = 6.022E+23 ; % Avogadro Number

N_Cl_targ = (ppm_targ(61)./A_k(61))*Avogadro*1e-6 ; % Concentrations in atom/g

N_k_targ = (ppm_targ(:,1:61)./A_k)*Avogadro*1e-6 ; % Concentrations in atom/g
N_k_targ(56) = N_k_targ(56)/rho_rock ; % divided by bulk-rock density according to CHLOE for H

N_k_bulk = (ppm_bulk(:,1:61)./A_k)*Avogadro*1e-6 ; % Concentrations in atom/g
N_k_bulk(56) = N_k_bulk(56)/rho_rock ; % divided by bulk-rock density according to CHLOE for H

% -------------------------------- PRODUCTION RATES --------------------------------

% ------------------------------------ Spallation ------------------------------------ 

Psi_Cl36_Ca_0 = Const_cosmo.Psi_Cl36_Ca_0 ;% spallation production rate for Ca, SLHL
Psi_Cl36_Ca_0_uncert = Const_cosmo.Psi_Cl36_Ca_0_uncert ;

Psi_Cl36_K_0 = Const_cosmo.Psi_Cl36_K_0 ;% Spallation production rate at surface of 39K
Psi_Cl36_K_0_uncert = Const_cosmo.Psi_Cl36_K_0_uncert ;% Spallation production rate at surface of 39K

Psi_Cl36_Ti_0 = Const_cosmo.Psi_Cl36_Ti_0; % Spallation production rate at surface of Ti
Psi_Cl36_Ti_0_uncert = Const_cosmo.Psi_Cl36_Ti_0_uncert; % Spallation production rate at surface of Ti

Psi_Cl36_Fe_0 = Const_cosmo.Psi_Cl36_Fe_0; % Spallation production rate at surface of Fe
Psi_Cl36_Fe_0_uncert = Const_cosmo.Psi_Cl36_Fe_0_uncert; % Spallation production rate at surface of Fe

% Psi_Cl36_Ca_0 ; % Spallation production rate at surface of 40Ca
% (at of Cl36 /g of Ca per yr) [48.8 ? 3.4 Stone et al. 1996 and Evans et al. 1997] 
% Stone 2000: 48.8 +/- 3.5; Dunai 2001: 53.7 +/- 3.9; Pigati and Lifton 2004 (Desilets and Zreda, 2003): 53.1 +/- 3.8
% Lifton et al., 2005: 59.4 +/- 4.3; Pigati and Lifton 2004 (Desilets et al., 2006): 54.7 +/- 4.0
% Lifton et al., 2008: 58.9 +/- 4.3
C_Ca = ppm_targ(62)*1e-6;  % Mass concentration of Ca (g of Ca per g of rock) % from ICP
P_sp_Ca = Psi_Cl36_Ca_0*C_Ca;  % result unscaled 36Cl production by spallation of 40Ca (atoms 36Cl g-1 yr-1)

C_K = ppm_targ(51)*1e-6 ; % Mass concentration of K (g of K per g of rock)
P_sp_K = Psi_Cl36_K_0*C_K ; % result unscaled 36Cl production by spallation of 39K (atoms 36Cl g-1 yr-1)

C_Ti = ppm_targ(52)*1e-6 ; % Mass concentration of Ti (g of Ti per g of rock)
P_sp_Ti = Psi_Cl36_Ti_0*C_Ti ; % result unscaled 36Cl production by spallation of Ti (atoms 36Cl g-1 yr-1)


C_Fe = ppm_targ(46)*1e-6 ; % Mass concentration of Fe (g of Fe per g of rock)
P_sp_Fe = Psi_Cl36_Fe_0*C_Fe ; % result unscaled 36Cl production by spallation of Fe (atoms 36Cl g-1 yr-1)

P_sp = (P_sp_Ca + P_sp_K + P_sp_Ti + P_sp_Fe)*exp(-e/Lambda_e) ; % Unscaled Spallation production rate (atoms 36Cl g-1 yr-1)

% Uncertainties
   % on Ca
    C_Ca_uncer = ppm_CaO_uncer .* 1e-6;
    P_sp_Ca_uncert = ((Psi_Cl36_Ca_0.*C_Ca_uncer)^2+(C_Ca.*Psi_Cl36_Ca_0_uncert)^2)^.5;
   % on K
    C_K_uncer = ppm_K2O_uncer .* 1e-6;
    P_sp_K_uncert = ((Psi_Cl36_K_0.*C_K_uncer)^2+(C_K.*Psi_Cl36_K_0_uncert)^2)^.5;
   % on Ti 
    C_Ti_uncer = ppm_TiO2_uncer .* 1e-6;
    P_sp_Ti_uncert = ((Psi_Cl36_Ti_0.*C_Ti_uncer)^2+(C_Ti.*Psi_Cl36_Ti_0_uncert)^2)^.5;
   % on Fe 
    C_Fe_uncer = ppm_Fe2O3_uncer .* 1e-6;
    P_sp_Fe_uncert = ((Psi_Cl36_Fe_0.*C_Fe_uncer)^2+(C_Fe.*Psi_Cl36_Fe_0_uncert)^2)^.5;

% -------------------- Direct capture of slow negative muons ---------------------
% -------------------- by target elements Ca and K ------------------------------- 

%f_n_K = 0.02 ; % Fabryka-Martin (1988)
%f_n_Ca = 0.062 ; % Fabryka-Martin (1988)
f_n_Ca = 0.045 ;  % +/- 0.005 Heisinger et al. (2002)
f_n_K = 0.035 ; % +/- 0.005 Heisinger et al. (2002)
f_i_Ca = 0.969 ; % Fabryka-Martin (1988)
f_i_K = 0.933 ; % Fabryka-Martin (1988)
f_d_Ca = 0.864 ; % Fabryka-Martin (1988)
f_d_K = 0.83 ; % Fabryka-Martin (1988)
Psi_mu_0 = Const_cosmo.Psi_mu_0; % slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

MjZj_bulk = Num_k.*ppm_bulk(:,1:61)./A_k.*1e-6;

MjZj_bulk_sum = sum(MjZj_bulk) ;
MjZj_targ_sum = sum(Num_k.*ppm_targ(:,1:61)./A_k)*1e-6 ;

f_c_k_bulk = MjZj_bulk./MjZj_bulk_sum;

f_c_Ca = (Num_k(49)*ppm_targ(62)*1e-6/A_k(49))/(sum(Num_k.*ppm_targ(:,1:61)./A_k)*1e-6) ; % for Ca (ICP)
f_c_K = (Num_k(51)*ppm_targ(51)*1e-6/A_k(51))/(sum(Num_k.*ppm_targ(:,1:61)./A_k)*1e-6) ; % for K

Y_Sigma_Ca = f_c_Ca*f_i_Ca*f_d_Ca*f_n_Ca ; % 36Cl production per stopped muon 
% Y_Sigma_Ca DEPENDS ON CHEMICAL COMPOSITION
Y_Sigma_K = f_c_K*f_i_K*f_d_K*f_n_K ; % 36Cl production per stopped muon 
% Y_Sigma_K DEPENDS ON CHEMICAL COMPOSITION

Y_Sigma = Y_Sigma_Ca + Y_Sigma_K ;

P_mu = Y_Sigma*Psi_mu_0*exp(-e/Lambda_mu);  % Unscaled slow negative muon production rate (atoms 36Cl g-1 yr-1)


% ------------------------------------ Epithermal neutrons ------------------------------------ 

B = sum(Xi_k.*sigma_sc_k.*N_k_bulk)*1e-24 ; % Scattering rate parameter
% B DEPENDS ON CHEMICAL COMPOSITION

I_eff = sum(I_a_k.*N_k_bulk)*1e-24 ; % (Eq 3.9, Gosse & Phillips, 2001)
% Effective macroscopic resonance integral for absorbtion of epith neutrons (cm2.g-1)
% I_eff DEPENDS ON CHEMICAL COMPOSITION

f_eth = N_Cl_targ*I_a_k(61)*(1e-24)/I_eff ; % (Eq 3.17, Gosse & Phillips, 2001)
% Fraction of epith neutrons absorbed by Cl35
% f_eth DEPENDS ON CHEMICAL COMPOSITION

p_E_th = exp(-I_eff/B) ; % (Eq 3.8, Gosse & Phillips, 2001)
% Resonance escape probability of a neutron from the epith energy range in subsurface
% p_E_th DEPENDS ON CHEMICAL COMPOSITION

A = sum(A_k.*N_k_bulk)/sum(N_k_bulk) ;
% Average atomic weight (g/mol)

A_a = 14.5 ; % Average atomic weight of air

R_eth = sqrt(A/A_a) ; % (Eq 3.24, Gosse & Phillips, 2001)
% Ratio of epithermal neutron production in subsurface to that in atm
% R_eth DEPENDS ON CHEMICAL COMPOSITION
R_eth_a = 1 ;

Sigma_sc = sum(sigma_sc_k.*N_k_bulk)*1e-24 ; % (Eq 3.22, Gosse & Phillips, 2001)
% Macroscopic neutron scattering cross-section (cm2.g-1)
% Sigma_sc DEPENDS ON CHEMICAL COMPOSITION

Sigma_sc_a = 0.3773 ;% macroscopic neutron scaterring cross section of the atmosphere (cm2.g-1)

Xi = B/Sigma_sc ;  % Eq 3.19 Goss and Phillips
% Average log decrement energy loss per neutron collision
% Xi DEPENDS ON CHEMICAL COMPOSITION

Sigma_eth = Xi*(I_eff + Sigma_sc) ; % (Eq 3.18, Gosse & Phillips, 2001)
% Effective epithermal loss cross-section (cm2.g-1)
% Sigma_eth DEPENDS ON CHEMICAL COMPOSITION

Lambda_eth = 1/Sigma_eth ; % (Eq 3.18,Gosse & Phillips, 2001)
% Attenuation length for absorbtion and moderation of epith neutrons flux (g.cm-2)
% Lambda_eth DEPENDS ON CHEMICAL COMPOSITION

D_eth = 1/(3*Sigma_sc*(1 - 2/(3*A))) ; % (Eq 3.21, Gosse & Phillips, 2001)
% Epithermal neutron diffusion coefficient (g.cm-2)
% D_eth DEPENDS ON CHEMICAL COMPOSITION

D_eth_a = 1/(3*Sigma_sc_a*(1 - 2/(3*A_a))) ; % (Eq 3.21, Gosse & Phillips, 2001)
% Epithermal neutron diffusion coefficient in atmosphere (g.cm-2)

P_f_0 = Const_cosmo.P_f_0 ; % Production rate of epithermal neutrons from fast neutrons in atm at land/atm interface (n cm-2 yr-1), Gosse & Philipps, 2001.
P_f_0_uncert = Const_cosmo.P_f_0_uncert;


phi_star_eth = P_f_0*R_eth/(Sigma_eth - (D_eth/(Lambda_e^2))) ; % Epithermal neutron flux at land/atmosphere
% interface that would be observed in ss if interface was not present (n cm-2 yr-1)

%Y_s = sum(f_d_k.*Y_n.*ppm_bulk(:,1:61).*N_k_bulk./A_k)/sum(ppm_bulk(:,1:61).*N_k_bulk./A_k) ;
Y_s = sum(f_c_k_bulk.*f_d_k.*Y_n);
% Average neutron yield per stopped negative muon
% Y_s DEPENDS ON CHEMICAL COMPOSITION

Sigma_eth_a = 0.0548 ; % Macroscopic absorption and moderation x-section in atm. (cm2 g-1) - Constant (Chloe)
D_th_a = 0.9260472 ; % Thermal neutron diffusion coeff in atm. (g*cm-2) - Constant (Chloe)
Sigma_sc_a = 0.3773 ; % Macroscopic neutron scaterring cross section of atmosphere (cm2.g-1) - Constant (Chloe)

phi_star_eth_a = P_f_0*R_eth_a/(Sigma_eth_a - (D_eth_a/(Lambda_e^2))) ; % Epithermal neutron flux at land/atmosphere interface that would be observed in atm 

% if interface was not present (n cm-2 yr-1)
phi_mu_f_0 = Const_cosmo.phi_mu_f_0 ; % Fast muon flux at land surface, sea level, high latitude, Gosse & Phillips, 2001 (? cm-2 yr-1)
P_n_mu_0 = (Y_s*Psi_mu_0 + 5.8e-6*phi_mu_f_0); % Fast muon flux at land surface SLHL, Eq.3.49 Gosse & Phillips, 2001 (n cm-2 yr-1)
R_mu = S_EL_mu*P_n_mu_0/(S_EL_f*P_f_0*R_eth) ; %Ratio of muon production rate to epithermal neutron production rate

Deltaphi_2star_eth_a = phi_star_eth - D_eth_a*phi_star_eth_a/D_eth ; % Adjusted difference between hypothetical equilibrium epithermal neutron fluxes in atm and ss (n cm-2 yr-1)

L_eth = 1/sqrt(3*Sigma_sc*Sigma_eth); % Epithermal neutron diffusion length (g cm-2)
% L_eth DEPENDS ON CHEMICAL COMPOSITION

L_eth_a = 1/sqrt(3*Sigma_sc_a*Sigma_eth_a); % Epithermal neutron diffusion length in atm (g cm-2)
FDeltaphi_star_eth = ((D_eth_a/L_eth_a)*(phi_star_eth_a - phi_star_eth) - ...
    Deltaphi_2star_eth_a*(D_eth/Lambda_e))/...
    ((D_eth_a/L_eth_a) + (D_eth/L_eth)) ; % EQ. 3.28 Gosse & Phillips, 2001

% Difference between phi_star_eth,ss and actual epithermal neutron flux at land surface

phi_eth_total = phi_star_eth*exp(-e/Lambda_e) + ...
    (1 + R_mu*R_eth)*FDeltaphi_star_eth*exp(-e/L_eth) + ...
    R_mu*phi_star_eth*exp(-e/Lambda_mu) ; % Epithermal neutron flux (concentration) (n cm-2 yr-1)

P_eth = (f_eth/Lambda_eth)*phi_eth_total*(1 - p_E_th) ;

% ------------------------------------ Thermal neutrons ------------------------------------ 

Sigma_th = sum(N_k_bulk.*sigma_th_k)*1e-24 ; % Eq 3.6 de Gosse and Phillips, 2001
% macroscopic thermal neutron absorbtion cross-section 
% Sigma_th DEPENDS ON CHEMICAL COMPOSITION

f_th = sigma_th_k(61)*N_Cl_targ*1e-24/Sigma_th ; % Eq 3.32 de Gosse and Phillips, 2001
% fraction of thermal neutrons absorbed by Cl35
% f_th DEPENDS ON CHEMICAL COMPOSITION

Lambda_th = 1/Sigma_th ; % Eq 3.35 Gosse anf Phillips, 2001
% Attenuation length for absorbtion of thermal neutrons flux (g.cm-2)
% Lambda_th DEPENDS ON CHEMICAL COMPOSITION

p_E_th_a = 0.56 ; % Resonance escape probability of the atmosphere - Constant (Chloe)
R_th = p_E_th/p_E_th_a ; % Ratio of thermal neutron production in ss to that in atm ; Eq 3.34 Gosse and Phillips, 2001
DD_th = D_eth ; % D_th = 2.99
R_th_a = 1 ;
Deltaphi_star_eth_a = phi_star_eth - phi_star_eth_a ; % difference in equilibrium epithermal neutron fluxes between atm and ss

FDeltaphi_star_eth_a = (D_eth*Deltaphi_star_eth_a/L_eth - D_eth*Deltaphi_2star_eth_a/Lambda_e)/ ...
    (D_eth_a / L_eth_a + D_eth / L_eth );

Sigma_th_a = 0.060241 ; % Constant from Chloe - macroscopic thermal neutron cross section of atm (cm2 g-1)
phi_star_th = (p_E_th_a*R_th*phi_star_eth)/(Lambda_eth*(Sigma_th - DD_th/(Lambda_e^2))) ;

% thermal neutron flux at land/atm interface that would be observed in atm if interface not present (n.cm_2.a-1)
R_prime_mu = (p_E_th_a/p_E_th)*R_mu ; % ratio of muon production rate to thermal neutron production rate

JDeltaphi_star_eth = (p_E_th_a*R_th*FDeltaphi_star_eth)/(Lambda_eth*(Sigma_th - DD_th/(L_eth^2))) ; % Eq. 3.39 Gosse & Phillips, 2001

% Portion of difference between phi_star_eth,ss and actual flux due to epithermal flux profile
JDeltaphi_star_eth_a = (p_E_th_a*R_th_a*FDeltaphi_star_eth_a)/((1/Sigma_eth_a)*(Sigma_th_a - D_th_a/(L_eth_a^2))) ;

% Portion of difference between pi_star_eth,a and actual flux due to epithermal flux profile

L_th = sqrt(DD_th/Sigma_th) ;
L_th_a = sqrt(D_th_a/Sigma_th_a) ; % thermal neutron diffusion length in atm (g cm-2)
phi_star_th_a = (p_E_th_a*R_th_a*phi_star_eth_a)/(1/Sigma_eth_a*(Sigma_th_a - D_th_a/(Lambda_e^2))) ; 

% thermal neutron flux at land/atmosphere interface that would be observed in atm if interface was not present (n cm-2 yr-1)
Deltaphi_star_th = phi_star_th_a - phi_star_th ; % difference between hypothetical equilibrium thermal neutron fluxes in atmosphere and ss

JDeltaphi_star_th = (D_th_a*(phi_star_th_a/Lambda_e - JDeltaphi_star_eth_a/L_eth_a) ...
   - DD_th*(phi_star_th/Lambda_e + JDeltaphi_star_eth/L_eth) ...
   + (D_th_a/L_th_a)*(Deltaphi_star_th + JDeltaphi_star_eth_a - JDeltaphi_star_eth)) ...
   / ((DD_th/L_th) + (D_th_a/L_th_a)) ; % portion of difference between phi_star_th,ss and actual flux due to thermal flux profile

phi_th_total = phi_star_th *exp(-e/Lambda_e) + ...
                + JDeltaphi_star_eth *exp(-e/L_eth) ...
                + R_prime_mu*JDeltaphi_star_eth*exp(-e/L_eth) ...
                + JDeltaphi_star_th *exp(-e/L_th) ... 
                + R_prime_mu*R_th*JDeltaphi_star_th*exp(-e/L_th)... 
                + phi_star_th*R_prime_mu*exp(-e/Lambda_mu) ; % Thermal neutron flux (n.cm_2.a-1)


P_th = (f_th/Lambda_th)*phi_th_total ; % Result unscaled sample specific 36Cl production rate by capture of thermal neutrons (atoms 36Cl g-1 yr-1)


% ------------------------------------ Radiogenic production -----------------------------------------

X = (sum(ppm_bulk(:,1:61).*S_i.*Y_U_n))/(sum(S_i.*ppm_bulk(:,1:61))) ;
% X DEPENDS ON CHEMICAL COMPOSITION

Y = (sum(ppm_bulk(:,1:61).*S_i.*Y_Th_n))/(sum(S_i.*ppm_bulk(:,1:61))) ;
% Y DEPENDS ON CHEMICAL COMPOSITION

U = ppm_bulk(37) ; % Concentration en Uranium (ppm)
Th = ppm_bulk(35) ; % Concentration en Thorium (ppm)

P_n_alphan = X*U + Y*Th ; % alpha,n reactions
P_n_sf = 0.429*U ; % spontaneous fission
P_th_r = (P_n_alphan + P_n_sf)*p_E_th ; % total radiogenic thermal neutron production
P_eth_r = (P_n_alphan + P_n_sf)*(1 - p_E_th) ; % total radiogenic epithermal neutron production

P_rad = P_th_r*f_th + P_eth_r*f_eth ;


% Shielding factors

S_L_th = 1 ; % diffusion out of objects (poorly constrained)
S_L_eth = 1 ; % diffusion out of objects (poorly constrained)


%% TOTAL PRODUCTION RATE at surface

% Cosmogenic production:
%P_cosmo = so_e*S_EL_f*(Q_sp.*P_sp + S_L_th*Q_th*P_th + S_L_eth*Q_eth*P_eth) + so_mu*S_EL_mu*Q_mu.*P_mu ;
P_cosmo = S_EL_f*(P_sp + S_L_th*P_th + S_L_eth*P_eth) + S_EL_mu*P_mu ;

%%% scaled sources of production %%%
P = P_cosmo+P_rad;

%% Inherited 36Cl
Param_cosmo.N36Cl.inh = Data.N36Cl.inh_0 * exp(-Const_cosmo.lambda36*Param_site.t_expo_estim); 

%% Muon calculation (Sato/Heisinger)
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later

%store the output fluxes that we need
Param_cosmo.negflux=flux_muon.R;
Param_cosmo.totalflux=flux_muon.phi;

%Also store the muons production rates from the code
Param_cosmo.muon36(1,:)=depthvector;
Param_cosmo.depthvector=depthvector;

%calculating fast muon contribution to chlorine (individually for Ca and K)
z=depthvector;

% fast muon production
% Sigma0 Ca and K
Param_cosmo.consts.sigma0_Ca = 7.3e-30; % Heisinger
Param_cosmo.consts.sigma0_K = 9.4e-30; % Marrero

aalpha = 1.0;
Beta = 1.0;
Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));

P_fast_K = flux_muon.phi .* Beta .* (Ebar.^aalpha) .* Param_cosmo.consts.sigma0_K .* N_k_targ(51);
P_fast_Ca = flux_muon.phi .* Beta .* (Ebar.^aalpha) .* Param_cosmo.consts.sigma0_Ca .* N_k_targ(49);
P_fast_total = P_fast_K + P_fast_Ca;

% negative muon capture
P_neg_K = flux_muon.R .* Y_Sigma_K; 
P_neg_Ca = flux_muon.R .* Y_Sigma_Ca;

%sum parts of the muon production
Param_cosmo.Prodmu=P_neg_K+P_neg_Ca+P_fast_total;
Param_cosmo.muon36(2,:)=P_neg_Ca+P_fast_Ca;
Param_cosmo.muon36(3,:)=P_neg_K+P_fast_K;
Param_cosmo.muon36(4,:)=P_neg_Ca;
Param_cosmo.muon36(5,:)=P_neg_K;
Param_cosmo.muon36(6,:)=P_fast_Ca;
Param_cosmo.muon36(7,:)=P_fast_K;


%% Passing parameters through the structure Param_cosmo.

Param_cosmo.P_sp_Ca = P_sp_Ca;
Param_cosmo.P_sp_K = P_sp_K;
Param_cosmo.P_sp_Ti = P_sp_Ti;
Param_cosmo.P_sp_Fe = P_sp_Fe;

Param_cosmo.P_cosmo = P_cosmo;
Param_cosmo.P = P;
Param_cosmo.P_sp = P_sp;
Param_cosmo.P_rad = P_rad;
Param_cosmo.P_th = P_th;
Param_cosmo.P_eth = P_eth;
Param_cosmo.P_mu = P_mu;

Param_cosmo.L_eth = L_eth; % thermal neutron diffusion length in subsurface
Param_cosmo.L_th = L_th; % epithermal neutron diffusion length in subsurface
Param_cosmo.L_sp = Lambda_e; % effective fast neutron attenuation coefficient
Param_cosmo.L_mu = Lambda_mu; % slow muon attenuation length

Param_cosmo.C_Ca = C_Ca;
Param_cosmo.C_K = C_K;
Param_cosmo.C_Ti = C_Ti;
Param_cosmo.C_Fe = C_Fe;

Param_cosmo.A_k = A_k;
Param_cosmo.I_eff = I_eff;
Param_cosmo.I_a_k = I_a_k;
Param_cosmo.sigma_th_k = sigma_th_k;
Param_cosmo.Sigma_th = Sigma_th;
Param_cosmo.p_E_th = p_E_th;
Param_cosmo.f_eth = f_eth;
Param_cosmo.Lambda_eth = Lambda_eth;
Param_cosmo.phi_star_eth = phi_star_eth;
Param_cosmo.f_th = f_th;
Param_cosmo.R_th = R_th;
Param_cosmo.Lambda_th = Lambda_th;
Param_cosmo.phi_star_th = phi_star_th;
Param_cosmo.FDeltaphi_star_eth = FDeltaphi_star_eth;
Param_cosmo.R_mu = R_mu;
Param_cosmo.R_eth = R_eth;
Param_cosmo.JDeltaphi_star_eth = JDeltaphi_star_eth;
Param_cosmo.JDeltaphi_star_th = JDeltaphi_star_th;
Param_cosmo.R_prime_mu = R_prime_mu;
Param_cosmo.Y_Sigma = Y_Sigma;
Param_cosmo.P_f_0 = P_f_0;
Param_cosmo.P_f_0_uncert = P_f_0_uncert;
Param_cosmo.p_E_th_a = p_E_th_a;
Param_cosmo.P_th_r = P_th_r ;
Param_cosmo.P_eth_r = P_eth_r;
Param_cosmo.Y_s = Y_s;
Param_cosmo.Sigma_eth = Sigma_eth;
Param_cosmo.D_eth = D_eth;
Param_cosmo.Sigma_eth_a = Sigma_eth_a;
Param_cosmo.D_eth_a = D_eth_a;
Param_cosmo.R_eth_a = R_eth_a;
Param_cosmo.D_eth_a = D_eth_a;
Param_cosmo.L_eth_a = L_eth_a;
Param_cosmo.Sigma_th = Sigma_th;
Param_cosmo.DD_th = DD_th;
Param_cosmo.R_th_a = R_th_a;
Param_cosmo.Sigma_th_a = Sigma_th_a;
Param_cosmo.D_th_a = D_th_a;
Param_cosmo.L_th_a = L_th_a;

Param_cosmo.ppm_targ = ppm_targ;
Param_cosmo.ppm_bulk = ppm_bulk;


% uncertainties
Param_cosmo.uncert.P_sp_Ca_uncert = P_sp_Ca_uncert;
Param_cosmo.uncert.P_sp_Fe_uncert = P_sp_Fe_uncert;
Param_cosmo.uncert.P_sp_Ti_uncert = P_sp_Ti_uncert;
Param_cosmo.uncert.P_sp_K_uncert = P_sp_K_uncert;

Param_cosmo_out{i} = Param_cosmo;
end

