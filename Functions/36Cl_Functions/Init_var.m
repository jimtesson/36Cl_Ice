function [Data,Param_site,Const_cosmo] = Init_var(DATA_CREP,flag)
% This function allows initialize the data, the cosmogenic constants, and 
%   site-specific variables.
% INPUT : 
%       DATA_CREP : sample data matrix from the user, using the CREp format.
% 
%           1: sample name, 2: latitude (deg), 
%           3: longitude (deg), 4: altitude (masl)
%           5: shielding, 6: density (g/cm3), 
%           7: tickness (cm), 8: erosion (cm/an)
%           9: formation age (age), 10: +/- formation age (a),

%           11: 36Cl conc (at/g), 12: +/- 36Cl conc (at/g),
%           13: Cl_ppm_target, 14: +/- Cl_ppm_target
%           15: SiO2_%_target, 
%           16: Al2O3_%_target, 
%           17: Fe2O3_%_target, 18:+/- Fe2O3_%_target
%           19: MnO_%_target,	
%           20: MgO_%_target, 
%           21: CaO_%_target, 22:	+/- CaO_%_target, 
%           23: Na2O_%_target, 
%           24: K2O_%_target, 25: +/- K2O_%_target, 
%           26: TiO2_%_target, 27: +/- TiO2_%_target, 
%           28: P2O5_%_target, 
%           29: LOI_%_target, 
%           30: H2O_%_target,
%           31: CO2_%_target.
%
%           32: SiO2_%_bulk, 
%           33: Al2O3_%_bulk, 
%           34: Fe2O3_%_bulk, 
%           35: MnO_%_bulk, 
%           36: MgO_%_bulk, 
%           37: CaO_%_bulk, 
%           38: Na2O_%_bulk, 
%           39: K2O_%_bulk, 
%           40: TiO2_%_bulk, 
%           41: P2O5_%_bulk,	
%           42: H2O_%_bulk,	
%           43: CO2_%_bulk, 
%           44: S_%_bulk, 
%           45: LOI_%_bulk, 
%           46: Li_ppm_bulk,
%           47: B_ppm_bulk, 
%           48: Cl_ppm_bulk,
%           49: Cr_ppm_bulk, 
%           50: Co_ppm_bulk, 
%           51: Sm_ppm_bulk,
%           52: Gd_ppm_bulk, 
%           53: Th_ppm_bulk, 54: +/- Th_ppm_bulk, 
%           55: U_ppm_bulk, 56: +/- U_ppm_bulk
%
%       Z : vector of the samples
%       EL : User-provided scaling factors for neutrons and muons.
%
% OUTPUT :
%       Data : Data matrix formatted following Schlagenhauf (2010)
%           from 1 to 10  : CaO K2O TiO2 Fe2O3 Cl SiO2 Na2O MgO Al2O3 MnO
%           from 11 to 20 : P2O5 CO2 H2O Li B Sm Gd Th U Cr
%           from 21 to 30 : Be N F S Sc Co Ni Se Br Rb
%           from 31 to 40 : Sr Zr Mo Ru Te I Cs Ba
%
%       Param_site : site-specific variables
%       Const_cosmo : cosmogenic constants
%

% number of samples
[~,N_samples] = size(DATA_CREP.Z);

% Constants
Const_cosmo.lambda36 = DATA_CREP.lambda36;
Const_cosmo.lambda36_uncert = DATA_CREP.lambda36_uncert;

% Attenuation length
Const_cosmo.Lambda_f_e = DATA_CREP.Lambda_f_e; % effective fast neutron attenuation coefficient (g.cm-2)
Const_cosmo.Lambda_mu = DATA_CREP.Lambda_mu ; % slow muon attenuation length (g.cm-2)
Const_cosmo.Psi_mu_0 = DATA_CREP.Psi_mu_0 ; % slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

% Unscaled sample specific 36Cl production rate by spallation of target elements
Const_cosmo.Psi_Cl36_Ca_0 = DATA_CREP.Psi_Cl36_Ca_0 ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_K_0 = DATA_CREP.Psi_Cl36_K_0 ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_Ti_0 = DATA_CREP.Psi_Cl36_Ti_0 ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_Fe_0 = DATA_CREP.Psi_Cl36_Fe_0 ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr)

% P_f_0: production rate of epithermal neutrons from fast neutrons in atmosphere at land/atm interface
%           P_f_0 is choosen following the scaling model (Lal/Stone, or
%           LSD) From Marrero 2015.
    if(strcmp(flag.scaling_model,'st'))
        Const_cosmo.P_f_0 = 691;  
        Const_cosmo.P_f_0_uncert = 186;
    elseif(strcmp(flag.scaling_model,'sa'))
        Const_cosmo.P_f_0 = 759; 
        Const_cosmo.P_f_0_uncert = 180;   
    end    
% phi_mu_f_0
Const_cosmo.phi_mu_f_0 = 7.9e+5 ; % Fast muon flux at land surface, sea level, high latitude, Gosse & Phillips, 2001 (? cm-2 yr-1)

% uncertainties
Const_cosmo.Psi_Cl36_Ca_0_uncert = DATA_CREP.Psi_Cl36_Ca_0_uncert ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_K_0_uncert = DATA_CREP.Psi_Cl36_K_0_uncert ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_Ti_0_uncert = DATA_CREP.Psi_Cl36_Ti_0_uncert ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
Const_cosmo.Psi_Cl36_Fe_0_uncert = DATA_CREP.Psi_Cl36_Fe_0_uncert ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr)


for i=1:N_samples

    clear site N36Cl tmp;
    data_bulk = zeros(1,66);
    data_target = zeros(1,66);
    uncert_target = zeros(1,5);  
% Chemical data
    
% DATA_CREP input data file format:
% 1: sample name, 2: latitude (deg), 3: longitude (deg), 4: altitude (masl)
% 5: shielding, 6: density (g/cm3), 7: tickness (cm), 8: erosion (cm/an)
% 9: formation age (age), 10: +/- formation age (a),

% 11: 36Cl conc (at/g), 12: +/- 36Cl conc (at/g),
% 13: Cl_ppm_target, 14: +/- Cl_ppm_target
% 15: SiO2_%_target, 
% 16: Al2O3_%_target, 
% 17: Fe2O3_%_target, 18:+/- Fe2O3_%_target
% 19: MnO_%_target,	
% 20: MgO_%_target, 
% 21: CaO_%_target, 22:	+/- CaO_%_target, 
% 23: Na2O_%_target, 
% 24: K2O_%_target, 25: +/- K2O_%_target, 
% 26: TiO2_%_target, 27: +/- TiO2_%_target, 
% 28: P2O5_%_target, 
% 29: LOI_%_target, 
% 30: H2O_%_target,
% 31: CO2_%_target.
%
% 32: SiO2_%_bulk, 
% 33: Al2O3_%_bulk, 
% 34: Fe2O3_%_bulk, 
% 35: MnO_%_bulk, 
% 36: MgO_%_bulk, 
% 37: CaO_%_bulk, 
% 38: Na2O_%_bulk, 
% 39: K2O_%_bulk, 
% 40: TiO2_%_bulk, 
% 41: P2O5_%_bulk,	
% 42: H2O_%_bulk,	
% 43: CO2_%_bulk, 
% 44: S_%_bulk, 
% 45: LOI_%_bulk, 
% 46: Li_ppm_bulk,
% 47: B_ppm_bulk, 
% 48: Cl_ppm_bulk,
% 49: Cr_ppm_bulk, 
% 50: Co_ppm_bulk, 
% 51: Sm_ppm_bulk,
% 52: Gd_ppm_bulk, 
% 53: Th_ppm_bulk, 54: +/- Th_ppm_bulk, 
% 55: U_ppm_bulk, 56: +/- U_ppm_bulk

% Fomarting data file: converting from the Schimmelfenning XCELL SPREADSHEET format to the
% Schlagenhauf matlab program format
%
% from 1 to 10  : CaO K2O TiO2 Fe2O3 Cl SiO2 Na2O MgO Al2O3 MnO
% from 11 to 20 : P2O5 CO2 H2O Li B Sm Gd Th U Cr
% from 21 to 30 : Be N F S Sc Co Ni Se Br Rb
% from 31 to 40 : Sr Zr Mo Ru Te I Cs Ba

% BULK CHEMICAL ELEMENTS

data_bulk(1,49) = DATA_CREP.CaO_bulk(i) ; % CaO (%)
data_bulk(1,51) = DATA_CREP.K2O_bulk(i) ; % K2O (%)
data_bulk(1,52) = DATA_CREP.TiO2_bulk(i) ; % TiO2 (%)
data_bulk(1,46) = DATA_CREP.Fe2O3_bulk(i) ; % Fe2O3 (%)
data_bulk(1,61) = DATA_CREP.Cl_bulk(i) ; % Cl (ppm)
data_bulk(1,44) = DATA_CREP.SiO2_bulk(i) ; % SiO2 (%)
data_bulk(1,50) = DATA_CREP.NaO_bulk(i) ; % Na2O (%)
data_bulk(1,48) = DATA_CREP.MgO_bulk(i) ; % MgO (%)
data_bulk(1,45) = DATA_CREP.Al2O3_bulk(i) ; % Al2O3 (%)
data_bulk(1,47) = DATA_CREP.MnO_bulk(i) ; % MnO (%)
data_bulk(1,53) = DATA_CREP.P2O5_bulk(i) ; % P2O5 (%)
data_bulk(1,58) = DATA_CREP.CO2_bulk(i) ; % CO2 (%)
data_bulk(1,56) = DATA_CREP.H2O_bulk(i) ; % H20 (%)
data_bulk(1,55) = DATA_CREP.Li_bulk(i) ; % Li (ppm)
data_bulk(1,54) = DATA_CREP.B_bulk(i) ; % B (ppm)
data_bulk(1,30) = DATA_CREP.Sm_bulk(i) ; % Sm (ppm)
data_bulk(1,15) = DATA_CREP.Gd_bulk(i) ; % Gd (ppm)
data_bulk(1,35) = DATA_CREP.Th_bulk(i) ; % Th (ppm)
data_bulk(1,37) = DATA_CREP.U_bulk(i) ; % U (ppm)
data_bulk(1,8) = DATA_CREP.Cr_bulk(i) ; % Cr (ppm)
%data_bulk(1,3) = DATA_BULK_SC(1,21) ; % Be (ppm) is not considered in the CREp program.
% data_bulk(1,) = DATA_BULK_SC(1,22) ; % N is not considered in the Schlagenhauf (2010) program. 
% data_bulk(1,) = DATA_BULK_SC(1,23) ; % F is not considered in the Schlagenhauf (2010) program.
data_bulk(1,57) = DATA_CREP.S_bulk(i) ; % S (in %)
%data_bulk(1,) = DATA_BULK_SC(1,25) ; % Sc is not considered in the Schlagenhauf (2010) program.
data_bulk(1,7) = DATA_CREP.Co_bulk(i) ; % Co (ppm)
%data_bulk(1,25) = DATA_BULK_SC(1,27) ; % Ni (ppm) is not considered in the CREp program.
%data_bulk(1,) = DATA_BULK_SC(1,28) ; % Se is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,) = DATA_BULK_SC(1,29) ; % Br is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,) = DATA_BULK_SC(1,30) ; % Rb is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,32) = DATA_BULK_SC(1,31) ; % Sr (ppm) is not considered in the CREp program.
%data_bulk(1,43) = DATA_BULK_SC(1,32) ; % Zr (ppm) is not considered in the CREp program.
%data_bulk(1,22) = DATA_BULK_SC(1,33) ; % Mo (ppm) is not considered in the CREp program.
%data_bulk(1,) = DATA_BULK_SC(1,34) ; % Ru is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,) = DATA_BULK_SC(1,35) ; % Te is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,) = DATA_BULK_SC(1,36) ; % I is not considered in the Schlagenhauf (2010) program.
%data_bulk(1,9) = DATA_BULK_SC(1,37) ; % Cs (ppm) is not considered in the CREp program.
%data_bulk(1,2) = DATA_BULK_SC(1,38) ; % Ba (ppm) is not considered in the CREp program.
data_bulk(1,62) = data_bulk(1,49)*0.715*10000 ; % Ca (ppm) obtained from the CaO

% TARGET CHEMICAL ELEMENTS

data_target(1,49) = DATA_CREP.CaO_targ(i) ; % CaO (%)
data_target(1,51) = DATA_CREP.K2O_targ(i) ; % K2O (%)
data_target(1,52) = DATA_CREP.TiO2_targ(i) ; % TiO2 (%)
data_target(1,46) = DATA_CREP.Fe2O3_targ(i) ; % Fe2O3 (%)
data_target(1,61) = DATA_CREP.Cl_targ(i) ; % Cl (ppm)
data_target(1,44) = DATA_CREP.SiO2_targ(i) ; % SiO2 (%)
data_target(1,50) = DATA_CREP.Na2O_targ(i) ; % Na2O (%)
data_target(1,48) = DATA_CREP.MgO_targ(i) ; % MgO (%)
data_target(1,45) = DATA_CREP.Al2O3_targ(i) ; % Al2O3 (%)
data_target(1,47) = DATA_CREP.MnO_targ(i) ; % MnO (%)
data_target(1,53) = DATA_CREP.P2O5_targ(i) ; % P2O5 (%)
data_target(1,58) = DATA_CREP.CO2_targ(i) ; % CO2 (%)
data_target(1,56) = DATA_CREP.H2O_targ(i) ; % H20 (%)
%data_target(1,55) = DATA_TARGET_SC(1,14) ; % Li (ppm) is not considered in the CREp program.
%data_target(1,54) = DATA_TARGET_SC(1,15) ; % B (ppm) is not considered in the CREp program.
%data_target(1,30) = DATA_TARGET_SC(1,16) ; % Sm (ppm) is not considered in the CREp program.
%data_target(1,15) = DATA_TARGET_SC(1,17) ; % Gd (ppm) is not considered in the CREp program.
%data_target(1,35) = DATA_TARGET_SC(1,18) ; % Th (ppm) is not considered in the CREp program.
%data_target(1,37) = DATA_TARGET_SC(1,19) ; % U (ppm) is not considered in the CREp program.
%data_target(1,8) = DATA_TARGET_SC(1,20) ; % Cr (ppm) is not considered in the CREp program.
%data_target(1,3) = DATA_TARGET_SC(1,21) ; % Be (ppm) is not considered in the CREp program.
% data_target(1,) = DATA_TARGET_SC(1,22) ; % N is not considered in the Schlagenhauf (2010) program. 
% data_target(1,) = DATA_TARGET_SC(1,23) ; % F is not considered in the Schlagenhauf (2010) program.
% data_target(1,57) = DATA_TARGET_SC(1,24) ; % S (in %) is not considered in the CREp program.
%data_target(1,) = DATA_TARGET_SC(1,25) ; % Sc is not considered in the Schlagenhauf (2010) program.
%data_target(1,7) = DATA_TARGET_SC(1,26) ; % Co (ppm)  is not considered in the CREp program.
%data_target(1,25) = DATA_TARGET_SC(1,27) ; % Ni (ppm)  is not considered in the CREp program.
%data_target(1,) = DATA_TARGET_SC(1,28) ; % Se is not considered in the Schlagenhauf (2010) program.
%data_target(1,) = DATA_TARGET_SC(1,29) ; % Br is not considered in the Schlagenhauf (2010) program.
%data_target(1,) = DATA_TARGET_SC(1,30) ; % Rb is not considered in the Schlagenhauf (2010) program.
%data_target(1,32) = DATA_TARGET_SC(1,31) ; % Sr (ppm)  is not considered in the CREp program.
%data_target(1,43) = DATA_TARGET_SC(1,32) ; % Zr (ppm)  is not considered in the CREp program.
%data_target(1,22) = DATA_TARGET_SC(1,33) ; % Mo (ppm)  is not considered in the CREp program.
%data_target(1,) = DATA_TARGET_SC(1,34) ; % Ru is not considered in the Schlagenhauf (2010) program.
%data_target(1,) = DATA_TARGET_SC(1,35) ; % Te is not considered in the Schlagenhauf (2010) program.
%data_target(1,) = DATA_TARGET_SC(1,36) ; % I is not considered in the Schlagenhauf (2010) program.
%data_target(1,9) = DATA_TARGET_SC(1,37) ; % Cs (ppm) is not considered in the CREp program.
%data_target(1,2) = DATA_TARGET_SC(1,38) ; % Ba (ppm) is not considered in the CREp program.
data_target(1,62) = data_target(1,49)*0.715*10000 ; % Ca (ppm) obtained from the CaO
N36Cl.meas = DATA_CREP.NuclCon(i); % 36Cl conc (at/g)
N36Cl.meas_uncert = DATA_CREP.NuclErr(i);% +/- 36Cl conc (at/g)
N36Cl.inh_0 = .0;
% Uncert_target: Five columns providing the uncertainties on CaO (%), K20 (%),
% TiO2 (%), Fe2O3 (%) and Cl (ppm).
uncert_target(1,1) = DATA_CREP.CaO_targ_er(i); % CaO (%)
uncert_target(1,2) = DATA_CREP.K2O_targ_er(i); % K2O (%)
uncert_target(1,3) = DATA_CREP.TiO2_targ_er(i); % TiO2 (%)
uncert_target(1,4) = DATA_CREP.Fe2O3_targ_er(i); % Fe203 (%)
uncert_target(1,5) = DATA_CREP.Cl_targ_er(i); % Cl (ppm)

%----- Various parameters for each sample

site.sample_name = DATA_CREP.Samples{1,i};
site.lat = DATA_CREP.Lat(i);
site.long = DATA_CREP.Lon(i);
site.alt = DATA_CREP.Alt(i);
site.shield = DATA_CREP.Shield(i);
site.erosion = DATA_CREP.Eros(i); % erosion rate (cm/an)
site.t_form = DATA_CREP.Age_form(i); % formation age (yr) of rock (independently determined or estimated) for radiogenic correction
site.t_form_uncert = DATA_CREP.Age_form_er(i); % formation age (yr) of rock (independently determined or estimated) for radiogenic correction
site.rho_rock = DATA_CREP.Dens(i); % density g/cm2
site.depth = DATA_CREP.Z(i);
site.thick = DATA_CREP.Thick(i);
site.mass_depth = site.depth .* site.rho_rock + (site.thick ./2) .* site.rho_rock;
site.Zs =  DATA_CREP.Z(i) .* site.rho_rock; % mass thickness

site.t_expo_estim = 0; % exposure duration (independently determined or estimated) (yr)

% transfering to structure
Param_site{i} = site;

tmp.bulk = data_bulk;
tmp.target = data_target;
tmp.N36Cl = N36Cl;
tmp.uncert_target = uncert_target;

Data{i} = tmp;
end


