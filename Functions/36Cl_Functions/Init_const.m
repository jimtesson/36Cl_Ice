function [Const_cosmo] = Init_const(DATA_CREP,flag)
% This function initialize the cosmogenic constants

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


end


