function [P36_tot,P36_s,P36_th,P36_eth,P36_mu] = prodz36_speed(Const_cosmo,Param_cosmo,sf,flag,e)
% The prodz36_speed function computes the 36Cl production rate for a  
% vector of samples as function of the depth. Follows the Gosse & Philipps
% (2001) formulation. Fast muon follow the Heisinger et al. (2002b,a) formulation 
% implemented by Balco et al. (2008).
%
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation.
%           Param_cosmo: Struct. containing the samples specific 
%                            variable previously calculated.                          
%           e : mass-depth (gr/cm2).
%           sf : Scaling factors, previously calculated as function of time
%                   and depth.
%        
% Output : 
%             P36_tot : Total production rate of 36Cl
%             P36_s : Production rate of 36Cl by spallation
%             P36_th : Production rate of 36Cl by thermal neutrons
%             P36_eth : Production rate of 36Cl by epithermal neutrons
%             P36_mu : Production rate of 36Cl by muons
%
% version 01/08/2018, written by TESSON J.

if (max(e) > Param_cosmo.maxdepth)
  disp('maxdepth is');
  Param_cosmo.maxdepth
  disp('z is ');
  e
  error('Prodz called for depth greater than maxdepth.');
end

    expfactor=exp(-e/Const_cosmo.Lambda_f_e);
% Spalation
    P36_s_Ca = sf.currentsf.Sel36Ca .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Ca .* expfactor;
    P36_s_K = sf.currentsf.Sel36K .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_K .* expfactor;
    P36_s_Ti = sf.currentsf.Sel36Ti .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Ti .* expfactor;
    P36_s_Fe = sf.currentsf.Sel36Fe .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Fe .* expfactor;
    P36_s = P36_s_Ca + P36_s_K + P36_s_Ti + P36_s_Fe;

        
if(flag.muon == 'exp') % if muon attenuation at depth are scaled by exponential approximation
   % Muon production rate
    P36_mu = sf.currentsf.SFmu .* Param_cosmo.Y_Sigma.*Const_cosmo.Psi_mu_0.*exp(-e/Const_cosmo.Lambda_mu);
   %depth independent variables for muon-produced neutrons
    P_mu = sf.currentsf.SFmu.*(Param_cosmo.Y_s*Const_cosmo.Psi_mu_0 + 5.8e-6*Const_cosmo.phi_mu_f_0);
   %depth dependent variables for muon-produced neutrons     
   	P_mu_depth = sf.currentsf.SFmu.*(Param_cosmo.Y_s*Const_cosmo.Psi_mu_0 + 5.8e-6*Const_cosmo.phi_mu_f_0) ...
                        .*exp(-e/Const_cosmo.Lambda_mu);
                    
elseif(flag.muon == 'num') % if muons attenuation are computed following Balco 2008.
    % Interpolation muon fluxes (neg , total, and 36Cl muon production rate): all together is faster!!       
        TMP = [Param_cosmo.negflux' Param_cosmo.totalflux' Param_cosmo.Prodmu(1,:)'];
        TMP_interp = interp1(Param_cosmo.depthvector,TMP,e,'linear')';
        negfluxdepth = TMP_interp(1,:);
        totalfluxdepth =  TMP_interp(2,:);
        P36_mu = sf.S_T .* TMP_interp(3,:);
    %depth dependent variables for muon-produced neutrons
    P_mu_depth=Param_cosmo.Y_s.*negfluxdepth+0.0000058.*totalfluxdepth;
    %depth independent variables for muon-produced neutrons
    P_mu=Param_cosmo.Y_s.*Param_cosmo.negflux(1)+0.0000058.*Param_cosmo.totalflux(1);  
end    
    
    %depth dependent variables for muon-produced neutrons   
    R_mu_depth=P_mu_depth./(sf.currentsf.SFeth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);
    R_prime_mu_depth=R_mu_depth.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);

    %depth independent variables for muon-produced neutrons
    R_mu=P_mu./(sf.currentsf.SFth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);
    R_prime_mu=R_mu.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);

% Thermal neutron production 
    phi_th = sf.currentsf.SFth .* sf.S_T .* (Param_cosmo.phi_star_th.*exp(-e./Const_cosmo.Lambda_f_e) + ...
    (1 + R_prime_mu).*Param_cosmo.JDeltaphi_star_eth.*exp(-e./Param_cosmo.L_eth) + ...
    (1 + R_prime_mu.*Param_cosmo.R_th).*Param_cosmo.JDeltaphi_star_th.*exp(-e./Param_cosmo.L_th) + ...
    R_prime_mu_depth.*Param_cosmo.phi_star_th) ; % Thermal neutron flux (n.cm_2.a-1)

% Get thermal production using the flux
    P36_th =(Param_cosmo.f_th./Param_cosmo.Lambda_th).*phi_th;

% Epithermal neutron production
    phi_eth = sf.currentsf.SFeth .* sf.S_T .* (Param_cosmo.phi_star_eth .* exp(-e./Const_cosmo.Lambda_f_e) + ...
     (1 + R_mu .* Param_cosmo.R_eth) .* Param_cosmo.FDeltaphi_star_eth .* exp(-e./Param_cosmo.L_eth) + ...
     R_mu_depth .* Param_cosmo.phi_star_eth) ; % Epithermal neutron flux (concentration) (n cm-2 yr-1)
% Get epithermal production using the flux 
    P36_eth = (Param_cosmo.f_eth ./ Param_cosmo.Lambda_eth) .* phi_eth .* (1 - Param_cosmo.p_E_th);

% Radiogenic production rate
    P36_rad = Param_cosmo.P_rad ;
    
% Total production rate
    P36_tot = P36_s_Ca + P36_s_K + P36_s_Ti + P36_s_Fe + P36_eth + P36_th + P36_mu + P36_rad;
    
end

