function [ N_calc_tot,N_calc_tot_uncert] = depth_profile_speed(Const_cosmo, Param_cosmo_in,Param_site_in,sf,erosion,t_expo,z_vector,flag)
% This function computes the theoretical 36Cl for a given vector of samples
% , an fixed erosion rate and exposure duration. User can choose between 2 
%   models: 1) the Schimmelfenning 2009 model, 2) numerical approach using 
%   the LSD model.
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation.
%           Param_cosmo_in : Struct. containing the samples specific 
%                            variable previously calculated.                          
%           Param_site_in : Site specific parameters.
%           sf : Scaling factors, previously calculated as function of time
%                   and depth.
%           z_vector : depth vector of the samples
%           t_expo : exposure age of the samples
%           erosion : user-provided erosion-rate in m/Myr
%           flag : user-provided flag to specify the model used.
%
%                   flag.model = 'exp' for the Schimmelpfennig 2009 model
%                   (attenuation length approximated by exponential)
%
%                   flag.model = 'num' for the numerical approach (Marrero
%                   et al. 2018, attenuation length are calculated from the
%                   energy flux))
%
%                   flag.scaling_model = 'st' for the Stone 1999 scaling 
%                                               scheme
%
%                   flag.scaling_model = 'sa' for the Lifton-Sato scaling 
%                                               scheme
%
%
%           
% Output : 
%             N_calc_tot: 36Cl concentration for each sample
%             N_calc_tot_uncert: uncertainties on the 36Cl.
% version 01/08/2018, written by TESSON J.

erosion = erosion*1E-4; % mm/kyr to cm/yr
N_samples = length(z_vector);


%% time factors
% if(strcmp(flag.model,'exp'))
%     for i=1:length(z_vector)
%         z = z_vector(i);
%         Param_cosmo = Param_cosmo_in{i};
%         Param_site = Param_site_in{i};
%         % time factors
%         Param_cosmo.tcosm_sp = (1-exp(-(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e)*t_expo))...
%                         /(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e);
%         Param_cosmo.tcosm_sp_uncert = ((t_expo_uncert*exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e)))^2 ...
%                                 +(erosion_uncert*((-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e)) ... 
%                                 *(-t_expo*Param_site.rho_rock/Const_cosmo.Lambda_f_e)* ... 
%                                 (Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e))...
%                                 -(1-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e)))...
%                                 *Param_site.rho_rock/Const_cosmo.Lambda_f_e)/(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_f_e)^2)^2)^.5;
% 
%         Param_cosmo.tcosm_th = (1-exp(-(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th)*t_expo))...
%                     /(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th);
%         Param_cosmo.tcosm_th_uncert = ((t_expo_uncert*exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th)))^2+(erosion_uncert*((-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th))*(-t_expo*Param_site.rho_rock/Param_cosmo.L_th)*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th))-(1-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th)))*Param_site.rho_rock/Param_cosmo.L_th)/(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_th)^2)^2)^.5;
%                 
%         Param_cosmo.tcosm_eth = (1-exp(-(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth)*t_expo))...
%                     /(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth);
%         Param_cosmo.tcosm_eth_uncert = ((t_expo_uncert*exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth)))^2+(erosion_uncert*((-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth))*(-t_expo*Param_site.rho_rock/Param_cosmo.L_eth)*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth))-(1-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth)))*Param_site.rho_rock/Param_cosmo.L_eth)/(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Param_cosmo.L_eth)^2)^2)^0.5;
% 
%         Param_cosmo.tcosm_mu = (1-exp(-(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu)*t_expo))...
%                     /(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu);
%         Param_cosmo.tcosm_mu_uncert = ((t_expo_uncert*exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu)))^2+(erosion_uncert*((-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu))*(-t_expo*Param_site.rho_rock/Const_cosmo.Lambda_mu)*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu))-(1-exp(-t_expo*(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu)))*Param_site.rho_rock/Const_cosmo.Lambda_mu)/(Const_cosmo.lambda36+Param_site.rho_rock*erosion/Const_cosmo.Lambda_mu)^2)^2)^.5;
% 
%         Param_cosmo.tcosm = (1-exp(-Const_cosmo.lambda36*t_expo))/Const_cosmo.lambda36;
%         t_expo_er = 0.1.*t_expo; % uncertainties on t_expo
%         Param_cosmo.tcosm_uncert = ((Const_cosmo.lambda36_uncert*((-t_expo*(-exp(-Const_cosmo.lambda36*t_expo_er))*Const_cosmo.lambda36)-(1-exp(-Const_cosmo.lambda36*t_expo)))/Const_cosmo.lambda36^2)^2+(t_expo_er*(exp(-Const_cosmo.lambda36*t_expo)))^2)^.5;
%         % production rate coefficient corrected for sample thickness
%         Q_sp = Param_cosmo.Q_sp;
%         Q_th = Param_cosmo.Q_th;
%         Q_eth = Param_cosmo.Q_eth;
%         Q_mu = Param_cosmo.Q_mu;
%         S_shape = sf.S_shape;
%         S_snow = sf.S_snow;
%         S_T = sf.S_T;
%         C_Ca = Param_cosmo.C_Ca;
%         C_K = Param_cosmo.C_K;
%         C_Ti = Param_cosmo.C_Ti;
%         C_Fe = Param_cosmo.C_Fe;
%         Psi_Cl36_Ca_0 = Const_cosmo.Psi_Cl36_Ca_0 ;% spallation production rate for Ca, SLHL
%         Psi_Cl36_K_0 = Const_cosmo.Psi_Cl36_K_0 ;% Spallation production rate at surface of 39K
%         Psi_Cl36_Ti_0 = Const_cosmo.Psi_Cl36_Ti_0; % Spallation production rate at surface of Ti
%         Psi_Cl36_Fe_0 = Const_cosmo.Psi_Cl36_Fe_0; % Spallation production rate at surface of Fe
% 
%         p_E_th = Param_cosmo.p_E_th;
%         f_eth = Param_cosmo.f_eth;
%         Lambda_eth = Param_cosmo.Lambda_eth;
%         phi_star_eth = Param_cosmo.phi_star_eth;
%         f_th = Param_cosmo.f_th;
%         Lambda_th = Param_cosmo.Lambda_th;
%         phi_star_th = Param_cosmo.phi_star_th;
%         R_mu = Param_cosmo.R_mu;
%         R_eth = Param_cosmo.R_eth;
%         FDeltaphi_star_eth = Param_cosmo.FDeltaphi_star_eth;
%         JDeltaphi_star_eth = Param_cosmo.JDeltaphi_star_eth;
%         R_prime_mu = Param_cosmo.R_prime_mu;
%         R_th = Param_cosmo.R_th;
%         JDeltaphi_star_th = Param_cosmo.JDeltaphi_star_th;
%         S_el_mu = sf.SF_St_mu;
%         S_el_sp_n = sf.SF_St_sp;
%         Y_Sigma = Param_cosmo.Y_Sigma;
%         Psi_mu_0 = Const_cosmo.Psi_mu_0;
% 
%         % Production rate
%              %P_tot = S_el_sp_n .* S_T .* (Param_cosmo.Q_sp .* S_shape .* S_snow .* (Param_cosmo.P_sp_Fe+Param_cosmo.P_sp_Ti+Param_cosmo.P_sp_K+Param_cosmo.P_sp_Ca)...
% %                                         + Param_cosmo.Q_th .* Param_cosmo.P_th + Param_cosmo.Q_eth .* Param_cosmo.P_eth)...
% %                                         + S_el_mu .* S_T .* Param_cosmo.Q_mu .* Param_cosmo.P_mu;
%             % Theoretical 36Cl                 
%             N_calc_tot(i) = S_el_sp_n .* S_T .*( Param_cosmo.tcosm_sp .* Param_cosmo.Df_s .* Param_cosmo.JQ_s ...
%                             + Param_cosmo.tcosm_th .* Param_cosmo.Df_th .* Param_cosmo.JQ_th ...
%                             + Param_cosmo.tcosm_eth .* Param_cosmo.Df_eth .* Param_cosmo.JQ_eth ...
%                             + Param_cosmo.tcosm_mu .* Param_cosmo.Df_mu .* Param_cosmo.JQ_mu) ...
%                             +Param_cosmo.N36Cl.rad + Param_cosmo.N36Cl.inh;
%             N_calc_tot_uncert(i) = ((1*Param_cosmo.N36Cl.rad_uncert).^2+(Param_cosmo.P_uncert*Param_cosmo.tcosm).^2+(Param_cosmo.tcosm_uncert*Param_cosmo.P_cosmo).^2).^0.5;
% 
% %             N_s(i) = S_el_sp_n*Param_cosmo.tcosm_sp * Param_cosmo.J_s *exp(-z*Param_site.rho_rock/Const_cosmo.Lambda_f_e);
% %             N_eth(i) = S_el_sp_n*Param_cosmo.tcosm_eth * Param_cosmo.J_eth *exp(-z*Param_site.rho_rock/Param_cosmo.Lambda_eth);
% %             N_th(i) = S_el_sp_n*Param_cosmo.tcosm_th * Param_cosmo.J_th *exp(-z*Param_site.rho_rock/Param_cosmo.Lambda_th);
% %             N_mu(i) = S_el_sp_n*Param_cosmo.tcosm_mu * Param_cosmo.J_mu *exp(-z*Param_site.rho_rock/Const_cosmo.Lambda_mu);
% %             %N_mu2(i)= N_mu_Balco(Const_cosmo, Param_cosmo,z.*Param_site.rho_rock,erosion.*Param_site.rho_rock,t_expo);
% %             N_calc_tot(i) = S_T*(N_s(i)+N_eth(i)+N_th(i)+N_mu(i))+Param_cosmo.N36Cl.rad+Param_cosmo.N36Cl.inh;
%         % Uncertainties
% %             N_s_uncert(i) = ((Param_cosmo.tcosm_sp_uncert * Param_cosmo.J_s *exp(-z*Param_site.rho_rock/Const_cosmo.Lambda_f_e))^2)^.5;
% %             N_eth_uncert(i) = ((Param_cosmo.tcosm_eth_uncert * Param_cosmo.J_eth *exp(-z*Param_site.rho_rock/Param_cosmo.Lambda_eth))^2)^.5;
% %             N_th_uncert(i) = ((Param_cosmo.tcosm_th_uncert * Param_cosmo.J_th *exp(-z*Param_site.rho_rock/Param_cosmo.Lambda_th))^2)^.5;
% %             N_mu_uncert(i) = ((Param_cosmo.tcosm_mu_uncert * Param_cosmo.J_mu *exp(-z*Param_site.rho_rock/Const_cosmo.Lambda_mu))^2)^.5;
% %           N_calc_tot_uncert(i) = (S_el_sp_n^2*S_T^2*...
% %                    (N_s_uncert(i)^2+N_eth_uncert(i)^2+N_th_uncert(i)^2+N_mu_uncert(i)^2)...
% %                    +Param_cosmo.N36Cl.rad_uncert^2)^.5;
%     end
%     
%elseif(strcmp(flag.model,'num')) % sato atmospheric and depth muonic flux calculation
    
    % Time and position increment
    deltat = 100; % (yr)

    % time vector
    t_vector = [0:100:t_expo];
    
    % Theoretical 36Cl vector
    N_calc_tot=zeros(1,N_samples);
    N_calc_tot_uncert=zeros(1,N_samples);
    
    % get current scaling factors.
    currentsf=getcurrentsf(sf,t_vector,flag); 
    sf.currentsf = currentsf;
    
    for i=1:length(z_vector)
        
            depth_time_vector = t_vector .* erosion + z_vector(i);
            
            Param_cosmo = Param_cosmo_in{i};
            Param_site = Param_site_in{i};
            
            
         % get production rate over the whole depth time vector, and
         % averaging over the whole sample.
         
            n_thick = 10;
            d_thick = Param_site.thick/n_thick.*[0:1:n_thick];      
            P_tot = zeros(1,length(depth_time_vector));
            d_integ_samp = z_vector(i)-Param_site.thick/2+d_thick; % z_vector = center of the sample
            d_integ_samp=d_integ_samp(d_integ_samp >= 0); % exclude negative value in the case z of the sample is zero, with a thickness > 0.

         for j=1:length(d_integ_samp) % loop over the thickness
             depth_time_vector_thick = t_vector .* erosion + d_integ_samp(j);
             P_tot_tmp = prodz36_speed(Const_cosmo,Param_cosmo,sf,flag,depth_time_vector_thick*Param_site.rho_rock);
             P_tot = P_tot + P_tot_tmp;
         end
            P_tot=P_tot./length(d_integ_samp); % average production rate
           
          %P_tot = prodz36_speed(Const_cosmo,Param_cosmo,sf,depth_time_vector*Param_site.rho_rock);
           
        % Get the amount of 36Cl for each time step
          N36_prod = P_tot.*(1.0-exp(-Const_cosmo.lambda36*deltat))./Const_cosmo.lambda36;
          N36_cum = flip(cumsum(flip(N36_prod)));
        % radiogenic decay of the 36Cl stock
          N36_rad = N36_cum .* (1-exp(-Const_cosmo.lambda36*deltat));
        % Total 36Cl 
          N_calc_tot(i) = sum(N36_prod) - sum(N36_rad);
          N_calc_tot_uncert(i) = .0;
    end
%end       


