function [ N_calc_tot,N_calc_tot_uncert] = depth_profile_speed(...
                   t_vector,z_history,Const_cosmo,...
                   Param_cosmo_in,Param_cosmo_ICE_in,...
                   Param_site_in,Param_site_ICE_in,...
                   sf,sf_ICE,erosion,t_expo,z_vector,flag)
               
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

    % Time and position increment
    deltat = 100; % (yr)

    % time vector
    %t_vector = [0:100:t_expo];
    
    % Theoretical 36Cl vector
    N_calc_tot=zeros(1,N_samples);
    N_calc_tot_uncert=zeros(1,N_samples);
    
    % get current scaling factors.
    currentsf=getcurrentsf(sf,t_vector,flag); 
    sf.currentsf = currentsf;
    
    % get current scaling factors for ice covering
    currentsf=getcurrentsf(sf,t_vector,flag); 
    sf_ICE.currentsf = currentsf;
    
    for i=1:length(z_vector)
            
            % depth-time vector
            depth_time_vector = z_history + z_vector(i); %t_vector .* erosion + z_vector(i); % evolution of the sample depth (cm) over the time
            % depth-time vector for surface sample
            depth_time_vector_surf = t_vector .* 0.0;
            
            % Rock parameters
            Param_cosmo = Param_cosmo_in{i};
            Param_site = Param_site_in{i};
            
            % Ice parameters
            Param_cosmo_ICE = Param_cosmo_ICE_in{i};
            Param_site_ICE = Param_site_ICE_in{i};
            
            % get Production rate within the sample at surface z=.0;
            
            
                % get production rate over the whole depth time vector, and
                % averaging over the whole sample.
         
                n_thick = 10;
                d_thick = Param_site.thick/n_thick.*[0:1:n_thick];      
          
                d_integ_samp = z_vector(i)-Param_site.thick/2+d_thick; % z_vector = center of the sample
                d_integ_samp=d_integ_samp(d_integ_samp >= 0); % exclude negative value in the case z of the sample is zero, with a thickness > 0.

                P_cosmo = zeros(1,length(depth_time_vector));
                P36_s = zeros(1,length(depth_time_vector));
                P36_th = zeros(1,length(depth_time_vector));
                P36_eth = zeros(1,length(depth_time_vector));
                P36_mu = zeros(1,length(depth_time_vector));
                P36_rad = zeros(1,length(depth_time_vector));
                
                for j=1:length(d_integ_samp) % loop over the thickness
                    %depth_time_vector_thick = t_vector .* erosion + d_integ_samp(j);
                    depth_time_vector_thick = t_vector .* erosion .* 0.0 + d_integ_samp(j);
                    % get production rates
                    [~,P36_s_tmp,P36_th_tmp,P36_eth_tmp,P36_mu_tmp,P36_rad_tmp] = prodz36_speed(Const_cosmo,Param_cosmo,sf,flag,depth_time_vector_thick*Param_site.rho_rock);
                    % total production rate within the sample at surface
                    P_cosmo_tmp = P36_s_tmp+P36_th_tmp+P36_eth_tmp+P36_mu_tmp;
                    
                    % Sum production rate over the sample thickness
                    P_cosmo = P_cosmo + P_cosmo_tmp;
                    
                    % details for each pathway
                    P36_s = P36_s + P36_s_tmp;
                    P36_th = P36_th + P36_th_tmp;
                    P36_eth = P36_eth + P36_eth_tmp;
                    P36_mu = P36_mu + P36_mu_tmp;
                    P36_rad = P36_rad + P36_rad_tmp;
                    
                end
                % Average production rate over the sample thickness
                P_cosmo = P_cosmo./length(d_integ_samp); 
                P36_s = P36_s./length(d_integ_samp);
                P36_th = P36_th./length(d_integ_samp);
                P36_eth = P36_eth./length(d_integ_samp);
                P36_rad = P36_rad./length(d_integ_samp);
                P36_mu = P36_mu./length(d_integ_samp);

            % Compute scaling due to the ice cover
                % Production rate within ice at surface
                [~,P36_s_ice_0,P36_th_ice_0,P36_eth_ice_0,P36_mu_ice_0,P36_rad_ice_0] = prodz36_speed(Const_cosmo,Param_cosmo_ICE,sf_ICE,flag,depth_time_vector_surf);
                P_cosmo_ice_zero = P36_s_ice_0 + P36_th_ice_0 + P36_eth_ice_0 + P36_mu_ice_0;
                
                %  Production rate within ice at depth
                [~,P36_s_ice,P36_th_ice,P36_eth_ice,P36_mu_ice,P36_rad_ice] = prodz36_speed(Const_cosmo,Param_cosmo_ICE,sf_ICE,flag,depth_time_vector*Param_site_ICE.rho_rock);
                P_cosmo_ice_depth = P36_s_ice + P36_th_ice + P36_eth_ice + P36_mu_ice;
                
                % Get scaling factor accounting for the ice thickness
                Scaling_ice = P_cosmo_ice_depth./P_cosmo_ice_zero;
                
            % Get the scaled Production rate
            %P_tot_cor = P_cosmo .* Scaling_ice + P36_rad;
            P36_s_cor = P36_s .* P36_s_ice./P36_s_ice_0;
            P36_eth_cor = P36_eth .* P36_eth_ice./P36_eth_ice_0;
            P36_th_cor = P36_th .* P36_th_ice./P36_th_ice_0;
            P36_mu_cor = P36_mu .* P36_mu_ice./P36_mu_ice_0;
            
            P_tot_cor = P36_s_cor+P36_eth_cor+P36_th_cor+P36_th_cor+ P36_rad;
            
        % Get the amount of 36Cl for each time step
          N36_prod = P_tot_cor.*(1.0-exp(-Const_cosmo.lambda36*deltat))./Const_cosmo.lambda36;
          N36_cum = flip(cumsum(flip(N36_prod)));
        % radiogenic decay of the 36Cl stock
          N36_rad = N36_cum .* (1-exp(-Const_cosmo.lambda36*deltat));
        % Total 36Cl 
          N_calc_tot(i) = sum(N36_prod) - sum(N36_rad);
          N_calc_tot_uncert(i) = .0;
          
        if(flag.plot)
            figure(1);  
      
            subplot(2,2,1); 
            %plot(t_vector,P_cosmo); hold on;
            plot(-t_vector,P_tot_cor); hold on;
%              plot(-t_vector,P36_s_cor); 
%              plot(-t_vector,P36_eth_cor); 
%              plot(-t_vector,P36_th_cor); 
%              plot(-t_vector,P36_mu_cor); 
%             
            
            subplot(2,2,2); 
            plot(-t_vector,depth_time_vector./100); hold on;
            
            subplot(2,2,3);  
            plot(-t_vector,(N36_cum-flip(cumsum(flip(N36_rad)))));hold on;
            
            subplot(2,2,4);
            plot(-t_vector,(N36_cum-flip(cumsum(flip(N36_rad))))./N_calc_tot(i).*100);hold on;
            
            figure(2);
            subplot(2,1,1)
            plot(-t_vector,P_cosmo_ice_depth); hold on;
            plot(-t_vector,P36_s_ice); 
            plot(-t_vector,P36_th_ice); 
            plot(-t_vector,P36_eth_ice); 
            plot(-t_vector,P36_mu_ice); 
            
            subplot(2,1,2)
            plot(-t_vector,P_cosmo_ice_depth./P_cosmo_ice_depth); hold on;
            plot(-t_vector,P36_s_ice./P_cosmo_ice_depth);
            plot(-t_vector,P36_th_ice./P_cosmo_ice_depth); 
            plot(-t_vector,P36_eth_ice./P_cosmo_ice_depth); 
            plot(-t_vector,P36_mu_ice./P_cosmo_ice_depth); 
        end
    end
%end       


