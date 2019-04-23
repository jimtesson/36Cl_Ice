function [P36_tot_uncert] = get_36Cl_uncert(Data_in,Const_cosmo,Param_cosmo,Param_site_in,sf,e,flag)



    
    %% Uncertainty on the spallation production rate
        expfactor=exp(-e/Const_cosmo.Lambda_f_e);
    % Spalation
        P36_s_Ca = sf.currentsf.Sel36Ca .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Ca .* expfactor;
        P36_s_K = sf.currentsf.Sel36K .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_K .* expfactor;
        P36_s_Ti = sf.currentsf.Sel36Ti .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Ti .* expfactor;
        P36_s_Fe = sf.currentsf.Sel36Fe .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.P_sp_Fe .* expfactor;
        P36_s = P36_s_Ca + P36_s_K + P36_s_Ti + P36_s_Fe;
    
        % uncertainties
        P36_s_Ca_uncert = ((sf.currentsf.Sel36Ca .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.uncert.P_sp_Ca_uncert .* expfactor ).^2).^.5; 
        P36_s_K_uncert =  ((sf.currentsf.Sel36K .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.uncert.P_sp_K_uncert .* expfactor).^2).^.5;
        P36_s_Ti_uncert = ((sf.currentsf.Sel36Ti .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.uncert.P_sp_Ti_uncert .* expfactor).^2).^.5; 
        P36_s_Fe_uncert = ((sf.currentsf.Sel36Fe .*sf.S_T .* sf.S_snow .* sf.S_shape .* Param_cosmo.uncert.P_sp_Fe_uncert .* expfactor).^2).^.5; 
        P36_s_uncert =((P36_s_Ca_uncert).^2+(P36_s_K_uncert).^2+(P36_s_Ti_uncert).^2+(P36_s_Fe_uncert).^2).^.5;
        
        % P36_s_uncert = P36_s_uncert./P36_s; % get uncertainty in percent 
        
    % Radiogenic production rate
        P36_rad = Param_cosmo.P_rad;

    
    %% Get Muon fluxes
    if(flag.muon == 'exp') % if muon attenuation at depth are scaled by exponential approximation
    % Muon production rate
        P36_mu = sf.currentsf.SFmu .* Param_cosmo.Y_Sigma.*Const_cosmo.Psi_mu_0*exp(-e/Const_cosmo.Lambda_mu);
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
    
    
    %% Uncertainty on the Thermal and Epithermal neutrons from P_f_0
    
    if(strcmp(flag.uncert.pf0_mc,'yes')) % Monte Carlo sampling of P_f_0; yes or no ?
        n_test = 500; % number of test
        % initialization results vectors
        xx=zeros(n_test);
        yy=xx;
        zz=xx;
     for i=1:n_test % Monte carlo sampling for uncertainty on P_th and P_eth
        
        % make vary P_f_0 within its uncertainty
        P_f_0=.0;
        while(P_f_0<=0) %P_f_0 cannot be negative!
            P_f_0=Param_cosmo.P_f_0 + Param_cosmo.P_f_0_uncert*randn(1);
        end
        % Epi-Thermal
        phi_star_eth = P_f_0*Param_cosmo.R_eth/(Param_cosmo.Sigma_eth - (Param_cosmo.D_eth/(Const_cosmo.Lambda_f_e^2))) ; % Epithermal neutron flux at land/atmosphere
        
        phi_star_eth_a = P_f_0*Param_cosmo.R_eth_a/(Param_cosmo.Sigma_eth_a - (Param_cosmo.D_eth_a/(Const_cosmo.Lambda_f_e^2))) ; % Epithermal neutron flux at land/atmosphere interface that would be observed in atm 

        R_mu=P_mu./(sf.currentsf.SFth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);
        
        R_mu_depth=P_mu_depth./(sf.currentsf.SFeth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);        

        Deltaphi_2star_eth_a = phi_star_eth - Param_cosmo.D_eth_a*phi_star_eth_a/Param_cosmo.D_eth ; % Adjusted difference between hypothetical equilibrium epithermal neutron fluxes in atm and ss (n cm-2 yr-1)

        FDeltaphi_star_eth = ((Param_cosmo.D_eth_a/Param_cosmo.L_eth_a)*(phi_star_eth_a - phi_star_eth) - ...
                                Deltaphi_2star_eth_a*(Param_cosmo.D_eth/Const_cosmo.Lambda_f_e))/...
                                ((Param_cosmo.D_eth_a/Param_cosmo.L_eth_a) + (Param_cosmo.D_eth/Param_cosmo.L_eth)) ; % EQ. 3.28 Gosse & Phillips, 2001

        phi_eth = sf.currentsf.SFeth .* sf.S_T .* (phi_star_eth .* exp(-e./Const_cosmo.Lambda_f_e) + ...
                    (1 + R_mu .* Param_cosmo.R_eth) .* FDeltaphi_star_eth .* exp(-e./Param_cosmo.L_eth) + ...
                    R_mu_depth .* phi_star_eth) ; % Epithermal neutron flux (concentration) (n cm-2 yr-1)
        
        P_eth = (Param_cosmo.f_eth ./ Param_cosmo.Lambda_eth) .* phi_eth .* (1 - Param_cosmo.p_E_th);
                    
        % Thermal neutrons
    
        Deltaphi_star_eth_a = phi_star_eth - phi_star_eth_a ; % difference in equilibrium epithermal neutron fluxes between atm and ss

        FDeltaphi_star_eth_a = (Param_cosmo.D_eth*Deltaphi_star_eth_a/Param_cosmo.L_eth - Param_cosmo.D_eth*Deltaphi_2star_eth_a/Const_cosmo.Lambda_f_e)/ ...
                                (Param_cosmo.D_eth_a / Param_cosmo.L_eth_a + Param_cosmo.D_eth / Param_cosmo.L_eth );
                            
        phi_star_th = (Param_cosmo.p_E_th_a*Param_cosmo.R_th*phi_star_eth)/...
                      (Param_cosmo.Lambda_eth*(Param_cosmo.Sigma_th - Param_cosmo.DD_th/(Const_cosmo.Lambda_f_e^2))) ;
%        phi_star_th = (p_E_th_a*R_th*phi_star_eth)/(Lambda_eth*(Sigma_th - DD_th/(Lambda_e^2))) ;

        R_prime_mu=R_mu.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);
        
        R_prime_mu_depth = R_mu_depth.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);
    
        JDeltaphi_star_eth = (Param_cosmo.p_E_th_a*Param_cosmo.R_th*FDeltaphi_star_eth)...
                            /(Param_cosmo.Lambda_eth*(Param_cosmo.Sigma_th ...
                            - Param_cosmo.DD_th/(Param_cosmo.L_eth^2))) ; % Eq. 3.39 Gosse & Phillips, 2001

        JDeltaphi_star_eth_a = (Param_cosmo.p_E_th_a*Param_cosmo.R_th_a*FDeltaphi_star_eth_a)/((1/Param_cosmo.Sigma_eth_a)*(Param_cosmo.Sigma_th_a - Param_cosmo.D_th_a/(Param_cosmo.L_eth_a^2))) ;
        
        phi_star_th_a = (Param_cosmo.p_E_th_a*Param_cosmo.R_th_a*phi_star_eth_a)/(1/Param_cosmo.Sigma_eth_a*(Param_cosmo.Sigma_th_a - Param_cosmo.D_th_a/(Const_cosmo.Lambda_f_e^2))) ; 

        Deltaphi_star_th = phi_star_th_a - phi_star_th ; % difference between hypothetical equilibrium thermal neutron fluxes in atmosphere and ss

        JDeltaphi_star_th = (Param_cosmo.D_th_a*(phi_star_th_a/Const_cosmo.Lambda_f_e - JDeltaphi_star_eth_a/Param_cosmo.L_eth_a) ...
                             - Param_cosmo.DD_th*(phi_star_th/Const_cosmo.Lambda_f_e + JDeltaphi_star_eth/Param_cosmo.L_eth) ...
                                + (Param_cosmo.D_th_a/Param_cosmo.L_th_a)*(Deltaphi_star_th + JDeltaphi_star_eth_a - JDeltaphi_star_eth)) ...
                            / ((Param_cosmo.DD_th/Param_cosmo.L_th) + (Param_cosmo.D_th_a/Param_cosmo.L_th_a)) ; % portion of difference between phi_star_th,ss and actual flux due to thermal flux profile

        % Thermal neutron production 
        phi_th = sf.currentsf.SFth .* sf.S_T .* (phi_star_th.*exp(-e./Const_cosmo.Lambda_f_e) + ...
                (1 + R_prime_mu).*JDeltaphi_star_eth.*exp(-e./Param_cosmo.L_eth) + ...
                (1 + R_prime_mu.*Param_cosmo.R_th).*JDeltaphi_star_th.*exp(-e./Param_cosmo.L_th) + ...
                R_prime_mu_depth.*phi_star_th) ; % Thermal neutron flux (n.cm_2.a-1)

        % Get thermal production using the flux
        P_th =(Param_cosmo.f_th./Param_cosmo.Lambda_th).*phi_th;                
        xx(i) = P_eth(1);
        yy(i) = P_th(1);
    end
     
             %P_eth_uncert = std(xx)./mean(xx); % in percent
             %P_th_uncert = std(yy)./mean(yy); % in percent
             P_eth_uncert = std(xx); % in percent
             P_th_uncert = std(yy); % in percent

    else % if not using the MC test.
        
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

        P_eth_uncert = P36_eth.*0.25; % uncertainty is 25%, found by Monte-Carlo test
        P_th_uncert = P36_th.*0.25;  % uncertainty is 25%, found by Monte-Carlo test
    end        
    
 %% Propagate uncertainty from the 35Cl concentration
 
     % Monte carlo for uncertainty on P_th and P_eth
     n_test = 500; % number of test
     
    for i=1:n_test
        
        Avogadro = 6.022E+23 ; % Avogadro Number
        
        % generate the random 35Cl concentration using its uncertainty
        tmp_Cl=.0;
        while(tmp_Cl<=0)
            tmp_Cl = Data_in{1}.target(1,61) + randn(1).*Data_in{1}.uncert_target(1,5);
            Param_cosmo.ppm_targ(61) = tmp_Cl;
        end
    
        % Compute all variables required to compute P_th and P_eth
        
        N_Cl_targ = (Param_cosmo.ppm_targ(61)./Param_cosmo.A_k(61))*Avogadro*1e-6;  % Concentrations in atom/g
        N_k_targ = (Param_cosmo.ppm_targ(:,1:61)./Param_cosmo.A_k)*Avogadro*1e-6 ; % Concentrations in atom/g
        N_k_targ(56) = N_k_targ(56)/Param_site_in.rho_rock ; % divided by bulk-rock density according to CHLOE for H

        % Fraction of epith neutrons absorbed by Cl35
        f_eth = N_Cl_targ*Param_cosmo.I_a_k(61)*(1e-24)/Param_cosmo.I_eff;  % (Eq 3.17, Gosse & Phillips, 2001)
        f_th = Param_cosmo.sigma_th_k(61)*N_Cl_targ*1e-24/Param_cosmo.Sigma_th ; % Eq 3.32 de Gosse and Phillips, 2001

        P_f_0=Param_cosmo.P_f_0;
    
        % Epi-Thermal
        phi_star_eth = P_f_0*Param_cosmo.R_eth/(Param_cosmo.Sigma_eth - (Param_cosmo.D_eth/(Const_cosmo.Lambda_f_e^2))) ; % Epithermal neutron flux at land/atmosphere
        
        phi_star_eth_a = P_f_0*Param_cosmo.R_eth_a/(Param_cosmo.Sigma_eth_a - (Param_cosmo.D_eth_a/(Const_cosmo.Lambda_f_e^2))) ; % Epithermal neutron flux at land/atmosphere interface that would be observed in atm 

        R_mu=P_mu./(sf.currentsf.SFth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);
        
        R_mu_depth=P_mu_depth./(sf.currentsf.SFeth.*Param_cosmo.P_f_0.*Param_cosmo.R_eth);        

        Deltaphi_2star_eth_a = phi_star_eth - Param_cosmo.D_eth_a*phi_star_eth_a/Param_cosmo.D_eth ; % Adjusted difference between hypothetical equilibrium epithermal neutron fluxes in atm and ss (n cm-2 yr-1)

        FDeltaphi_star_eth = ((Param_cosmo.D_eth_a/Param_cosmo.L_eth_a)*(phi_star_eth_a - phi_star_eth) - ...
                                Deltaphi_2star_eth_a*(Param_cosmo.D_eth/Const_cosmo.Lambda_f_e))/...
                                ((Param_cosmo.D_eth_a/Param_cosmo.L_eth_a) + (Param_cosmo.D_eth/Param_cosmo.L_eth)) ; % EQ. 3.28 Gosse & Phillips, 2001

        phi_eth = sf.currentsf.SFeth .* sf.S_T .* (phi_star_eth .* exp(-e./Const_cosmo.Lambda_f_e) + ...
                    (1 + R_mu .* Param_cosmo.R_eth) .* FDeltaphi_star_eth .* exp(-e./Param_cosmo.L_eth) + ...
                    R_mu_depth .* phi_star_eth) ; % Epithermal neutron flux (concentration) (n cm-2 yr-1)
        
                    
        % Thermal neutrons
        Deltaphi_star_eth_a = phi_star_eth - phi_star_eth_a ; % difference in equilibrium epithermal neutron fluxes between atm and ss

        FDeltaphi_star_eth_a = (Param_cosmo.D_eth*Deltaphi_star_eth_a/Param_cosmo.L_eth - Param_cosmo.D_eth*Deltaphi_2star_eth_a/Const_cosmo.Lambda_f_e)/ ...
                                (Param_cosmo.D_eth_a / Param_cosmo.L_eth_a + Param_cosmo.D_eth / Param_cosmo.L_eth );
                            
        phi_star_th = (Param_cosmo.p_E_th_a*Param_cosmo.R_th*phi_star_eth)/...
                      (Param_cosmo.Lambda_eth*(Param_cosmo.Sigma_th - Param_cosmo.DD_th/(Const_cosmo.Lambda_f_e^2))) ;

        R_prime_mu=R_mu.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);
        
        R_prime_mu_depth = R_mu_depth.*(Param_cosmo.p_E_th_a./Param_cosmo.p_E_th);
    
        JDeltaphi_star_eth = (Param_cosmo.p_E_th_a*Param_cosmo.R_th*FDeltaphi_star_eth)...
                            /(Param_cosmo.Lambda_eth*(Param_cosmo.Sigma_th ...
                            - Param_cosmo.DD_th/(Param_cosmo.L_eth^2))) ; % Eq. 3.39 Gosse & Phillips, 2001

        JDeltaphi_star_eth_a = (Param_cosmo.p_E_th_a*Param_cosmo.R_th_a*FDeltaphi_star_eth_a)/((1/Param_cosmo.Sigma_eth_a)*(Param_cosmo.Sigma_th_a - Param_cosmo.D_th_a/(Param_cosmo.L_eth_a^2))) ;
        
        phi_star_th_a = (Param_cosmo.p_E_th_a*Param_cosmo.R_th_a*phi_star_eth_a)/(1/Param_cosmo.Sigma_eth_a*(Param_cosmo.Sigma_th_a - Param_cosmo.D_th_a/(Const_cosmo.Lambda_f_e^2))) ; 

        Deltaphi_star_th = phi_star_th_a - phi_star_th ; % difference between hypothetical equilibrium thermal neutron fluxes in atmosphere and ss

        JDeltaphi_star_th = (Param_cosmo.D_th_a*(phi_star_th_a/Const_cosmo.Lambda_f_e - JDeltaphi_star_eth_a/Param_cosmo.L_eth_a) ...
                             - Param_cosmo.DD_th*(phi_star_th/Const_cosmo.Lambda_f_e + JDeltaphi_star_eth/Param_cosmo.L_eth) ...
                                + (Param_cosmo.D_th_a/Param_cosmo.L_th_a)*(Deltaphi_star_th + JDeltaphi_star_eth_a - JDeltaphi_star_eth)) ...
                            / ((Param_cosmo.DD_th/Param_cosmo.L_th) + (Param_cosmo.D_th_a/Param_cosmo.L_th_a)) ; % portion of difference between phi_star_th,ss and actual flux due to thermal flux profile

        % Thermal neutron production 
        phi_th = sf.currentsf.SFth .* sf.S_T .* (phi_star_th.*exp(-e./Const_cosmo.Lambda_f_e) + ...
                (1 + R_prime_mu).*JDeltaphi_star_eth.*exp(-e./Param_cosmo.L_eth) + ...
                (1 + R_prime_mu.*Param_cosmo.R_th).*JDeltaphi_star_th.*exp(-e./Param_cosmo.L_th) + ...
                R_prime_mu_depth.*phi_star_th) ; % Thermal neutron flux (n.cm_2.a-1)
            
        % Radiogenic production
        P_rad = Param_cosmo.P_th_r .* f_th + Param_cosmo.P_eth_r .* f_eth ;
            
        % Production-rate
            P_eth = (f_eth ./ Param_cosmo.Lambda_eth) .* phi_eth .* (1 - Param_cosmo.p_E_th);
            P_th =(f_th./Param_cosmo.Lambda_th).*phi_th;       
            
        % keep it in memory     
            xx(i) = P_eth(1);
            yy(i) = P_th(1);
            zz(i) = P_rad(1);
    end
            % Propagated uncertainty from Cl35 on the final P_eth and P_th
             P_eth_35Cl_uncert = std(xx); % uncertainty on P_eth from the [Cl] uncertainty
             P_th_35Cl_uncert = std(yy); % uncertainty on P_eth from the [Cl] uncertainty
             P_rad_uncert = std(zz); % uncertainty on P_rad from the [Cl] uncertainty
    
             
   % Compare uncertainties obtain from 35Cl and from P_f_0, and take the largest        
   if(P_eth_35Cl_uncert>P_eth_uncert); P_eth_uncert=P_eth_35Cl_uncert;end % pick up the largest uncertainty on P_eth
   if(P_th_35Cl_uncert>P_th_uncert); P_th_uncert=P_th_35Cl_uncert; end % pick up the largest uncertainty on P_th

     
 %% Total uncertainty
    % total production rate
    P36_tot = P36_s_Ca + P36_s_K + P36_s_Ti + P36_s_Fe + P36_eth + P36_th + P36_mu + P36_rad;
    % total production rate uncertainty
    P36_tot_uncert = ( (P36_s_uncert(1)).^2 + (P_eth_uncert(1)).^2 ...
                       + (P_th_uncert(1)).^2 + (P36_mu(1).*0.3).^2 ...
                       + (P_rad_uncert(1)).^2 ).^.5;
    % get uncertainty in percent relatively to the production rate
    P36_tot_uncert = P36_tot_uncert(1)./P36_tot(1); 
    

end

