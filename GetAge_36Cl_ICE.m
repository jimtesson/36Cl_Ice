function [ ] = GetAge_36Cl_ICE()
% Function of the crep program that calculates exposure ages

% Loading data
    % ice
    load Data_36_ice;
    Data_ICE = Data;
    clear Data;
    % rock
    load Data_36;

% !!!!!!!!!!!!!!!!! to be given by the user: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Constants
     Data.lambda36 = 2.303e-6 ;
     Data.lambda36_uncert = Data.lambda36 .* 0.0066;

    % Attenuation length
     Data.Lambda_f_e = 160; % effective fast neutron attenuation coefficient (g.cm-2)
     Data.Lambda_mu = 1510 ; % slow muon attenuation length (g.cm-2)
     Data.Psi_mu_0 = 190 ; % slow negative muon stopping rate at land surface (muon/g/an), Heisinger et al. (2002)

    % Unscaled sample specific 36Cl production rate by spallation of target elements
     Data.Psi_Cl36_Ca_0 = 42.2 ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_K_0 = 148.1 ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Ti_0 = 13 ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Fe_0 = 1.9 ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr)
    %uncertainties
     Data.Psi_Cl36_Ca_0_uncert = 4.8 ;% spallation production rate for Ca, SLHL (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_K_0_uncert = 7.8 ;% Spallation production rate at surface of 39K (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Ti_0_uncert = 3 ; % Spallation production rate at surface of Ti (at of Cl36 /g of Ca per yr)
     Data.Psi_Cl36_Fe_0_uncert = 0.2 ; % Spallation production rate at surface of Fe (at of Cl36 /g of Ca per yr) 
% !!!!!!!!!!!!!!!!! to be given by the user: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Open paths and load cst
        addpath Functions
        addpath Constants
        load('Constants/GMDB.mat');
        load('Constants/OtherCst.mat');

        Atm=Data.Atm; % atmospheric model
        
    % Geomagnetic data base
    if length(Data.GMDB)==1;
        %NumGMDB=Data.GMDB;
        NumGMDB = 2; % 1: Mush; 2: GLOPIS; 3: LSD; 4: own user geomagnetic db
        if NumGMDB==1;
            SelGMDB=GMDB.Musch;
        elseif NumGMDB==2;
            SelGMDB=GMDB.GLOPIS;
        else %  (NumGMDB=3)
            SelGMDB=GMDB.LSD;
        end
    else
        NumGMDB=4;
        SelGMDB=Data.GMDB;
    end

    % Get number of samples
    VecLat=Data.Lat';
    [NbSpl,~]=size(VecLat);

    % Formatting Outputs
    ExitMat=zeros(NbSpl,4);
    CellPDF=cell(1,2*NbSpl);
    ErrCol=zeros(NbSpl,1);
    StatCell=cell(NbSpl,1);

% Get 36Cl parameters 
    addpath(genpath('Functions/36Cl_Functions'));
    
   % Model used for the modeling of 36Cl concentrations
   
     Scheme=1; % scaling model (1: LAL-STONE, 2: LSD)
     Muon_model = 2; % 1: Exponential, 2: numeric integration of the flux (Balco 2008)
            % Muon
            if(Muon_model == 1)
                flag.muon = 'exp';  % Muon attenuation length approximated by exponential (Schimmelfenning 2009).  
           elseif(Muon_model == 2)
                flag.muon = 'num'; % Muon attenuation length calculated following Heisinger (2002a,b) and Balco et al. (2008)
            else
                Mess=sprintf('Invlid muon model');
            end
            % Scaling model
            if(Scheme == 1)
                flag.scaling_model = 'st'; % LAL-STONE scheme
            elseif(Scheme == 2)
                flag.scaling_model = 'sa'; % LSD scheme  
            else
                Mess=sprintf('Invlid scaling model');
            end  
    % Constants initialization
        [Const_cosmo] = Init_const(Data,flag);   
        
    % Variable and Data Initialization 
        % rock
        [Data_formatted,Param_site] = Init_var(Data);
        % ice
        [Data_formatted_ICE,Param_site_ICE] = Init_var(Data_ICE);

    % Scaling factors initialization
    %w = 0.2; % water content for Sato & Niita (2006)
    w = -1; % water content =  default value
    w_ICE = 1;
            
    % geomagnetic database
    flag.NumGMDB = NumGMDB;
    Sf = Func_Sf(Param_site,Data,Atm,w,SelGMDB,flag);
    Sf_ICE = Func_Sf(Param_site_ICE,Data_ICE,Atm,w_ICE,SelGMDB,flag);
    
    % Production rates and constants initialization
    Param_cosmo = clrock(Data_formatted,Param_site,Const_cosmo,Sf);
    Param_cosmo_ICE = clrock(Data_formatted_ICE,Param_site_ICE,Const_cosmo,Sf_ICE);
    
    %% Lets calculate ages
    for is=1:NbSpl
        % Create de dataset
        dataset(1,:) = Data.NuclCon(:);
        dataset(2,:) = Data.NuclErr(:);
        dataset(3,:) = Data.Z(:);
        
        % Fix erosion rate
        erosion = 1000; Data.Eros;

        % Search parameters
        flag.min_bounds = 0.0; % minimum bound for the search
        flag.max_bounds = 400000; % maximum bound for the search
        flag.search = 'fminsearch'; %flag.search = 'nmsmax';
        flag.plot = 1;
        
        t_expo = 300000;
        t_vector = [0:100:t_expo];
        t_degla = 10000; % LGM
        thick_test = [30 60 100 150 200];
        
        for i=1:length(thick_test)
        ice_thickness = thick_test(i).*100; % in cm
        
        z_history=ones(1,length(t_vector)).*ice_thickness;
        i_degla = find(t_vector<t_degla);
        z_history(i_degla)=t_vector(i_degla).*ice_thickness/t_degla;

                 
        depth_profile_speed(t_vector,z_history,...
                            Const_cosmo,Param_cosmo,Param_cosmo_ICE, ...
                            Param_site,Param_site_ICE,...
                            Sf{is},Sf_ICE{is},erosion,t_expo,Data.Z(:),flag);
        end
                        
                        
        % First guess considering the sample belongs to a surface without erosion.
        %t_guess = -log(1-(Data.NuclCon(is)-Param_cosmo{is}.N36Cl.rad-Param_cosmo{is}.N36Cl.inh)*Const_cosmo.lambda36/Param_cosmo{is}.P_cosmo(1))/Const_cosmo.lambda36;
        t_guess = -log(1-(Data.NuclCon(is))*Const_cosmo.lambda36/Param_cosmo{is}.P_cosmo(1))/Const_cosmo.lambda36;

        %t_guess = 10000;randi([flag.min_bounds flag.max_bounds],1,1);

            n_trial = 100;
            % preparing dataset integrating data uncertainty
            for k=1:n_trial
                for j=1:length(NbSpl)                    
                    dataset(1,j) = Data.NuclCon(j) + Data.NuclErr(j)*randn(1);
                end
                d{k} = dataset;
            end
            % anonymous function
            func = @(ik) Find_age(Const_cosmo, Param_cosmo, Param_site, Sf{is}, d{ik}, erosion, t_guess, flag);
            age = zeros(1,n_trial);
            
            resolution = 'serial';
            % serial resolution
            if(strcmp(resolution,'serial')==1)
                for k=1:n_trial
                    %Best_age = Find_age(Const_cosmo, Param_cosmo, Param_site, Sf, dataset, erosion, t_guess, flag);
                    age(k) = func(k);
                    t_guess = age(k); % starting age for the next search
                    h = waitbar(k/n_trial)
                end
                    close(h)
            elseif(strcmp(resolution,'parallel')==1)
            % parallel resolution
            age = [];           
           if isOctave
               %OCTAVE
                numCores = nproc();
                [age] = pararrayfun (numCores, func, [1:1:n_trial]);
           else
                %Matlab
                parfor k=1:n_trial
                    age(k) =feval(func,k);
                    
                end    
                   
           end
            end
           
           Age_MC = mean(age)./1000
           Err_MC = std(age)./1000
           
           % best age
           dataset(1,1) = Data.NuclCon(is);
           dataset(2,1) = Data.NuclErr(is);
           dataset(3,1) = Data.Z(is);
           AgeCosmo = Find_age(Const_cosmo, Param_cosmo, Param_site, Sf{is}, dataset, erosion, t_guess, flag);
           
           if isreal(AgeCosmo)==0; % check error
               ErrPRflag=1;
               AgeCosmo=-100;
           end
           
           
            % 2. Assume there is some value of P_effective that satisfies the simple exposure age equation. N = (P_effective/lambda)*(1-exp(-lambda*t))
            
                P_effective = Data.NuclCon(is).*Const_cosmo.lambda36./(1-exp(-Const_cosmo.lambda36.*AgeCosmo));
                
            % 3. Estimate the total uncertainty on P_effective by adding up uncertainties on the individual production rates from K, Ca, Cl, and then scaling to whatever the value of P_effective is. Each one of those includes both an uncertainty in the element concentration and an uncertainty in the reference production rate. 
                    % get current scaling factors.
                    t_vector = [0:100:AgeCosmo];
                    currentsf=getcurrentsf(Sf{is},t_vector./1000,flag);
                    Sf{is}.currentsf = currentsf;
                    
                    % Settings to compute uncertainties
                    
                        % Monte Carlo sampling of p_f_0 to obtain the propagated uncertainty on P_th
                        % and P_eth ? yes or no. If no, the uncertainty is 25%
                        flag.uncert.pf0_mc='no';
                        
                    % uncertainty on the Production rate , in %   
                    P36_uncert = get_36Cl_uncert(Data_formatted,Const_cosmo,Param_cosmo{1},Param_site{1},Sf{is},Data.Z(1).*Param_site{1}.rho_rock,flag);
                    
            % 4. Propagate the uncertainty on N (Cl-36 measurement uncertainty) and P_effective (from above) using the normal linear error propagation formula.        
                    ErrCosmo = (((-1/Const_cosmo.lambda36)./(1 - Data_formatted{1}.N36Cl.meas*Const_cosmo.lambda36/P_effective) ...
                        .*(Data_formatted{1}.N36Cl.meas*Const_cosmo.lambda36/(P_effective.^2)))^2).^.5  ...
                        .*P36_uncert.*P_effective; % in year
                
             % yr to kyr
                AgeCosmo = AgeCosmo./1000
                ErrCosmo = ErrCosmo./1000
             
                [XPDFCosmo,PDFCosmo] = PDF_from_Age( AgeCosmo,ErrCosmo );
                [XPDF_MC,PDF_MC] = PDF_from_Age( Age_MC,Err_MC );
                
                figure(1)
                plot(XPDFCosmo,PDFCosmo); hold on
                plot(XPDF_MC,PDF_MC); hold on
                
                figure(2)
                if isOctave
                    hist(age);
                else
                    histogram(age);
                end
                
             % OUTPUT   
                ExitMat(is,1) = .0;
                ExitMat(is,2)=AgeCosmo;
                ExitMat(is,3)=ErrCosmo;
                ExitMat(is,4)=ErrCosmo;
                CellPDF{2*is-1}=XPDFCosmo;
                CellPDF{2*is}=PDFCosmo;
                
              % Look if there is a problem with the length of the database
                if isempty(XPDFCosmo)==1;
                    ErrCol(is)=1;
                    if ErrPRflag==1;
                        Mess=sprintf('Production rate too low');
                    elseif isnan(AgeCosmo)==1;
                        Mess=sprintf('Sample too old for the Geomagnetic database');
                %         elseif Age<2; Chope if to young ?
                %             Mess=sprintf('Sample too young to provide probability density distribution');
                    else
                        Mess=sprintf('Relative error bar too large to provide probability density function (excursion in negative ages or ages older than the Geomagnetic database)');
                    end
                else
                    Mess=sprintf('Ok');
                end
                
                StatCell{is}=Mess;
                
                    
    end % end loop over samples    


