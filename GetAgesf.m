function [ ] = GetAgesf( name )
% Function of the crep program that calculates exposure ages

% Load packages
%pkg load parallel

% Open paths and load data
addpath Functions
addpath Constants
addpath jsonlab
load('Constants/GMDB.mat');
load('Constants/OtherCst.mat');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load parallel
    pkg load statistics
    name=strrep(name,"\\","/"); %OCTAVE !!!!!
else
    name=strrep(name,'\\','/'); %MATLAB !!!!!
end

% Retreive json file
%Data=loadjson(name);
% test the structure
%if(isfield(Data,'Name')==0); Data = Data{1,1}; end;
load Data_36;
%load Data_be;
%load data_fault_1s;
% Load the parameters
Nucl=Data.Nucl;
Scheme=Data.Scheme;
Atm=Data.Atm;
if Atm==0;
    load ERA40
else
    meanP=[];
    meanT=[];
    ERA40lat=[];
    ERA40lon=[];
end


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
% mu=OtherCst.Dens/OtherCst.Attlgth; % Mu becomes a vector with the user
% input
Attlg=OtherCst.Attlgth;
tBe=OtherCst.tBe;
Lambda10Be=log(2)/tBe;
SelPR=Data.PR;

% Set the waitbar
% h = waitbar(0,'Calculating');

% Get user data in right form
VecName=Data.Samples';
VecLat=Data.Lat';
[NbSpl,~]=size(VecLat);
VecLon=Data.Lon';
VecAlt=Data.Alt';
VecConc=Data.NuclCon';
VecErrConc=Data.NuclErr';
VecThick=Data.Thick';
VecShield=Data.Shield';
VecErosion=Data.Eros';
VecDens=Data.Dens';
VecMu=VecDens./Attlg;

% Formatting Outputs
ExitMat=zeros(NbSpl,4);
CellPDF=cell(1,2*NbSpl);
ErrCol=zeros(NbSpl,1);
StatCell=cell(NbSpl,1);


% waitbar(1/(NbSpl+3))

%% 36 Cl Stuff
% Get 36Cl parameters 
if Nucl==36
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

    
    % Variable and Data Initialization 
    [Data_formatted,Param_site,Const_cosmo] = Init_var(Data,flag);

    % Scaling factors initialization
    %w = 0.2; % water content for Sato & Niita (2006)
    w = -1; % water content =  default value
            
    % geomagnetic database
    flag.NumGMDB = NumGMDB;
    Sf = Func_Sf(Param_site,Data,Atm,w,SelGMDB,flag);
    
    % Production rates and constants initialization
    Param_cosmo = clrock(Data_formatted,Param_site,Const_cosmo,Sf);
    
    [check36] = sanity_check_36Cl(Param_cosmo,Sf);
    
    %% Lets calculate ages
    for is=1:NbSpl
        % Create de dataset
        dataset(1,:) = Data.NuclCon(:);
        dataset(2,:) = Data.NuclErr(:);
        dataset(3,:) = Data.Z(:);
        
        % Fix erosion rate
        erosion = Data.Eros;

        % Search parameters
        flag.min_bounds = 0.0; % minimum bound for the search
        flag.max_bounds = 400000; % maximum bound for the search
        flag.search = 'fminsearch'; %flag.search = 'nmsmax';
                 
        % First guess considering the sample belongs to a surface without erosion.
        %t_guess = -log(1-(Data.NuclCon(is)-Param_cosmo{is}.N36Cl.rad-Param_cosmo{is}.N36Cl.inh)*Const_cosmo.lambda36/Param_cosmo{is}.P_cosmo(1))/Const_cosmo.lambda36;
        t_guess = -log(1-(Data.NuclCon(is))*Const_cosmo.lambda36/Param_cosmo{is}.P_cosmo(1))/Const_cosmo.lambda36;

        %t_guess = 10000;randi([flag.min_bounds flag.max_bounds],1,1);

            n_trial = 20;
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
            
            resolution = 'parallel'
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
                    h = waitbar(k/n_trial)
                end    
                    close(h)
           end
            end
           
           Age_MC = mean(age)./1000;
           Err_MC = std(age)./1000;
           
           % best age
           dataset(1,1) = Data.NuclCon(is);
           dataset(2,1) = Data.NuclErr(is);
           dataset(3,1) = Data.Z(is);
           AgeCosmo = Find_age(Const_cosmo, Param_cosmo, Param_site, Sf{is}, dataset, erosion, t_guess, flag)
           
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
                AgeCosmo = AgeCosmo./1000;
                ErrCosmo = ErrCosmo./1000;
             
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
else

%% Calculating
for is=1:NbSpl;
    % waitbar((i+1)/(NbSpl+3))
    ErrPRflag=0;
    % Preliminary corrections
    ThickCorr=(Attlg/(VecThick(is)*VecDens(is)))*(1-exp(-1*(VecThick(is)*VecDens(is))/Attlg));
    CorrConc=VecConc(is)/(ThickCorr*VecShield(is));
    
    if Scheme==1; % Lal-Stone Scaling
        
        % CosmoAge calculation
        StoneFactor=StoneFactV2(VecLat(is),VecLon(is),VecAlt(is),Atm);
        if Nucl==3;
            if VecErosion(is)==0;
                CosmoAge=CorrConc/(SelPR(1)*1000*StoneFactor);
            else
                CosmoAge=(-1/(VecMu(is)*VecErosion(is)))*log(1-(VecMu(is)*VecErosion(is)*CorrConc/(SelPR(1)*StoneFactor)))/1000;
                if isreal(CosmoAge)==0;
                    ErrPRflag=1;
                    CosmoAge=-100;
                end
            end
        elseif Nucl==10;
            CosmoAge=(-1/(Lambda10Be+VecMu(is)*VecErosion(is)))*log(1-((Lambda10Be+VecMu(is)*VecErosion(is))*CorrConc/(SelPR(1)*StoneFactor)))/1000;
            if isreal(CosmoAge)==0;
                ErrPRflag=1;
                CosmoAge=-100;
            end
        end
        ErrCosmo=CosmoAge*VecErrConc(is)/VecConc(is);
        
        % Time-dependent integration-correction
        [Age,Err,Err2,X,Y]=AgeCosmoAgeReelV22(CosmoAge,ErrCosmo,(SelPR(2)/SelPR(1)),VecLat(is),VecLon(is),VecAlt(is),SelGMDB,Atm,ERA40lat,ERA40lon,meanP,meanT);
        
    else % LSD Scaling
        
        if NumGMDB==3;
            % CosmoAge calculation
            [~,SF]=LSDv9(VecLat(is),VecLon(is),VecAlt(is),Atm,0,-1,Nucl,NumGMDB);
            if Nucl==3;
                if VecErosion(is)==0;
                    CosmoAge=CorrConc/(SelPR(1)*1000*SF);
                else
                    CosmoAge=(-1/(VecMu(is)*VecErosion(is)))*log(1-(VecMu(is)*VecErosion(is)*CorrConc/(SelPR(1)*SF)))/1000;
                    if isreal(CosmoAge)==0;
                        ErrPRflag=1;
                        CosmoAge=-100;
                    end
                end
            else
                CosmoAge=(-1/(Lambda10Be+VecMu(is)*VecErosion(is)))*log(1-((Lambda10Be+VecMu(is)*VecErosion(is))*CorrConc/(SelPR(1)*SF)))/1000;
                if isreal(CosmoAge)==0;
                    ErrPRflag=1;
                    CosmoAge=-100;
                end
            end
            ErrCosmo=CosmoAge*VecErrConc(is)/VecConc(is);
            % Time corrected age calculation
            [Age,Err,Err2,X,Y] = AgeCosmoAgeReelLSD(CosmoAge,ErrCosmo,(SelPR(2)/SelPR(1)),VecLat(is),VecLon(is),VecAlt(is),NumGMDB,Atm,Nucl);
        else
            % CosmoAge calculation
            [~,SF]=LSDv9(VecLat(is),VecLon(is),VecAlt(is),Atm,0,-1,Nucl,SelGMDB);
            if Nucl==3;
                if VecErosion(is)==0;
                    CosmoAge=CorrConc/(SelPR(1)*1000*SF);
                else
                    CosmoAge=(-1/(VecMu(is)*VecErosion(is)))*log(1-(VecMu(is)*VecErosion(is)*CorrConc/(SelPR(1)*SF)))/1000;
                    if isreal(CosmoAge)==0;
                        ErrPRflag=1;
                        CosmoAge=-100;
                    end
                end
            else
                CosmoAge=(-1/(Lambda10Be+VecMu(is)*VecErosion(is)))*log(1-((Lambda10Be+VecMu(is)*VecErosion(is))*CorrConc/(SelPR(1)*SF)))/1000;
                if isreal(CosmoAge)==0;
                    ErrPRflag=1;
                    CosmoAge=-100;
                end
            end
            ErrCosmo=CosmoAge*VecErrConc(is)/VecConc(is);
            % Time corrected age calculation
            [Age,Err,Err2,X,Y] = AgeCosmoAgeReelLSD(CosmoAge,ErrCosmo,(SelPR(2)/SelPR(1)),VecLat(is),VecLon(is),VecAlt(is),SelGMDB,Atm,Nucl);
        end
    end
    
    ExitMat(is,2)=Age;
    ExitMat(is,3)=Err;
    ExitMat(is,4)=Err2;
    CellPDF{2*is-1}=X;
    CellPDF{2*is}=Y;
    
    % Look if there is a problem with the length of the database
    if isempty(X)==1;
        ErrCol(is)=1;
        if ErrPRflag==1;
            Mess=sprintf('Production rate too low');
        elseif isnan(Age)==1;
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
    
    % Store the time correctd scaling factor
    if Nucl==3;
        SFtc=CorrConc/(SelPR(1)*Age*1000);
    else
        SFtc=CorrConc*Lambda10Be/(SelPR(1)*(1-exp(-Lambda10Be*Age*1000)));
    end
    ExitMat(is,1)=SFtc;
end
end

MatPDF=CreeMatCP4(CellPDF,0.01);
ColSum=sum(MatPDF(:,2:end),2);
ToExit= ColSum==0;
MatPDF(ToExit,:)=[];
% waitbar((NbSpl+2)/(NbSpl+3))

% Fill PDFtable
VecName2=cell(1,NbSpl);
ifill=1;
for is=1:NbSpl;
    if ErrCol(is)==0;
        VecName2{ifill}=VecName{is};
        ifill=ifill+1;
    end
end
VecName2(ifill:end)=[];

% Density functions to be displayed
PDFdisp=cell(1,(ifill-1));
nbsigdisp=6;
nbptdisp=500;

if ifill>=2;
    for is=1:ifill-1;
        [~,~,Mediane,Err1SigInf,Err1SigSup,~,~]=ParamDistribV11(CellPDF{2*is-1},CellPDF{2*is});
        Binf=max(0,max((Mediane-nbsigdisp*Err1SigInf),CellPDF{2*is-1}(1)));
        Bsup=min(Mediane+nbsigdisp*Err1SigSup,CellPDF{2*is-1}(end));
        ti=linspace(Binf,Bsup,nbptdisp);
        pdfi=interp1(CellPDF{2*is-1},CellPDF{2*is},ti);
        PDFdisp{is}=[ti' pdfi'];
    end
else
    PDFdisp=[];
end

% Prepare Output
DataOut.Samples=Data.Samples;
DataOut.ScalFact=ExitMat(:,1)';
DataOut.Ages=ExitMat(:,2)';
DataOut.AgesErr=ExitMat(:,3)';
DataOut.AgesErr2=ExitMat(:,4)';
DataOut.Status=StatCell';
DataOut.TimeV=MatPDF(:,1)';
DataOut.SamplePDF=VecName2;
DataOut.PDF=MatPDF(:,2:end)';
DataOut.PDFdisp=PDFdisp;

% Write json
DataOut=savejson(name,DataOut);
NameOut=strcat(name(1:end-2),'out');
fileID=fopen(NameOut,'w');
fprintf(fileID,'%s',DataOut);
fclose(fileID);

% Closing
%waitbar(1/1)
% close(h)
end

