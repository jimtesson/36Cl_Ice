function [ Sf_out ] = Func_Sf(Param_site_in,Data,atm,w,Paleomag,flag)
% Func_Sf: Function to calculate the scaling factors over a time vector,
% based on Marrero et al. 2016. The spacing of the time vector is variable:
% in the Holocene, it follows the spacing of the paleomagnetic data records, 
% mostly at 500 yrs. Between 12000 and 800,000 yr, the spacing is 1000 yr 
% to agree with the SINT-800 record. After 800,000 yt there is logarithmic
% spacing out to 10 million years. See the hard-copy documentation for 
% Marrero et al. 2016 for more details. 

% INPUT: 
%       scaling_model :
%           scaling_model = 'user' : user-provided scaling factors for
%                                   neutrons and muons.
%           scaling_model = 'st' : stone 2000 scheme scaling factors
%           scaling_model = 'sa' : Lifton-Sato scheme scaling factors
%           scaling_model = 'all' : computes stone 2000 and Lifton-Sato
%
%       Param_site :  site-specific parameters.
% OUTPUT:
%       out : struct. containing the scaling factors.
%

NbSpl = numel(Param_site_in); % number of samples

for i=1:NbSpl;
    Param_site = Param_site_in{i};
    % site parameters
    lat = Param_site.lat;
    long = Param_site.long;
    alt = Param_site.alt;

    % Correction factors
    out.S_T = Data.Shield(1); % correction factor for shielding of a sample of arbitrary orientation by surrounding topography
    %Sf.S_T_er = 0.01*Data.Shield(1); % uncertainties
    out.S_snow = 1; % correction factor for snow shielding for spallogenic production
    %Sf.S_snow_er = 0.01; % uncertainties
    out.S_shape = 1; % correction factor for geometry effects on spallogenic production
    %Sf.S_shape_er = 0.03; % uncertainties

    NumGMDB = flag.NumGMDB;

%% Lal-Stone Scaling using cutoff rigidity 

%if(strcmp(flag.scaling_model,'st'))
        %--------VDM Importation--------------------------------------------------
        VecAgeReel1=Paleomag(1,:);
        PaleoVDM=Paleomag(2,:);

        %--------Niishizumi et al.,(1989) correction------------------------------
        % Computation of cutoff rigidity Rc(t) (Dunai 2001, equation 1 and Lifton et al 2008):
        Rc = (PaleoVDM*1e22*4*1e-7*3*1e8)/(16*1e9*(6.3712*1e6)^2)*(cos(pi*lat/180))^4;

        %--------Use Stone (2000)-------------------------------------------------
        % The implementation the atmosphere models used here are from Lifton et al.
        % (2014).

        if atm==1;
                gmr = -0.03417;
                dtdz = 0.0065;
                P=1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (alt.*dtdz)) ) );

            [SF_St,SF_St_sp,SF_St_mu]=StoneFactorL(Rc,P,1013.25);

        elseif atm==0;
            load ERA40
            if long<0;
                long=long+360;
            end
            site_slp = interp2(ERA40lon,ERA40lat,meanP,long,lat);
            site_T = interp2(ERA40lon,ERA40lat,meanT,long,lat);
            gmr = -0.03417;
            lr = [-6.1517E-03 -3.1831E-06 -1.5014E-07 1.8097E-09 1.1791E-10 ...
                    -6.5359E-14 -9.5209E-15];
            dtdz = lr(1) + lr(2).*lat + lr(3).*lat.^2 ...
                + lr(4).*lat.^3 + lr(5).*lat.^4 + lr(6).* lat.^5 ...
                + lr(7).*lat.^6;
            dtdz = -dtdz;
            P=site_slp .* exp( (gmr./dtdz) .* ( log(site_T) - log(site_T - (alt.*dtdz)) ) );
     
            [SF_St,SF_St_sp,SF_St_mu]=StoneFactorL(Rc,P,site_slp);
            
        end
        
        %[SF_St2,SF_St_sp2,SF_St_mu2] =StoneFactV2(lat,long,alt,atm);
        out.SF_St = SF_St;
        out.SF_St_sp = SF_St_sp;
        %out.SF_St_sp_er = SF_St_sp .* 0.05;
        out.SF_St_mu = SF_St_mu;
        %out.SF_St_mu_er = SF_St_mu .* 0.05;
        out.SF_Vectime = VecAgeReel1.*1000; % GMDB in ka but LSD in a
        
%%  Lifton-Sato        
%elseif(strcmp(flag.scaling_model,'sa') )      

    % load constants
    load('pmag_consts.mat')
    load('Constants/consts_LSD.mat')
    consts = pmag_consts; clear pmag_consts;
    
    % Pressure Model
    if atm == 1
        stdatm = 1;
        gmr = -0.03417; % Assorted constants
        dtdz = 0.0065; % Lapse rate from standard atmosphere
    else
        stdatm = 0;
    end
    % Pressure correction
    if stdatm == 1
        % Calculate site pressure using the Standard Atmosphere parameters with the
        % standard atmosphere equation.
        pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (alt.*dtdz)) ) );
    else    
        pressure = ERA40atm(lat,long,alt);
    end
    
    % What type of paleomagnetic data
    if NumGMDB==3
        % Particular case of LSD framework
        GeoMDBNorm(1,:)=consts.t_M;
        GeoMDBNorm(2,:)=consts.M;
        Mt0=[0;1];
        GeoMDBNorm=[Mt0 GeoMDBNorm];
    else
        % Case of classical VDM data
        Paleomag(2,:)=Paleomag(2,:)./Paleomag(2,1);
        GeoMDBNorm=Paleomag;
        GeoMDBNorm(1,:)=GeoMDBNorm(1,:).*1000; % GMDB in ka but LSD in a
    end

    % catch for negative longs before Rc interpolation
    if long < 0; long = long + 360;end;

    %% Parameters used to compute scaling factors over time 
    
    % Time vector for time dependent scaling factors: Age Relative to t0=2010
    tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % Need solar modulation for Lifton SF's
    this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
    this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)
    
    % water content
    if w < 0
        w = 0.066; % default gravimetric water content for Sato et al. (2008)
    end
    
    % interpolate an M for tv > 7000...
    %first sort the values to be interpolated for
    %[sorttv,indextv]=sort(tv(77:end));
    
    %   New Equation - unpublished -  Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average. This is the one I'm using now...
    dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
    %use "interpolate" to interpolate for the sorted values
    VecM = interp1(GeoMDBNorm(1,:),GeoMDBNorm(2,:),tv);
    
    % Make up the Rc vectors.
    LSDRc = zeros(1,length(tv));
    LSDRc(1:end) = VecM.*(dd(1)*cosd(lat) + ...
       dd(2)*(cosd(lat)).^2 + ...
       dd(3)*(cosd(lat)).^3 + ...
       dd(4)*(cosd(lat)).^4 + ...
       dd(5)*(cosd(lat)).^5 + ...
       dd(6)*(cosd(lat)).^6);    
    % Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
    if Paleomag==3;
        [longi,lati,tvi] = meshgrid(long,lat,tv(1:76));
        LSDRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,longi,lati,tvi);
    end

    if (min(LSDRc)<0)
        %        error('LiRc contains negative rigidity cutoffs');
        LSDRc(LSDRc<0)=0.0;
    end
    
    % Get scaling factors
	scaling=LiftonSatoSX(pressure,LSDRc,this_SPhi,w,consts);
	out.SF_Sf=scaling.sp;
	out.Rc_Sf=LSDRc;
	out.SF_Sa36Ca=scaling.ClCa;
	out.SF_Sa36K=scaling.ClK;
	out.SF_Sa36Ti=scaling.ClTi;
	out.SF_Sa36Fe=scaling.ClFe;
	out.SF_Saeth=scaling.eth;
	out.SF_Sath=scaling.th;
	out.Rc_Sa=LSDRc;
    out.tv = tv;
    
Sf_out{i} = out;
end      

