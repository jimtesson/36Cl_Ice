function [AgeReel,ErrReel,ErrReel2,XPDFReel,PDFReel] = AgeCosmoAgeReelV22(AgeCosmo,ErrCosmo,ErrRelPR,Latitude,Longitude,Altitude,PaleoMag,Atm,ERA40lat,ERA40lon,meanP,meanT)
%--------Description------------------------------------------------------
% This functions applies the geomagnetic time correction on a CRE ages
% using the 'Lal modified' model (Balco et al., 2008; Lal, 1991; Nishiizumi
% et al., 1989; Stone, 2000). For a given geomagnetic database it computes
% a calibration curve to convert the gaussian Probability Density Function
% (PDF) associated with the non-corrected ages into an indeterminate
% corrected PDF. This function return a corrected age +/- 1-sigma and the
% corresponding Probability Density Function (PDF). Yet the Age value
% include the production rate uncertainty used for preliminary calculations
% whereas the PDF does not. This is because the PDF are proposed for
% dispersion discussion on ages acounting for the same object.
%
% Input : 
%        AgeCosmo  : Preliminary age with no geomagnetic corrections (in ka)                
%        ErrCosmo  : 1-sigma associated uncertainty (in ka)
%        ErrRelPR  : 1-sigma uncertainty associated to the production rate used for preliminary calculation (in at.g-1.yr-1)
%        Latitude  : Sample latitude (in decimal degrees, <0 if S)
%        Longitude : Sample longitude (in decimal degrees, <0 if W)
%        Altitude  : Sample elevation (in masl)
%        Paleomag  : Geomagnetic Virtual Dipolar Moment (VDM) reconstruction (2-lines matrix: L1= ages in ka, L2=VDM in 10^2 A.m-2)
%        Atm       : Atmosphere model, 0 : ERA-40
%                                      1 : Atmosphere standard.
%        ERA40lat, ERA40lon, meanP, meanT : Input data to use the ERA40atmosphere model
%
% Output  : 
%        AgeReel   : Past geomagnetic activity corrected ages (ka)
%        ErrReel   : 1-sigma associated uncertainty (ka) (production rate uncertainty is considered)
%        XPDFReel  : Vector of time values for the associated PDF (ka)
%        PDFReel   : Vector of density function values for the associated PDF (production rate uncertainty is not considered)
%
% IMPORTANT : Requires Matlab 2009 or any more recent versions.

% Code written by LCP Martin, PH Blard and J Lav?
% Centre de Recherches P?trographiques et G?ochimiques (CRPG-CNRS), France
% blard@crpg.cnrs-nancy.fr
% Program desciprtion provided in Martin et al., (In Prep)
%
% This code contains portions of the code of N Lifton (Lifton et al., 2014,
% under GNU licence)
% 
% Copyright 2015, CNRS-Universit? de Lorraine
% All rights reserved
%
% This file is part of the CREp program.
% CREp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>
%-------------------------------------------------------------------------

%--------Basic Stop-------------------------------------------------------
if AgeCosmo<0;
    AgeReel=NaN;
    ErrReel=NaN;
    XPDFReel=[];
    PDFReel=[];
    ErrReel2=NaN;
    return
end

%--------VDM Importation--------------------------------------------------
VecAgeReel1=PaleoMag(1,:);
PaleoVDM=PaleoMag(2,:);
[~,Dates]=size(VecAgeReel1);

ErrCosmo2=AgeCosmo*sqrt((ErrCosmo/AgeCosmo)^2+(ErrRelPR)^2);
if VecAgeReel1(end)>(AgeCosmo+10*ErrCosmo2); % Chope if too long
    VecIndice=find(VecAgeReel1>(AgeCosmo+10*ErrCosmo2));
    VecAgeReel1=VecAgeReel1(1:VecIndice(1));
    PaleoVDM=PaleoVDM(1:VecIndice(1));
    Dates=length(VecAgeReel1);
end

%--------Niishizumi et al.,(1989) correction------------------------------
% Computation of cutoff rigidity Rc(t) (Dunai 2001, equation 1 and Lifton et al 2008):
Rc = (PaleoVDM*1e22*4*1e-7*3*1e8)/(16*1e9*(6.3712*1e6)^2)*(cos(pi*Latitude/180))^4;

% Cos(LambdaM)
%CosLambdaM=(PaleoVDM/PaleoVDM(1)).^0.25*cos(pi*Latitude/180);

% Find LambdaM
%LambdaM=zeros(1,Dates);
%for k=1:length(CosLambdaM);
 %   if CosLambdaM(k)<=1;
  %      LambdaM(1,k)=acos(CosLambdaM(k))*180/pi;
  %  else
  %      LambdaM(1,k)=0.01;
  %  end
%end

%--------Use Stone (2000)-------------------------------------------------
% The implementation the atmosphere models used here are from Lifton et al.
% (2014).
StoneFact=zeros(1,Dates);
if Atm==1;
    gmr = -0.03417;
    dtdz = 0.0065;
    P=1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (Altitude.*dtdz)) ) );

        StoneFact=StoneFactorL(Rc,P,1013.25);

elseif Atm==0;
    if Longitude<0;
        Longitude=Longitude+360;
    end
    site_slp = interp2(ERA40lon,ERA40lat,meanP,Longitude,Latitude);
    site_T = interp2(ERA40lon,ERA40lat,meanT,Longitude,Latitude);
    gmr = -0.03417;
    lr = [-6.1517E-03 -3.1831E-06 -1.5014E-07 1.8097E-09 1.1791E-10 ...
        -6.5359E-14 -9.5209E-15];
    dtdz = lr(1) + lr(2).*Latitude + lr(3).*Latitude.^2 ...
        + lr(4).*Latitude.^3 + lr(5).*Latitude.^4 + lr(6).* Latitude.^5 ...
        + lr(7).*Latitude.^6;
    dtdz = -dtdz;
    P=site_slp .* exp( (gmr./dtdz) .* ( log(site_T) - log(site_T - (Altitude.*dtdz)) ) );
     
    StoneFact=StoneFactorL(Rc,P,site_slp);
  
end

%--------Time average of the scaling factors------------------------------
StoneFactMoy=MoyenneIntegrV2(VecAgeReel1,StoneFact);


%--------Normalisation to present day scaling factor----------------------
StoneFactMoyNorm=zeros(1,Dates);
for k=1:length(StoneFactMoy);
    StoneFactMoyNorm(1,k)=StoneFactMoy(k)/StoneFact(1);
end


%Max(StoneFactMoyNorm(1,:)

%--------Calibration curve------------------------------------------------
CourbeAgeCosmo1=StoneFactMoyNorm.*VecAgeReel1;

%--------Interpolation of the created vectors------------------------------

% Determine time step
AgeDiscr=0.01;


%% Create corresponding vector
%CourbeAgeCosmo2 = logspace(0,CourbeAgeCosmo1(end),100); TRY THIS. MODIFY
%CMAP CREP4.
CourbeAgeCosmo2=0:AgeDiscr:CourbeAgeCosmo1(end);
% Interpolate

VecAgeReel2=interp1(CourbeAgeCosmo1,VecAgeReel1,CourbeAgeCosmo2,'spline');

%--------Basic stopping conditions----------------------------------------
if AgeCosmo>CourbeAgeCosmo2(end);
    AgeReel=NaN;
    ErrReel=NaN;
    XPDFReel=[];
    PDFReel=[];
    ErrReel2=[];
    return
end

%------------------Calibration of gaussian ages---------------------------

% For the age value, the production rate uncertainty is considered
[AgeReel,ErrReel,~,~]=CalibGauss(AgeCosmo,ErrCosmo2,VecAgeReel2,CourbeAgeCosmo2,AgeDiscr);
% For the PDF values, the production rate uncertainty is not considered
[~,ErrReel2,XPDFReel,PDFReel]=CalibGauss(AgeCosmo,ErrCosmo,VecAgeReel2,CourbeAgeCosmo2,AgeDiscr);

end