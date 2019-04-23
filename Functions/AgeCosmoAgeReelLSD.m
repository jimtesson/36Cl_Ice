function [AgeReel,ErrReel,ErrReel2,XPDFReel,PDFReel] = AgeCosmoAgeReelLSD(AgeCosmo,ErrCosmo,ErrRelPR,Latitude,Longitude,Altitude,PaleoMag,Atm,nuclide)
%--------Description------------------------------------------------------
% This functions applies the geomagnetic time correction on a CRE ages
% using the LSD model (Lifton Sato Dunai, 2014). For a given geomagnetic
% database it computes a calibration curve to convert the gaussian
% Probability Density Function (PDF) associated with the non-corrected ages
% into an indeterminate corrected age PDF. This function return a corrected
% age +/- 1-sigma and the corresponding Probability Density Function (PDF).
% Yet the Age value include the production rate uncertainty used for
% preliminary calculations whereas the PDF does not. This is because the
% PDF are proposed for dispersion discussion on ages acounting for the same
% object.
%
% Input : 
%        AgeCosmo  : Preliminary age with no geomagnetic corrections (in ka)                
%        ErrCosmo  : 1-sigma asociated uncertainty (in ka)
%        ErrRelPR  : 1-sigma uncertainty associated to the production rate used for preliminary calculation (in at.g-1.yr-1)
%        Latitude  : Sample latitude (in decimal degrees, <0 if S)
%        Longitude : Sample longitude (in decimal degrees, <0 if W)
%        Altitude  : Sample elevation (in masl)
%        Paleomag  : Geomagnetic Virtual Dipolar Moment (VDM) reconstruction (2-lines matrix: L1= ages in ka, L2=VDM in 10^2 A.m-2)
%        Atm       : Atmosphere model, 0 : ERA-40
%                                      1 : Atmosph?re standard.
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

%--------Creation of scaling factors vectors------------------------------
if sum(size(PaleoMag))>2;
   AgeFin=1000*PaleoMag(1,end);
else
   AgeFin=2000000;
end
[VecAgeReel1,StoneFact]=LSDv9(Latitude,Longitude,Altitude,Atm,AgeFin,-1,nuclide,PaleoMag);
VecAgeReel1=VecAgeReel1./1000;
VecAgeReel1(end-2:end)=[];
StoneFact(end-2:end)=[];

%--------Time average of the scaling factors------------------------------
StoneFactMoy=MoyenneIntegrV2(VecAgeReel1,StoneFact);

%--------Normalisation to present day scaling factor----------------------
Dates=length(VecAgeReel1);
StoneFactMoyNorm=zeros(1,Dates);
for k=1:length(StoneFactMoy);
    StoneFactMoyNorm(1,k)=StoneFactMoy(k)/StoneFact(1);
end

%--------Calibration curve------------------------------------------------
CourbeAgeCosmo1=StoneFactMoyNorm.*VecAgeReel1;

%--------Interpolation of the created vectors------------------------------

% Determine time step
AgeDiscr=0.01;
% Create corresponding vector
CourbeAgeCosmo2=0:AgeDiscr:CourbeAgeCosmo1(end);
% Interpolate
VecAgeReel2=interp1(CourbeAgeCosmo1,VecAgeReel1,CourbeAgeCosmo2,'spline');

%--------Basic stopping conditions----------------------------------------
if AgeCosmo>CourbeAgeCosmo2(end);
    AgeReel=NaN;
    ErrReel=NaN;
    ErrReel2=NaN;
    XPDFReel=[];
    PDFReel=[];
    return
end

%------------------Calibration of gaussian ages---------------------------

% For the age value, the production rate uncertainty is considered
ErrCosmo2=AgeCosmo*sqrt((ErrCosmo/AgeCosmo)^2+(ErrRelPR)^2);
[AgeReel,ErrReel,~,~]=CalibGauss(AgeCosmo,ErrCosmo2,VecAgeReel2,CourbeAgeCosmo2,AgeDiscr);
% For the PDF values, the production rate uncertainty is not considered
[~,ErrReel2,XPDFReel,PDFReel]=CalibGauss(AgeCosmo,ErrCosmo,VecAgeReel2,CourbeAgeCosmo2,AgeDiscr);

end