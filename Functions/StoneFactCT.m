function [FactStoneCorrige] = StoneFactCT(AgeFact,Latitude,Longitude,Altitude,PaleoMag,Atm,ERA40lat,ERA40lon,meanP,meanT)

%--------Description------------------------------------------------------
% This functions computes Lal-Stone scaling factors (Lal, 1991; Stone 2000)
% with the geomagnetic time correction proposed by Nishiizumi et al. (1989)
% as described in Balco et al. (2008). It is implemented for two atmosphere
% models. Different Virtual Dipolar Moment reconstruction can be input if
% they fit the required format.
%
% Input : 
%        AgeFact   : Age of the calibration object (in ka)      
%
%        Latitude  : Latitude of the calibration object (in decimal
%        degrees, <0 if S)-Note that latitude is latter transformed into a
%        cutoff rigidity to make sure that cutoff rigidites larger than
%        14.4 GeV (present day cutoff rigidity at equator) can be handled
%
%        Longitude : Longitude of the calibration object (in decimal degrees, <0 if W)
%        Altitude  : Elevation (in masl)
%        Paleomag  : Geomagnetic Virtual Dipolar Moment (VDM) reconstruction (2-lines matrix: L1= ages in ka, L2=VDM in 10^2 A.m-2)
%        Atm       : Atmosphere model, 0 : ERA-40 (Uppala et al., 2005)
%                                      1 : US Standard atmosphere (National Oceanic and Atmospheric Administration, 1976)
%        ERA40lat, ERA40lon, meanP, meanT : Input data to use the ERA40atmosphere model
%
% Output  : 
%        FactStoneCorrige : Lal-Stone scaling factor accounting for geomagnetic correction
%
% IMPORTANT : Requires Matlab 2009 or any more recent versions.

% Code written by LCP Martin, PH Blard and J Lavé
% Centre de Recherches Pétrographiques et Géochimiques (CRPG-CNRS), France
% blard@crpg.cnrs-nancy.fr
% Program desciprtion provided in Martin et al., (In Prep)
% 
% Copyright 2015, CNRS-Université de Lorraine
% All rights reserved
%
% This file is part of the CREp program.
% CREp is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. See <http://www.gnu.org/licenses/>
%-------------------------------------------------------------------------


%--------VDM Importation--------------------------------------------------
VecAgeReel1=PaleoMag(1,:);
PaleoVDM=PaleoMag(2,:);
[~,Dates]=size(VecAgeReel1);

% Block if the age is too old
if AgeFact>VecAgeReel1(end-2);
    FactStoneCorrige=NaN;
    return
end

% Chope if the geomagnetic database is too long
if VecAgeReel1(end)>(1.5*AgeFact); 
    VecIndice=find(VecAgeReel1>(1.5*AgeFact));
    VecAgeReel1=VecAgeReel1(1:VecIndice(1));
    PaleoVDM=PaleoVDM(1:VecIndice(1));
    Dates=length(VecAgeReel1);
end

%--------Niishizumi et al.,(1989) correction------------------------------

% Cos(LambdaM)
% Computation of cutoff rigidity Rc(t) (Dunai 2001, equation 1 and Lifton et al 2008):
Rc = (PaleoVDM*1e22*4*1e-7*3*1e8)/(16*1e9*(6.3712*1e6)^2)*(cos(pi*Latitude/180))^4;




% Find LambdaM
%LambdaM=zeros(1,Dates);
%for k=1:length(CosLambdaM);
 %   if CosLambdaM(k)<=1;
  %      LambdaM(1,k)=acos(CosLambdaM(k))*180/pi;
  %  else
   %     LambdaM(1,k)=0.01;
   % end
%end

%--------Use Stone (2000)-------------------------------------------------
% The implementation the atmosphere models used here are from Lifton et al.
% (2014).
StoneFact=zeros(1,Dates);
if Atm==1;
    gmr = -0.03417;
    dtdz = 0.0065;
    P = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (Altitude.*dtdz)) ) );
    
 
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
  StoneFact=StoneFactorL(Rc,P,1013.25);
  
end

%--------Time average of the scaling factors------------------------------
StoneFactMoy=MoyenneIntegrV2(VecAgeReel1,StoneFact);


%--------Find the appropriate factor--------------------------------------
FactStoneCorrige=interp1(VecAgeReel1,StoneFactMoy,AgeFact,'spline');

end