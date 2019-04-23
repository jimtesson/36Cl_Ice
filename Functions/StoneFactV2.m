function [StFact,StFact_sp,StFact_mu] = StoneFactV2( Lat,Lon,Alt,Atm )
% This function computes the time INDdependent scaling factors of Lal-Stone (Lal, 1991;
% Stone, 2000) using either the ERA40 atmosphere (Uppala et al., 2005) or
% the US Standard Atmosphere (National Oceanic and Atmospheric
% Administration, 1976). 
%
% Input : 
%          Lat    : Sample Latitude (decimal degree)
%          Lon    : Sample Longitude (decimal degree)
%          Alt    : Sample elevation (masl)
%          Atm    : 0 - ERA40 Atmosphere
%                   1 - Standard Atmosphere
%
% Output : 
%          StFact : Lal-Stone scaling factor
%
% Code written by LCP Martin, PH Blard and J Lav?
% Centre de Recherches P?trographiques et G?ochimiques (CRPG-CNRS), France
% blard@crpg.cnrs-nancy.fr
% Program desciprtion provided in Martin et al., (In Prep)
%
% This code contains portions of the code from N Lifton for the LSD model
% (Lifton et al., 2014, under GNU licence)
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
%Conversion of latitude in present-day cutoff rigidity, using the 2010 VDM
Rc = 14.3.*((cos(Lat.*pi/180)).^4);


if Atm==0;
    load ERA40
    if Lon<0;
        Lon=Lon+360;
    end
    site_slp = interp2(ERA40lon,ERA40lat,meanP,Lon,Lat);
    site_T = interp2(ERA40lon,ERA40lat,meanT,Lon,Lat);
    gmr = -0.03417;
    lr = [-6.1517E-03 -3.1831E-06 -1.5014E-07 1.8097E-09 1.1791E-10 ...
        -6.5359E-14 -9.5209E-15];
    dtdz = lr(1) + lr(2).*Lat + lr(3).*Lat.^2 ...
        + lr(4).*Lat.^3 + lr(5).*Lat.^4 + lr(6).* Lat.^5 ...
        + lr(7).*Lat.^6;
    dtdz = -dtdz;
    P=site_slp .* exp( (gmr./dtdz) .* ( log(site_T) - log(site_T - (Alt.*dtdz)) ) );
    [StFact,StFact_sp,StFact_mu] = StoneFactorL(Rc,P,site_slp);
elseif Atm==1;
    gmr = -0.03417;
    dtdz = 0.0065;
    SLP = 1013.25;
    P = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (Alt.*dtdz)) ) );
    [StFact,StFact_sp,StFact_mu] = StoneFactorL(Rc,P,SLP);
end
end

