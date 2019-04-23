function [ SFi ] = LSDtimeInteg( Iso,Time,SF )
%--------Description------------------------------------------------------
% This function computes a time integrated value from the scaling factor
% vectors produced by the LSD model.
% 
% Input : 
%           Iso  : Nuclide of interest (3He, 10Be, ...)
%           Time : time vector
%           SF   : Scaling factors vector
%
% Output :  
%           SFi  : Time integrated scaling factor
%
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

if (Iso~=3 && Iso~=10);
    fprintf('Isotope non reconnu')
    SFi=NaN;
    return
end

if Iso==10;
    tBe=1387000;
    LambdaBe=log(2)/tBe;
    SF=SF.*exp(-LambdaBe*Time);
end
SFi=trapz(Time,SF)/Time(end);

end

