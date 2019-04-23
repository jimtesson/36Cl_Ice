function [ ValeurDiscret ] = TrouvValeurDiscret( ValeurBrute,PasDiscret )
%--------Description------------------------------------------------------
% Find the closest value of ValeurBrute on a grid using PasDiscret as a
% discretisation time step. 
%
% Input : 
%          ValeurBrute   : Initial value
%          PasDiscret    : Time step
% Output : 
%          ValeurDiscret : Closest value considering the time step
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

reste=mod(ValeurBrute,PasDiscret);
if reste==0;
    ValeurDiscret=ValeurBrute;
elseif reste<(PasDiscret/2);
    ValeurDiscret=ValeurBrute-reste;
else
    ValeurDiscret=ValeurBrute-reste+PasDiscret;
end

end

