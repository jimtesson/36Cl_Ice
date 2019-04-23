function out = LiftonSatoSX(h,Rc,SPhi,w,consts)

% Implements the Lifton Sato et al scaling scheme for spallation.
%
% Syntax: scalingfactor = LiftonSatoSX(h,Rc,SPhi,w,consts);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%   
%
% Vectorized. Send in scalars or vectors of common length. 
%
% Modified by Shasta Marrero (NMT) to include the reactions for Ti & Fe to
% produce chlorine-36. July 2011.
%
% Written by Nat Lifton 2011, Purdue University
% Based on code by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% convert pressure to atmospheric depth

X = h.*1.019716;

% Sref = 509.6; %1950-1850 mean
% Sref = 416.518; %Long-term mean
Sref = SPhi(1); %587.4 MV mean 2001-2010 from Usoskin et al 2011. USE THIS!
% Sref = 400; %Idealized solar minimum value, similar to long-term mean (11.4 ka)
% Changed from SPhi(1) 12/8/11 - more consistent with other calculations
% assuming idealized reference states such as SLHL pressure and Rc = 0. BUT
% does not normalize flux at t0 to 1 at SLHL - more like 0.92 - results in
% really high production predictions. Need to use the value from the 2001-2010.
% Go back to SPhi(1)... 12/15/11

Rcref = 0;
Href = 1013.25;
Xref = Href.*1.019716;


% Full version with cross sections and includes thermal and
% epithermal fluxes
nuclide = 36;
tmp_N = Neutrons_36(Href,Rcref,Sref,w,consts,nuclide);
[ethflux,thflux] = NeutronsLowE_36(Href,Rcref,Sref,w);
tmp_P = Protons_36(Href,Rcref,Sref,consts,nuclide);

SpRef = tmp_N.nflux + tmp_P.pflux;% Sato et al. (2008) Reference hadron flux integral >1 MeV
ClCaRef = tmp_N.P36Can + tmp_P.P36Cap;
ClKRef = tmp_N.P36Kn + tmp_P.P36Kp;
ClTiRef = tmp_N.P36Tin + tmp_P.P36Tip;
ClFeRef = tmp_N.P36Fen + tmp_P.P36Fep;
EthRef = ethflux;
ThRef = thflux;

% Full version with cross sections and includes thermal and
% epithermal fluxes
tmp_N = Neutrons_36(h,Rc,SPhi,w,consts,nuclide);

[ethflux,thflux] = NeutronsLowE_36(h,Rc,SPhi,w);

tmp_P = Protons_36(h,Rc,SPhi,consts,nuclide);

Site.sp = ((tmp_N.nflux + tmp_P.pflux))./SpRef;
Site.ClCa = (tmp_N.P36Can + tmp_P.P36Cap)./ClCaRef;
Site.ClK = (tmp_N.P36Kn + tmp_P.P36Kp)./ClKRef;
Site.ClTi = (tmp_N.P36Tin + tmp_P.P36Tip)./ClTiRef;
Site.ClFe = (tmp_N.P36Fen + tmp_P.P36Fep)./ClFeRef;
Site.eth = ethflux./EthRef;
Site.th = thflux./ThRef;

out = Site;

