function [Age,Err,XPDF,PDF]=CalibGauss(AgeCosmo,ErrCosmo,VecTCal,VecYCal,AgeDiscr)
%--------Description------------------------------------------------------
% This function converts the gaussian Probability Density Function (PDF)
% associated with the non-corrected ages into an indeterminate corrected
% age PDF. The conversion needs the calibration curve given by vectors
% VecTCal and VecYCal.
%
% Input : 
%        AgeCosmo  : Preliminary age with no geomagnetic corrections (in ka)                
%        ErrCosmo  : 1-sigma asociated uncertainty (in ka)
%        VecTCal   : Time vector for non calibrated ages (in ka)
%        VecYCal   : Time vector for calibrated ages (in ka)
%        AgeDiscr  : Time step (ka)
% Output  : 
%        Age       : Past geomagnetic activity corrected ages (ka)
%        Err       : 1-sigma associated uncertainty (ka) (production rate uncertainty is considered)
%        XPDF      : Vector of time values for the associated PDF (ka)
%        PDF       : Vector of density function values for the associated PDF (production rate uncertainty is not considered)

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

% PDf boundaries
NbSigma=7;

AgeNSigmaInfBrut=AgeCosmo-NbSigma*ErrCosmo;
AgeNSigmaSupBrut=AgeCosmo+NbSigma*ErrCosmo;

if AgeNSigmaInfBrut<VecYCal(2);
    AgeNSigmaInfBrut=VecYCal(2);
end

if AgeNSigmaSupBrut>VecYCal(end-1);
    AgeNSigmaSupBrut=VecYCal(end-1);
end

% Set on the grid
AgeNSigmaInfDiscr=TrouvValeurDiscret(AgeNSigmaInfBrut,AgeDiscr);
AgeNSigmaSupDiscr=TrouvValeurDiscret(AgeNSigmaSupBrut,AgeDiscr);

% Create time vector for PDF
XPDFCosmo=AgeNSigmaInfDiscr:AgeDiscr:AgeNSigmaSupDiscr;

% Creat PDF
PDFCosmo=normpdf(XPDFCosmo,AgeCosmo,ErrCosmo);

% Stop if the age is too old
if trapz(XPDFCosmo,PDFCosmo)<0.975;
    AgeReporte=TrouvValeurDiscret(AgeCosmo,AgeDiscr);
    indice= fix(AgeReporte/AgeDiscr)+1;
    Age=VecTCal(indice);
    Err=ErrCosmo*Age/AgeCosmo;
    XPDF=[];
    PDF=[];
    return
end

% Find the begining and the end of the final vector
IndiceDebut=fix(XPDFCosmo(1)/AgeDiscr)+1;
IndiceFin=IndiceDebut+length(XPDFCosmo)-1;

% Find derivative of the calibration curve
DerYCal=zeros(1,length(XPDFCosmo));
for k=IndiceDebut:IndiceFin;
    DerYCal(1,k-IndiceDebut+1)=(VecYCal(k+1)-VecYCal(k-1))/(VecTCal(k+1)-VecTCal(k-1));
end

% Pick final time vector and PDF
XPDF=VecTCal(IndiceDebut:IndiceFin);
PDF=PDFCosmo.*DerYCal;

% Compute associated age and 1-sigma error 
[~,~,Age,Err1SigInf,Err1SigSup,~,~]=ParamDistribV11(XPDF,PDF);
Err=(Err1SigInf+Err1SigSup)/2;
end