function [ Max,Moyenne,Mediane,Err1SigInf,Err1SigSup,Err2SigInf,Err2SigSup ] = ParamDistribV11( Xabs,Ypdf ) 
%--------Description------------------------------------------------------
% This fuctions computes a set of parameters associated to a distribution
% Ypdf=f(Xabs).
%
% Input : 
%          Xabs       : absissa values of the Ypdf=f(Xabs) distribution (vector)
%          Ypdf       : ordinate values of the Ypdf=f(Xabs) distribution (vector)
%
% Output : 
%          Max        : Distribution mode
%          Moyenne    : Distribution mean
%          Mediane    : Distribution mediane
%          Err1SigInf : Associated 1 sigma inferior uncertainty
%          Err1SigSup : Associated 1 sigma superior uncertainty
%          Err2SigInf : Associated 2 sigma inferior uncertainty
%          Err2SigSup : Associated 2 sigma superior uncertainty
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

NbX=length(Xabs);
NbY=length(Ypdf);

% Check vector length
if NbX~=NbY;
    Max=NaN;
    Moyenne=NaN;
    Mediane=NaN;
    Err1SigInf=NaN;
    Err1SigSup=NaN;
    Err2SigInf=NaN;
    Err2SigSup=NaN;
    return
end

% Calculate mode
[~, IndiceMax]=max(Ypdf);
Max=Xabs(IndiceMax);


% Calculate mean
AireTrapz=trapz(Xabs,Ypdf);
Moyenne=trapz(Xabs,(Xabs.*Ypdf))/AireTrapz;

% Calculate other parameters

% Store elementary areas
NbX=length(Xabs);
VectAires=zeros(1,NbX);
for k=1:(NbX-1);
    DeltaX=Xabs(k+1)-Xabs(k);
    Trapeze=0.5*(Ypdf(k+1)+Ypdf(k))*DeltaX;
    VectAires(1,k+1)=Trapeze;
end

% Calculate total area
AireTot=sum(VectAires);

% Calculate cumulated area
VectAiresCum=zeros(1,NbX);
for k=1:NbX-1;
    VectAiresCum(1,k+1)=VectAiresCum(1,k)+VectAires(1,k+1);    
end

% Normalization
VectAiresCumNorm=VectAiresCum/AireTot;

% Calculate the distribution parameters
Bounds=[0.025 0.16 0.5 0.84 0.975];

% Mediane
tmp=abs(VectAiresCumNorm-Bounds(3));
[~,idx]=min(tmp);
Mediane=Xabs(idx);

% 2-Sigmas inferior error
tmp=abs(VectAiresCumNorm-Bounds(1));
[~,idx]=min(tmp);
Err2SigInf=Mediane-Xabs(idx);

% 1-Sigmas inferior error
tmp=abs(VectAiresCumNorm-Bounds(2));
[~,idx]=min(tmp);
Err1SigInf=Mediane-Xabs(idx);

% 1-Sigmas superior error
tmp=abs(VectAiresCumNorm-Bounds(4));
[~,idx]=min(tmp);
Err1SigSup=Xabs(idx)-Mediane;

% 2-Sigmas superior error
tmp=abs(VectAiresCumNorm-Bounds(5));
[~,idx]=min(tmp);
Err2SigSup=Xabs(idx)-Mediane;

end

