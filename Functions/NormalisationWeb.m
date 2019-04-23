function [ OutMat ] = NormalisationWeb( Scheme,Nucl,Agev,Atm,U,Paleomag)
% This function produce a local production rate from a group of sample with
% the Lal-Stone model or the LSD model. It normalise the nuclide concentration to a
% normalisation point.

% Load parameters
% Nucl: 3He ou 10Be
% Atm: Atmosphere model 0: ERA40, 1: Atmosphere Standard

Age=Agev(1);
ErrAge=Agev(2);
Lat=U(:,1);
Lon=U(:,2);
Alt=U(:,3);
Ep=U(:,4);
Masque=U(:,5);
Dens=U(:,6);
Ero=U(:,7);
Conc=U(:,8);
ErrConc=U(:,9);

NbSpl=length(Lat);

% longitudes Correction
for i=1:NbSpl;
    if Lon(i)>180;
        Lon(i)=Lon(i)-360;
    end
end

% Choice of atmosphere type  
if Atm==0; % ERA 40
    load ERA40
else
    % Standard atmosphere
    meanP=[];
    meanT=[];
    ERA40lat=[];
    ERA40lon=[];
end

% Normalisation point
LatNorm=mean(Lat);
LonNorm=mean(Lon);
AltNorm=mean(Alt);

% Computation of the normalization factor

% Choice of the scaling scheme

if Scheme == 1
    % Stone-time dependent
Fact_Norm=StoneFactCT(Age/1000,LatNorm,LonNorm,AltNorm,Paleomag,Atm,ERA40lat,ERA40lon,meanP,meanT);
    
else
    
    %LSD
[VecT,VecSF]=LSDv9(LatNorm,LonNorm,AltNorm,Atm,Age,-1,Nucl,Paleomag);
Fact_Norm=LSDtimeInteg(Nucl,VecT,VecSF);

end

% Thickness and shielding correction
Attlgth=160;
CorrEp=ones(NbSpl,1);
for i=1:NbSpl;
    if Ep(i)>0;
        Densit=Dens(i);
        CorrEp(i)=Attlgth/(Densit*Ep(i))*(1-exp(-Densit*Ep(i)/Attlgth));
    end
end
ConcCorr=Conc./(CorrEp.*Masque);
ErrConcCorr=ErrConc./(CorrEp.*Masque);


% Spacio-temporal correction for each sample

Fact=zeros(NbSpl,1);

if Scheme == 1
    % Normalization done with the Stone factor

for i=1:NbSpl;
    Fact(i)=StoneFactCT(Age/1000,Lat(i),Lon(i),Alt(i),Paleomag,Atm,ERA40lat,ERA40lon,meanP,meanT);
end
else
    % Normalization done with the LSD factor
for i=1:NbSpl;
    [VecTe,VecSFe]=LSDv9(Lat(i),Lon(i),Alt(i),Atm,Age,-1,Nucl,Paleomag);
    Fact(i)= LSDtimeInteg(Nucl,VecTe,VecSFe);
end
end


% Normalization at depth = 0 and at average latitude longitude and altitude
ConcNorm=(ConcCorr*Fact_Norm)./Fact;

ErrConcNorm=(ErrConcCorr*Fact_Norm)./Fact;


% W=1./(ErrConcNorm.^2);
% ConcNormLSD=(ConcCorr*LSDNorm)./LSDFact;
% ErrConcNormLSD=(ErrConcCorr*LSDNorm)./LSDFact;
% WLSD=1./(ErrConcNormLSD.^2);


% Outliers removal using Peirce criterion (weighted standard deviation 
% is used as the reference). Only if more than 2 samples have been measured.  

if length(ConcNorm) >2;

[Concf,Errf]=Peirce(ConcNorm, ErrConcNorm);

else
    Concf=ConcNorm;
    Errf=ErrConcNorm;
end


% Weighted mean
[~,Wmean,ErrWm]=CorrWM([Concf Errf]);


% Output
OutMat=zeros(1,8);
OutMat(1)=Age;
OutMat(2)=ErrAge;
OutMat(3)=Wmean;
OutMat(4)=ErrWm;
OutMat(5)=LatNorm;
OutMat(6)=LonNorm;
OutMat(7)=AltNorm;
OutMat(8)=mean(Ero);
% OutMat(9)=WmeanLSD;
% OutMat(10)=ErrWmLSD;

end