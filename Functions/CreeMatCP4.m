function [ MatCP ] = CreeMatCP4( u,Pas )
%--------Description------------------------------------------------------
% This function prepares the Probability Density Function (PDF) output
% data. It recieves a cell array containing a define number of PDF (one
% time vector, one density vector) and replace all the PDF on a common time
% vector for an easier exportation of the PDF.
%
% Input :   
%        u     : Cell array containing PDF data
%        Pas   : Time step (ka)
%
% Output :   
%        MatCP : Matrix bearing all the PDF on the same time axe
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


% Elimination of empty vectors
NbVide=0;
for i=1:length(u);
    if isempty(u{i});
        NbVide=NbVide+1;
    end
end

if NbVide>0;
    k=length(u);
    NbSuppr=0;
    while NbSuppr<NbVide;
        if isempty(u{k});
            u(k)=[];
            NbSuppr=NbSuppr+1;
        end
        k=k-1;
    end
end


% Interpolation of time vectors
NbEch=length(u)/2;
for i=1:NbEch;
    % current vector 
    Vdates=u{2*i-1};
    Vpdf=u{2*i};
    % initial and final values
    TInf=ceil(Vdates(1)/Pas)*Pas;
    TSup=floor(Vdates(end)/Pas)*Pas;
    % Interpolation
    Vdatesi=TInf:Pas:TSup;
    Vpdfi=interp1(Vdates,Vpdf,Vdatesi,'spline');
    % Remplacement in u
    u{2*i-1}=Vdatesi;
    u{2*i}=Vpdfi;
end


% Reference time vector
% Minimum age
VectMin=zeros(1,NbEch);
for i=1:NbEch;
    VectMin(i)=u{2*i-1}(1);
end
DateMin=min(VectMin);

% Maximum age
VectMax=zeros(1,NbEch);
for i=1:NbEch;
    VectMax(i)=u{2*i-1}(end);
end
DateMax=max(VectMax);

% time vector 
VectTempsTot=DateMin:Pas:DateMax;
NbDatesTot=length(VectTempsTot);
VectTempsTot=VectTempsTot';


% Matrix creation
MatCP=zeros(length(VectTempsTot),1+NbEch);
MatCP(:,1)=VectTempsTot;

for i=1:NbEch;
    VectPDFPropre=zeros(NbDatesTot,1);
    VectTemps=u{2*i-1};
    VectPDF=u{2*i};
    IDeb=1;
    while VectTempsTot(IDeb)-VectTemps(1)<(-1e-5*Pas);
        IDeb=IDeb+1;
    end
    IFin=NbDatesTot;
    while (VectTempsTot(IFin)-VectTemps(end))>(1e-5*Pas);
        IFin=IFin-1;
    end
    VectPDFPropre(1:IDeb-1)=zeros(IDeb-1,1); 
    VectPDFPropre(IFin+1:NbDatesTot)=zeros(NbDatesTot-IFin,1);
    VectPDFPropre(IDeb:IFin)=VectPDF';
    MatCP(:,i+1)=VectPDFPropre;
end

end