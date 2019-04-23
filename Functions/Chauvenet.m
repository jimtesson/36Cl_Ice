function [ Xout,Errout ] = Chauvenet( Xin,Errin,NbSig )
% Fonction qui applique le critre de chauvenet sur la srie Xin avec des
% carts types pondrs par Win pour un nombre de sigma NbSig

% Parameters
NbSpl=length(Xin);

if NbSpl==1;
    Xout=Xin;
    Errout=Errin;
    fprintf('Warning Chauvenet, one sample input\n')
    return
end

[Li,Co]=size(Xin);
Win=1./(Errin.^2);
Dist=NbSig*Wstd(Xin,Win);
Xout=zeros(Li,Co);
Errout=zeros(Li,Co);

% Mean
Mean=sum(Win.*Xin)/sum(Win);

% Selection
indice=1;
for i=1:NbSpl;
    if abs(Xin(i)-Mean)<=Dist;
        Xout(indice)=Xin(i);
        Errout(indice)=Errin(i);
        indice=indice+1;
    end
end

% Shortening
Xout(indice:end)=[];
Errout(indice:end)=[];
end

