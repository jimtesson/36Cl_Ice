function [WeightedStd] = Wstd (Xin, Win)

% Check input length
Nbr=length(Xin);
Nbr2=length(Win);
if Nbr~=Nbr2;
fprintf('Different length in Wstd')
return
end

% Calculate weighted standard deviation
Mean=sum(Win.*Xin)/sum(Win);
Numer=sum(((Xin-Mean).^2).*Win);
Denom=(Nbr-1)*sum(Win)/Nbr;
WeightedStd=sqrt(Numer/Denom);

end
