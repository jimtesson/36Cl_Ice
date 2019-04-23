function [ ConcNorm ] = StdNorm( Conc, Norm, Nucl )
% Fonction for 10Be and 26Al concentration normalization versus official
% standard. Output matrix has the dimension of the concentration input.
% Conc is suppose to be a matrix with for each sample a concentration and
% an error. 

% Check if 3He
if Nucl==3;
    ConcNorm=Conc;
else
    % Names and corrections
    if Nucl==26;
        STD = {'KNSTD'; 'ZAL94'; 'ZAL94N'; 'SMAL11'; 'Z92-0222'};
        facteur = [1; 0.9134; 1; 1.021; 1];
    else
        STD = {'07KNSTD'; 'KNSTD'; 'NIST_Certified'; 'NIST_30000'; 'NIST_30200'; 'NIST30_300'; 'NIST_30600'; 'NIST_27900'; 'BEST433'; 'S555'; 'S2007'; 'BEST433N'; 'S555N'; 'S2007N'; 'LLNL31000'; 'LLNL10000'; 'LLNL3000'; 'LLNL1000'; 'LLNL300'};
        facteur = [1; 0.9042; 1.0425; 0.9313; 0.9251; 0.9221; 0.9130; 1; 0.9124; 0.9124; 0.9124; 1; 1; 1; 0.8761; 0.9042; 0.8644; 0.9313; 0.8562];
    end
    % Prepare output
    [NbL,NbC]=size(Conc);
    ConcNorm=zeros(NbL,NbC);
    % Calculate
    if length(Conc)>1;
        for i=1:length(Conc);
            index=strcmp(STD,Norm{i});
            ConcNorm(i)=Conc(i)*facteur(index);
        end
    else
        index=strcmp(STD,Norm);
        ConcNorm=facteur(index)*Conc;
    end
end
end