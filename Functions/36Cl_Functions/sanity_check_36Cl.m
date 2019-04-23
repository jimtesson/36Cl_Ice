function [ check36 ] = sanity_check_36Cl( Param_cosmo,Sf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N_samples = numel(Param_cosmo);
Mess = cell(3);
for i=1:N_samples % loop over samples
    

    % check proportion of Thermal and Epithermal production rate
    P35Cl = Sf{i}.SF_St_sp(1).*(Param_cosmo{i}.P_th+Param_cosmo{i}.P_eth)...
            ./Param_cosmo{i}.P(1).*100;

    if(P35Cl>45)
        Mess{i}=sprintf('Epithermal and Thermal 36Cl production rate is high (> 45%)');
    end

    % check if the sample has < 100 ppm of Chlorine
    if(Param_cosmo{i}.ppm_targ(61)>100)
        Mess{i}=sprintf('Sample has [Cl] > 100 ppm');
    end

    % check if Radiogenic production rate is not too high
    if(Param_cosmo{i}.P_rad>(0.45.*Param_cosmo{i}.P))
        Mess{i}=sprintf('Radiogenic 36Cl is high (>45%)');
    end
    
    % check if the sample is at equilibrium ?
end

check36 = Mess;