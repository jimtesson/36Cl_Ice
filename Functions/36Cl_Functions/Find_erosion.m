function best_erosion = Find_erosion(Const_cosmo, Param_cosmo, Param_site, Sf, dataset, Z, ero_guess, t_expo, flag)
% Function to obtain the best age of a single sample and a fixed erosion
% rate. It uses the fminsearchbnd function to find the age associated with
% the lowest misfit (here RMSw).
% INPUT :   
%           Const_cosmo : Constants for cosmogenic's calculation
%           Param_cosmo : Sample specific cosmogenic's variable previously
%                           calculated
%           Param_site : Site specific parameters.
%           Sf : Scaling factors.
%           dataset : vector containing the data [36Cl, 36Cl uncertainties, depth of the sample] .
%           Z : depth of the sample
%           ero_guess : guess erosion-rate in m/Myr, = starting value of the search.
%           t_expo : user-provided exposure age.
%           flag : user-provided flag to specify the model used.
%
%                   flag.model = 'exp' for the Schimmelpfennig ?2009 model
%                   (attenuation length approximated by exponential)
%
%                   flag.model = 'num' for the numerical approach (Marrero
%                   et al. 2018, attenuation length are calculated from the
%                   energy flux))
%
%                   flag.scaling_model = 'st' for the Stone 1999 scaling 
%                                               scheme.
%
%                   flag.scaling_model = 'sa' for the Lifton-Sato scaling 
%                                               scheme.
%                   flag.search = 'fminsearch' for the 'fminsearch' method. 
%                   flag.search = 'nmsmax' for the Nelder-Mead simplex
%                                   method.
%                                   Required: 
%                                       - flag.min_bounds: minimum bound
%                                       - flag.max_bounds: maximun bound
%                                 
%                   flag.search = 'fminbnd' for the 'fminsearch' method
%                                   using boundaries.
%                                               
%
%
%           
% Output : 
%            best_erosion: output erosion rate for the sample
%
% version 01/08/2018, written by TESSON J.

% Function to obtain the modeled 36Cl as function of the erosion rate and the exposure age
f = @(t_expo,erosion) inv_misfit(Const_cosmo,Param_cosmo,Param_site,Sf,dataset,Z,erosion,t_expo,flag); 
guess = ero_guess; % 
% Search
if(strcmp(flag.search,'fminsearch')==1)
    best_erosion = fminsearch(@(x) f(t_expo,x),guess);
elseif(strcmp(flag.search,'nmsmax')==1)
    precx = 10^-4;%'Simplex' accuracy
    [best_erosion,Obj3,N3] = nmsmax(@(x) -f(t_expo,x),guess,[precx inf inf 0 0]);
elseif(strcmp(flag.search,'fminbnd')==1)
    min = flag.min_bounds; % minimum bound for the search
    max = flag.max_bounds; % maximum bound for the search
    best_erosion = fminbnd(@(x) f(t_expo,x),min,max);
else
    stop
end


end

function RMSw = inv_misfit(Const_cosmo,Param_cosmo,Param_site,Sf,dataset,Z,erosion,t_expo,flag)
% Function to calculate the sample concentration and return the RMSw.
[N_depth_profile, N_depth_profile_uncert] = depth_profile_speed(Const_cosmo,Param_cosmo,Param_site,Sf,erosion,t_expo,Z,flag);
RMSw = sum((((dataset(1,:)-N_depth_profile).^2).^.5)./dataset(2,:));
end

