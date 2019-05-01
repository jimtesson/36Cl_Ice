function [ Data ] = Load_data_36(file_name)
%Function in order to get data from the xsl file

% import data
[~,~,raw] = xlsread(file_name);

% number of Samples
ns = length(raw(:,1))-1;
nrowmax = ns+1;

% samples names
Data.Samples = raw(2:nrowmax,1);

% Location
Data.Lat = cell2mat(raw(2:nrowmax,2));
Data.Lon = cell2mat(raw(2:nrowmax,3));
Data.Alt = cell2mat(raw(2:nrowmax,4));

% Depth. Must be = .0 if taken at surface
Data.Z =  cell2mat(raw(2:nrowmax,5));

% [36Cl]
Data.NuclCon = cell2mat(raw(2:nrowmax,6));
Data.NuclErr = cell2mat(raw(2:nrowmax,7));

% Shielding
Data.Shield = cell2mat(raw(2:nrowmax,8));

% density
Data.Dens = cell2mat(raw(2:nrowmax,9));

% Thickness
Data.Thick = cell2mat(raw(2:nrowmax,10));

% Erosion rate (cm/yr)
Data.Eros = cell2mat(raw(2:nrowmax,11));


% Age of formation
Data.Age_form = cell2mat(raw(2:nrowmax,12));
Data.Age_form_er = cell2mat(raw(2:nrowmax,13));

%% Target Chemistry
Data.Cl_targ = cell2mat(raw(2:nrowmax,14)); % (ppm)
Data.Cl_targ_er = cell2mat(raw(2:nrowmax,15)); % (ppm)
Data.SiO2_targ = cell2mat(raw(2:nrowmax,16)); % (in %)
Data.Al2O3_targ = cell2mat(raw(2:nrowmax,17)); % (in %)
Data.Fe2O3_targ = cell2mat(raw(2:nrowmax,18)); % (in %)
Data.Fe2O3_targ_er = cell2mat(raw(2:nrowmax,19)); % (in %)
Data.MnO_targ = cell2mat(raw(2:nrowmax,20)); % (in %)
Data.MgO_targ = cell2mat(raw(2:nrowmax,21)); % (in %)
Data.CaO_targ = cell2mat(raw(2:nrowmax,22)); % (in %)
Data.CaO_targ_er = cell2mat(raw(2:nrowmax,23)); % (in %)
Data.Na2O_targ = cell2mat(raw(2:nrowmax,24)); % (in %)
Data.K2O_targ = cell2mat(raw(2:nrowmax,25)); % (in %)
Data.K2O_targ_er = cell2mat(raw(2:nrowmax,26)); % (in %)
Data.TiO2_targ = cell2mat(raw(2:nrowmax,27)); % (in %)
Data.TiO2_targ_er = cell2mat(raw(2:nrowmax,28)); % (in %)
Data.P2O5_targ = cell2mat(raw(2:nrowmax,29)); % (in %)
Data.LOI_targ = cell2mat(raw(2:nrowmax,30)); % (in %)
Data.H2O_targ = cell2mat(raw(2:nrowmax,31)); % (in %)
Data.CO2_targ = cell2mat(raw(2:nrowmax,32)); % (in %)


%% Bulk Chemistry
Data.SiO2_bulk = cell2mat(raw(2:nrowmax,33)); % (in %)
Data.Al2O3_bulk = cell2mat(raw(2:nrowmax,34)); % (in %)
Data.Fe2O3_bulk = cell2mat(raw(2:nrowmax,35)); % (in %)
Data.MnO_bulk = cell2mat(raw(2:nrowmax,36)); % (in %)
Data.MgO_bulk = cell2mat(raw(2:nrowmax,37)); % (in %)
Data.CaO_bulk = cell2mat(raw(2:nrowmax,38)); % (in %)
Data.NaO_bulk = cell2mat(raw(2:nrowmax,39)); % (in %)
Data.K2O_bulk = cell2mat(raw(2:nrowmax,40)); % (in %)
Data.TiO2_bulk = cell2mat(raw(2:nrowmax,41)); % (in %)
Data.P2O5_bulk = cell2mat(raw(2:nrowmax,42)); % (in %)
Data.H2O_bulk = cell2mat(raw(2:nrowmax,43)); % (in %)
Data.CO2_bulk = cell2mat(raw(2:nrowmax,44)); % (in %)
Data.S_bulk = cell2mat(raw(2:nrowmax,45)); % (in %)
Data.LOI_bulk = cell2mat(raw(2:nrowmax,46)); % (in %)
Data.Li_bulk = cell2mat(raw(2:nrowmax,47)); % (in %)
Data.B_bulk = cell2mat(raw(2:nrowmax,48)); % (in %)
Data.Cl_bulk = cell2mat(raw(2:nrowmax,49)); % (in %)
Data.Cr_bulk = cell2mat(raw(2:nrowmax,50)); % (in %)
Data.Co_bulk = cell2mat(raw(2:nrowmax,51)); % (in %)
Data.Sm_bulk = cell2mat(raw(2:nrowmax,52)); % (in %)
Data.Gd_bulk = cell2mat(raw(2:nrowmax,53)); % (in %)
Data.Th_bulk = cell2mat(raw(2:nrowmax,54)); % (in %)
Data.Th_bulk_er = cell2mat(raw(2:nrowmax,55)); % (in %)
Data.U_bulk = cell2mat(raw(2:nrowmax,56)); % (in %)
Data.U_bulk_er = cell2mat(raw(2:nrowmax,57)); % (in %)


end

