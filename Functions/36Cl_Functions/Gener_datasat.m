% Generate dataset
clear Z
load data.mat

% parameters
n_samples = 10;
dz = 30; % Sample every dz cm
Z = zeros(1,n_samples);

% declaration
tmp = DATA_CREP;
clear DATA_CREP;
DATA_CREP = cell(n_samples,56);

% create the dataset
for j=1:n_samples
    for i=1:56
        DATA_CREP{j,i} = tmp{1,i};
    end
    Z(j) =  (j-1) .* dz ;
end

clear j i dz n_samples tmp_data