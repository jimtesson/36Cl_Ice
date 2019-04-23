%function [out] = inv_MCMC(Const_cosmo, Param_cosmo,Param_site,sf,Data,Z,param_mcmc,flag)
function out = inv_MCMC(Const_cosmo, Param_cosmo, Param_site, Sf, dataset, Z, param_mcmc, flag)
%INV_MCMC Summary of this function goes here
%  Detailed explanation goes here

%misfit = @(N_36Cl_model) sum((((N36Cl.meas-N_36Cl_model)./N36Cl.meas_uncert).^2).^.5); % misfit function



pp = @(t_expo,erosion) inv_misfit(Const_cosmo,Param_cosmo,Param_site,Sf,dataset,Z,erosion,t_expo,flag); % Function to obtain the modeled 36Cl as function of the erosion rate and the exposure age 
p = @(x) 1/(2*pi)*exp(-.5*x^2);% likelyhood function

% search parameters for the erosion rate:
e_min = param_mcmc.search_bd(1,1); % min bound value (mm/kyr)
e_max = param_mcmc.search_bd(1,2); % max bound value (mm/kyr)
b = param_mcmc.search_std(1); % std on the gaussian proposal function (mm/kyr)

% search parameters for the time exposure:
t_min = param_mcmc.search_bd(2,1); % min bound value (yr)
t_max = param_mcmc.search_bd(2,2); % max bound value (yr)
a = param_mcmc.search_std(2); % std on the gaussian proposal function (yr)

n_it = param_mcmc.n_it;

% pre allocate memory
t_out = zeros(n_it,1);
e_out = zeros(n_it,1);
rmsw_out = zeros(n_it,1);
rej_flag = zeros(n_it,1);
t_rej = zeros(size(t_out));
e_rej = zeros(size(e_out));


% starting point
t_out(1,:) = rand*(t_max-t_min) + t_min; % exposure age
e_out(1,:) = rand*(e_max-e_min) + e_min; % erosion rate

 for i = 2:n_it
     if(rand>=0.5)
 	t_in = (rand-0.5)*a + t_out(i-1,:); % sample from q ~ N(x(n-1),b^2)
    e_in = e_out(i-1,:);
     else
    t_in = t_out(i-1,:);
 	e_in = (rand-0.5)*b + e_out(i-1,:); % sample from q ~ N(x(n-1),b^2)
     end
     
     if((t_in>=t_min)&&(t_in<=t_max)&&(e_in>=e_min)&&(e_in<=e_max))
          %RMSw from forward model
          RMSw_in = pp(t_in,e_in);
          RMSw_out = pp(t_out(i-1,:),e_out(i-1,:));
          % likelyhood
          L_in = p(RMSw_in);
          L_out = p(RMSw_out);
        
    	alpha = min(L_in/L_out,1);
        u = rand; % sample u from U[0,1)
        if(u<=alpha) % accept
            t_out(i,:) = t_in;
            e_out(i,:) = e_in;
            rmsw_out(i,:) = RMSw_in;
            rej_flag(i) = 1;
        else % reject
            t_out(i,:) = t_out(i-1,:);
            e_out(i,:) = e_out(i-1,:);
            e_rej(i,:) = e_in;
            t_rej(i,:) = t_in;
        end      
     else % if the proposed model is out of the bounds, it is directly rejected
         L_in = 0;
         t_out(i,:) = t_out(i-1,:);
         e_out(i,:) = e_out(i-1,:);
         e_rej(i,:) = e_in;
         t_rej(i,:) = t_in;
     end
        
    

 end
 
%% Plot 
burnin = 0.5; % burning period
it_min = round(burnin * length(t_out(:,1)));
it_max = length(t_out(:,1));

m(1,1,:)=t_out(it_min:it_max);
m(2,1,:)=e_out(it_min:it_max);

frac = 100; % fraction of the models saved
m_thin = m(:,:,1:frac:end);

figure(3);
  plot(t_out(it_min:it_max),e_out(it_min:it_max),'.')  

  
figure(4)
    [C,lags,ESS]=eacorr(m);
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
  
figure(5)
    ecornerplot(m,'names',{'expo (yr)' 'Erosion (mm/kyr)'})

figure(6)
    ecornerplot(m_thin,'names',{'expo (yr)' 'Erosion (mm/kyr)'})
 out = 1;

%figure(7)
%histogram(rmsw_out(it_min:it_max));
end

function RMSw = inv_misfit(Const_cosmo,Param_cosmo,Param_site,Sf,dataset,Z,erosion,t_expo,flag)
[N_depth_profile, N_depth_profile_uncert] = depth_profile_speed(Const_cosmo,Param_cosmo,Param_site,Sf,erosion,t_expo,Z,flag);
RMSw = sum((((dataset(1,:)-N_depth_profile).^2).^.5)./dataset(2,:));
end

