
function out=getcurrentsf(sf,t,flag)
% Function to get the elevation/latitude time dependent scaling factors for
% neutrons by nuclide (SelXX) and thermal (SFth) and epithermal neutrons (SFeth) 
% by interpolating the tables of time dependent scaling factors. 
% INPUT:
%       sf : structure containing the scaling factors.
%       t : time vector (yr) for which we want the scaling factors.
%       flag.scaling_model : scaling scheme for neutrons ('st' for stone 2000,
%       'sa' for Lifton-Sato).
%       nuclide : for LSD scheme, specify the nuclide: 'cl' for chlorine.
%
% OUTPUT:
%       out: struct containing the scaling factors 
%
% version 01/08/2019 TESSON J.
% Modified after Marrero 2016, Chronus code.
%

if (max(t)>max(sf.tv))
  t(length(t))=max(sf.tv);
  warning('extremely old time in getcurrentsf!');
end

%
% Look for times after 2110 AD.  This is most likely an error.
% 
if (t(1)<-300)
  t
  warning('extremely young time in getcurrentsf!');
end
%
% For all times after t=0 (2010), use the 2010 scaling factors by
% forcing t=0.  We don't let this get too far out of hand- this is
% limited to 2110 (see above.)
%
if (t(1) < 0)
  t(t<0)=0.0;
end


switch flag.scaling_model
case 'st' % for the Lal/Stone scheme.  
        SF = [sf.SF_St' sf.SF_St_sp' sf.SF_St_mu'];
        SF_interp = interp1(sf.SF_Vectime,SF,t,'linear')';
        
        SF_St_int = SF_interp(1,:);
        SF_St_sp_int = SF_interp(2,:);
        SF_St_mu_int = SF_interp(3,:);
        
        out.SelSF = SF_St_int;
        out.Sel36Ca = SF_St_sp_int;
        out.Sel36K = SF_St_sp_int;
        out.Sel36Ti = SF_St_sp_int;
        out.Sel36Fe = SF_St_sp_int;
        out.SFth = SF_St_sp_int;
        out.SFeth = SF_St_sp_int;
        out.SFmu = SF_St_mu_int;
        
case 'sa' % For LSD scheme
        SF = [sf.SF_Sa36Ca' sf.SF_Sa36K' sf.SF_Sa36Ti' sf.SF_Sa36Fe' sf.SF_Sath' sf.SF_Saeth' sf.SF_Sf'];
        SF_interp = interp1(sf.tv,SF,t,'linear')';
        
        out.Sel36Ca = SF_interp(1,:);
        out.Sel36K = SF_interp(2,:);
        out.Sel36Ti = SF_interp(3,:);
        out.Sel36Fe = SF_interp(4,:);
        out.SFth = SF_interp(5,:);
        out.SFeth = SF_interp(6,:);
        out.SelSF = SF_interp(7,:);
        
        if(strcmp(flag.muon,'exp')) % special case where the muon model is exponential, use of lal-stone scaling
            SF = sf.SF_St_mu';
            SF_interp = interp1(sf.SF_Vectime,SF,t,'linear')';
            out.SFmu =  SF_interp(1,:);
        end
    
end
