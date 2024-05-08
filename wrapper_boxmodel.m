function [fjord_out] = wrapper_boxmodel(X,Parameters,flag_debug)
% function [ohc_mean,osc_mean] = wrapper_boxmodel(X,Parameters,flag_debug)
% function [output_time,heat_content,salt_content,status,fjord_run] = wrapper_boxmodel(X,Parameters)
% Wrapper function for running the boxmodel in UQLab
% X contains all parameters we want to explore (10)
% Parameters contain all parameters "in common" that we use for all model
% runs
[p,a] = get_model_default_parameters(); % default params, standard initialisation
p.H = Parameters.H;
p.Hmin=5;
% p.M0=0;      % no icebergs if zero
p.C0=X(9);  % no shelf-exchange if zero
p.P0=X(8);  % no meltwater plume if zero

p.dt = Parameters.t(2)-Parameters.t(1);
%% Getting the parameters to be explored into variables that we can more easily recall
p.L         = X(1);
% p.W         = X(2);
p.H         = X(2);
% TODO: add H and bump up the indices of the rest
p.W         = X(1) ./ X(3); % W  = L /(L/W)
% p.W         = X(1)./X(2); % W  = L /(L/W)
p.silldepth = -X(4) .* X(2); % Zs = (Zs/H) * H
p.zgl       = -X(5) .* X(2); % Zg = (Zg/H) * H 
% p.silldepth = X(3); 
% p.zgl       = X(4); 

dQamp  = X(6);
dDanom = X(7); % ignoring solid-ice discharge (for now?)

if length(X) > 10
    dTanom = X(10);
    dSanom = X(11);
    % Xfrq   = X(12);
else
    i_cast=X(10);
end


%% Set up model forcings
% Xamp=0.92;  % if using anomaly to profiles
% Xamp=80; % if using isopycnal stretching; maximum range in isopycnal depth STD from Jackson & Straneo (2016; JPO)

% we modulate the "storm events" to only happen in winter
% doy_peak = 173; % june 22nd
% winter_wave = -sin(2*pi/365 .* (Parameters.t - (doy_peak - 365/4))); % 1 at peak winter, -1 at peak summer
% winter_wave(winter_wave < 0) = 0;
% Xper = (Xamp .* sin(-2*pi*Xfrq*Parameters.t)) .* winter_wave; 

% Ocean forcings
if length(X) > 10
    % f.Ts  = (Parameters.Tocn + (Parameters.Tdec .* (dTanom + dTanom .* Xper)))'; % if using anomaly to profiles
    % f.Ss  = (Parameters.Socn + (Parameters.Sdec .* (dSanom + dSanom .* Xper)))'; % if using anomaly to profiles
    f.Ts  = (Parameters.Tocn + (Parameters.Tdec .* dTanom))'; % if using isopycnal stretching
    f.Ss  = (Parameters.Socn + (Parameters.Sdec .* dSanom))'; % if using isopycnal stretching
    f.zs  = -Parameters.zs;
    
    % if ismember(Parameters.regID,[2,4,6]) % isopycnal stretching due to storms only applies to eastern Greenland    
    %     [f.Ts,f.Ss] = heave_profiles(f.Ts,f.Ss,f.zs,Xper,p.sigma_bnds(2)); % using Cowton et al. value for the pycnocline
    % end
    zs = flip(f.zs);
    Ts = flip(f.Ts,1); 
    Ss = flip(f.Ss,1);
else
    temp_profile  = Parameters.ocn(i_cast).T;
    salt_profile  = Parameters.ocn(i_cast).S;
    depth_profile = Parameters.ocn(i_cast).depth;

    f.Ts = repmat(temp_profile,length(Parameters.t),1);
    f.Ss = repmat(salt_profile,length(Parameters.t),1);
    f.zs = repmat(depth_profile,length(Parameters.t),1);
    Ts = f.Ts;
    Ss = f.Ss;
    zs = f.zs;
end



a.H0 = double(get_fjord_boxes_from_density(mean(Ts,2),mean(Ss,2),zs,p));
% Check sum of layer thicknesses is equal to fjord depth.
if abs(sum(a.H0)-p.H) > 1e-10
    disp('Error: box thicknesses must sum to fjord depth');
    return
end
% Alternative: evenly distribute layer thicknesses above the sill
% a.H0 = ones([1,p.N]).*(abs(p.silldepth.*p.sill)/p.N);
% if p.sill, a.H0 = [a.H0,p.H-abs(p.silldepth)]; end
[a.T0, a.S0] = bin_ocean_profiles(Ts(:,1),Ss(:,1),zs,a.H0,p); 

p.Snudge = get_interface_salinities(zs,mean(Ts,2),mean(Ss,2),p);

% Glacier forcings
f.Qsg = Parameters.Qglc .* dQamp;
f.D   = Parameters.Dglc + dDanom;
f.zi = Parameters.zi;
f.xi = Parameters.xi;
a.I0 = Parameters.I0;

%% For plotting purposes
if nargin > 2
    p.plot_runtime=flag_debug; % for debugging purposes
end
% run figure_forcings_summary.m

%% Preparing fjord encapsulating structure
fjord_run.t = Parameters.t;
fjord_run.p = p;
fjord_run.f = f;
fjord_run.a = a;

%% Run the model itself, postprocess, and return the desired metrics for evaluation

fjord_out = struct("time",[],"temp",[],"salt",[],"H",[],"ts",[],"ss",[],"zs",[],"p",[],"phi",[],"qvs",[],"qsg",[],"qts",[]);
try
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);

    fjord_out.temp  = fjord_run.s.T;
    fjord_out.salt  = fjord_run.s.S;
    fjord_out.H     = fjord_run.s.H;
    fjord_out.qvs   = fjord_run.s.QVs;
    fjord_out.qsg   = fjord_run.s.QSg;
    fjord_out.phi   = fjord_run.s.phi;
    fjord_out.qts   = fjord_run.s.QTs;
    fjord_out.ts    = fjord_run.s.Ts;
    fjord_out.ss    = fjord_run.s.Ss;
catch ME
    fjord_run.s.T  = NaN;
    fjord_run.s.S  = NaN;
    fjord_run.s.H  = NaN;
    fjord_out.qvs  = NaN;
    fjord_out.qsg  = NaN;
    fjord_out.qts  = NaN;
    fjord_out.phi  = NaN;
    fjord_run.s.Ts = NaN;
    fjord_run.s.Ss = NaN;
    
    fprintf('run failed: %s\n',ME.message)
end

% fjord_run.o = postprocess_boxmodel(fjord_run);
% heat_content = sum(fjord_run.o.hc)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
% salt_content = sum(fjord_run.o.sc)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
fjord_out.time  = fjord_run.t;
fjord_out.zs    = fjord_run.f.zs;
fjord_out.p     = fjord_run.p;

end