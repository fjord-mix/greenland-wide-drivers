% function ohc_mean = wrapper_boxmodel(X,Parameters)
function [output_time,heat_content,salt_content,status,fjord_run] = wrapper_boxmodel(X,Parameters)
% Wrapper function for running the boxmodel in UQLab
% X contains all parameters we want to explore (M=8)
% Parameters contain all parameters "in common" that we use for all model
% runs

[p,a] = get_model_default_parameters(); % default params, standard initialisation
p.H = Parameters.H;
p.trelax=365/2;
p.Hmin=5;
p.M0=0;
% p.K0 = 0;
p.dt = Parameters.t(2)-Parameters.t(1);
%% Getting the parameters to be explored into variables that we can more easily recall
p.L         = X(1);
p.W         = X(2);
p.silldepth = X(3);
p.zgl       = X(4);

dTanom = X(5);
dSanom = X(6);
Xfrq   = X(7);
dQamp  = X(8);
dDanom = X(9);

%% Set up model forcings
Xamp=0.92;  % up to 92% increase (from Fraser & Inall, 2017): avg. of heat transport increase in the three inner-fjord sections

% we modulate the "storm events" to only happen in winter
doy_peak = 173; % june 22nd
winter_wave = -sin(2*pi/365 .* (Parameters.t - (doy_peak - 365/4))); % 1 at peak winter, -1 at peak summer
winter_wave(winter_wave < 0) = 0;
Xper = (Xamp .* sin(-Xfrq*Parameters.t)) .* winter_wave; 

% Ocean forcings
f.Ts  = (Parameters.Tocn + (Parameters.Tdec .* (dTanom + dTanom .* Xper)))';
f.Ss  = (Parameters.Socn + (Parameters.Sdec .* (dSanom + dSanom .* Xper)))';
% f.Ts  = (Parameters.Tocn + (Parameters.Tdec .* dTanom))';
% f.Ss  = (Parameters.Socn + (Parameters.Sdec .* dSanom))';
f.zs  = -Parameters.zs;

zs = flip(f.zs);
Ts = flip(f.Ts,1); 
Ss = flip(f.Ss,1);

a.H0 = double(get_fjord_boxes_from_density(mean(Ts,2),mean(Ss,2),zs,p));

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
p.plot_runtime=0;
% run figure_forcings_summary.m

%% Preparing fjord encapsulating structure
fjord_run.t = Parameters.t;
fjord_run.p = p;
fjord_run.f = f;
fjord_run.a = a;

%% Run the model itself, postprocess, and return the desired metrics for evaluation

[fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);
status = fjord_run.s.status;
fjord_run.o               = postprocess_boxmodel(fjord_run);

% gets output quantities in "per unit volume"
heat_content = sum(fjord_run.o.hc)./(fjord_run.p.L.*fjord_run.p.W.*sum(fjord_run.p.H));
salt_content = sum(fjord_run.o.sc)./(fjord_run.p.L.*fjord_run.p.W.*sum(fjord_run.p.H));
output_time  = fjord_run.s.t;
ohc_mean = mean(heat_content);
% average over the last N years of the run - to be implemented when runs
% are more stable
% if status
%     plot_fluxes(fjord_run)
%     plot_outputs(fjord_run)
% end
end