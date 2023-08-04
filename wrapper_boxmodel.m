function [heat_content,salt_content,status] = wrapper_boxmodel(X,Parameters)
% Wrapper function for running the boxmodel in UQLab
% X contains all parameters we want to explore (M=8)
% Parameters contain all parameters "in common" that we use for all model
% runs

[p,a] = get_model_default_parameters(); % default params, standard initialisation
p.H = Parameters.H;

%% Getting the parameters to be explored into variables that we can more easily recall
p.L         = X(1);
p.W         = X(2);
p.silldepth = X(3);
p.zgl       = X(4);

dTanom = X(5);
dSanom = X(6);
Qamp   = X(7);
dDanom = X(8);

%% Set up model forcings

% Ocean forcings
f.Ts  = Parameters.Tocn' + dTanom;% .* Parameters.Tdecay;
f.Ss  = Parameters.Socn' + dSanom;% .* Parameters.Sdecay;
f.zs  = Parameters.zs;

zs = flip(-f.zs);
Ts = flip(f.Ts,1); 
Ss = flip(f.Ss,1);

a.H0 = double(get_fjord_boxes_from_density(mean(Ts,2),mean(Ss,2),zs,p));
[a.T0, a.S0] = bin_ocean_profiles(Ts(:,1),Ss(:,1),zs,a.H0,p); 
p.Snudge = get_interface_salinities(zs,mean(Ts,2),mean(Ss,2),p);

% Glacier forcings
f.Qsg = Parameters.Qglc .* Qamp;
f.D   = Parameters.Dglc + dDanom;
f.zi = Parameters.zi;
f.xi = Parameters.xi;
a.I0 = Parameters.I0;

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

% average over the last N years of the run - to be implemented when runs
% are more stable


end