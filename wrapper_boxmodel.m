function [heat_content,salt_content] = wrapper_boxmodel(X,Parameters)
% Wrapper function for running the boxmodel in UQLab
% X contains all parameters we want to explore (M=8)
% Parameters contain all parameters "in common" that we use for all model
% runs

[p,a] = get_model_default_parameters(); % default params, standard initialisation

%% Getting the parameters to be explored into variables that we can more easily recall
p.L         = X(1);
p.W         = X(2);
p.silldepth = X(3);
p.zgl       = X(4);

dTanom = X(5);
dSanom = X(6);
dQanom = X(7);
dDanom = X(8);

%% Set up model forcings
f.Ts  = Parameters.Tocn + Parameters.Tdecay .* dTanom;
f.Ss  = Parameters.Socn + Parameters.Sdecay .* dSanom;
f.Qsg = Parameters.Qglc + dQanom;
f.D   = Parameters.Dglc + dDanom;

fjord_run.f.zs = Parameters.zs;

% zs needs to be negative
a.H0 = double(get_fjord_boxes_from_density(f.Ts(:,1),f.Ss(:,1),-f.zs,p));
[a.T0, a.S0] = bin_ocean_profiles(f.Ts(:,1),f.Ss(:,1),-f.zs,a.H0,p); 

p.Snudge = get_interface_salinities(-f.zs,mean(f.Ts,2),mean(f.Ss,2),p);

fjord_run.p = p;
fjord_run.t = t;
fjord_run.f = f;
fjord_run.a = a;

%% Run the model itself, postprocess, and return the desired metrics for evaluation

[fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);
fjord_run.o               = postprocess_boxmodel(fjord_run);

% gets output quantities in "per unit volume"
heat_content = sum(fjord_run.o.hc)./(fjord.p.L.*fjord.p.W.*sum(fjord.p.H));
salt_content = sum(fjord_run.o.sc)./(fjord.p.L.*fjord.p.W.*sum(fjord.p.H));

end