function [params,iOpts,probs] = define_model_param_distrib(datasets,fjords_compilation,i_reg)
%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: W, L, Zg, Zs, Ta, Sa, Da, Qa

%% Get compiled fjord data in pre-processed structure
% also selects the desired period in time

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2014,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step in days

fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

%% Gets distribution of fjord geometry parameters
% since the geometry does not change with time, we can stick to the
% original fjords structure
% verbose.plot=1;
fjord_stats = print_fjord_statistics(fjords_compilation);%,verbose);


iOpts.Marginals(1).Type = 'ks';
iOpts.Marginals(1).Inference.Data = fjord_stats.W.total;
iOpts.Marginals(1).Parameters     = fjord_stats.W.total;
probs(1) = fjord_stats.W.pd;

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = fjord_stats.L.total;
iOpts.Marginals(end).Parameters     = fjord_stats.L.total;
probs(end+1) = fjord_stats.L.pd;

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = fjord_stats.Zs.total;
iOpts.Marginals(end).Parameters     = fjord_stats.Zs.total;
probs(end+1) = fjord_stats.Zs.pd;

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = fjord_stats.Zg.total;
iOpts.Marginals(end).Parameters     = fjord_stats.Zg.total;
probs(end+1) = fjord_stats.Zg.pd;

%% Gets the ocean forcing

[tocn_clim, tocn_anom, depths] = get_var_clim_by_region(fjords_processed,'Ts');
[socn_clim, socn_anom, ~]      = get_var_clim_by_region(fjords_processed,'Ss');

% repeat the climatology for the same time period
tocn_forcing = repmat(tocn_clim,size(tocn_anom,1)/size(tocn_clim,1),1,1);
socn_forcing = repmat(socn_clim,size(socn_anom,1)/size(socn_clim,1),1,1);

% reduce anomaly to a standardised/normalised profile and a "dT factor"
% Anomaly based on surface, multiplied by an "atenuation profile"?
inds_upper_sfc = depths > -500;
tocn_anom_upper = squeeze(mean(tocn_anom(:,inds_upper_sfc,:),2));
socn_anom_upper = squeeze(mean(socn_anom(:,inds_upper_sfc,:),2));


tocn_pd = fitdist(tocn_anom_upper(:,i_reg),'kernel');    
socn_pd = fitdist(socn_anom_upper(:,i_reg),'kernel');

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = tocn_anom(:,:,i_reg);
iOpts.Marginals(end).Parameters     = tocn_anom(:,:,i_reg);
probs(end+1) = tocn_pd;

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = socn_anom(:,:,i_reg);
iOpts.Marginals(end).Parameters     = socn_anom(:,:,i_reg);
probs(end+1) = socn_pd;

%Plot to check parameter space
% figure('Name','Tanom parameter space'); histogram(random(tocn_pd,1000));
% figure('Name','Sanom parameter space'); histogram(random(socn_pd,1000));

%% Gets the glacier forcing

% Subglacial discharge (same procedure as for T and S)
% TODO: We might need something amplitude-related here, not the anomalies
% themselves, as small changes from non-summer months will dominate the signal
[q_clim,q_anom,~] = get_var_clim_by_region(fjords_processed,'Qsg');
q_pd = fitdist(q_anom(:,i_reg),'kernel');
q_forcing = repmat(q_clim,size(q_anom,1)/size(q_clim,1),1,1);

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = q_anom(:,i_reg);
iOpts.Marginals(end).Parameters     = q_anom(:,i_reg);
probs(end+1) = q_pd;

% It makes no sense to have Qsg in terms of anomalies, because in the
% frequency distribution, the very small variations in non-summer months would
% dominate the signal. So we prescribe a uniform distribution of
% "amplifying factors" instead
iOpts.Marginals(end+1).Type = 'uniform';
iOpts.Marginals(end).Parameters     = [0.5 5];
probs(end+1) = q_pd; % TODO: add uniform distribution here

%Plot to check parameter space
% figure('Name','Qsg parameter space'); histogram(random(q_pd,1000));

% Solid-ice discharge (same procedure as for T and S)
[d_clim,d_anom,~] = get_var_clim_by_region(fjords_processed,'D');
d_pd = fitdist(d_anom(:,i_reg),'kernel');
d_forcing = repmat(d_clim,size(d_anom,1)/size(d_clim,1),1,1);

iOpts.Marginals(end+1).Type = 'ks';
iOpts.Marginals(end).Inference.Data = d_anom(:,i_reg);
iOpts.Marginals(end).Parameters     = d_anom(:,i_reg);
probs(end+1) = d_pd;

%Plot to check parameter space
% figure('Name','D parameter space'); histogram(random(d_pd,1000));

%% interpolates time series variables to the actual time steps used by the model
datasets.opts.dt            = 1; % time step in days
fjord_dummy = prepare_boxmodel_input(datasets,fjords_compilation(1));
time_axis = fjord_dummy.t;

tocn_forcing = interp1(fjords_processed(1).t,tocn_forcing,time_axis,'linear','extrap');
socn_forcing = interp1(fjords_processed(1).t,socn_forcing,time_axis,'linear','extrap');
q_forcing    = interp1(fjords_processed(1).t,q_forcing,time_axis,'linear','extrap');
d_forcing    = interp1(fjords_processed(1).t,d_forcing,time_axis,'linear','extrap');

%% Compiles the parameters and variables into manageable structures

params.H    = max(fjord_stats.H.total);
params.zs   = -depths;
params.Tocn = tocn_forcing(:,:,i_reg);
params.Socn = socn_forcing(:,:,i_reg);
params.Qglc = q_forcing(:,i_reg);
params.Dglc = d_forcing(:,i_reg);
params.t    = time_axis;
params.zi   = fjord_dummy.f.zi;
params.xi   = fjord_dummy.f.xi;
params.I0   = fjord_dummy.a.I0;


end


