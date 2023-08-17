function [params,iOpts,probs,fjords_processed] = define_model_param_distrib(datasets,fjords_compilation,i_reg)
%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: W, L, Zg, Zs, Ta, Sa, Da, Qa

%% Get compiled fjord data in pre-processed structure
% also selects the desired period in time
verbose.plot=0;  % change to 1 to produce plots of all parameter distributions
verbose.print=0; % change 1 to print basic fjord statistics

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
fjord_stats = print_fjord_statistics(fjords_compilation,verbose);

probs(1) = fjord_stats.L.pd;
iOpts.Marginals(1) = uq_KernelMarginals(fjord_stats.L.total', [0 max(fjord_stats.L.total)]);
iOpts.Marginals(1).Name = 'length';

probs(end+1) = fjord_stats.W.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.W.total',[0, max(fjord_stats.W.total)]);
iOpts.Marginals(end).Name = 'width';

probs(end+1) = fjord_stats.Zs.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zs.total', [min(fjord_stats.Zs.total) -50]);
iOpts.Marginals(end).Name = 'z_sill';

probs(end+1) = fjord_stats.Zg.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zg.total', [min(fjord_stats.Zg.total) -50]);
iOpts.Marginals(end).Name = 'z_gl';

% Xgeom = [fjord_stats.H.total;-fjord_stats.Zg.total;-fjord_stats.Zs.total];
% opts.Inference.Data = Xgeom';
% opts.Copula.Inference.Data = Xgeom;
% opts.Copula.Type='auto';
% InputHat = uq_createInput(opts);


%% Gets the ocean forcing

[tocn_forcing, tocn_anom, tocn_decay, depths] = get_var_clim_by_region(fjords_processed,'Ts');
[socn_forcing, socn_anom, socn_decay, ~]      = get_var_clim_by_region(fjords_processed,'Ss');


% reduce anomaly to a standardised/normalised profile and a "dT factor"
% Anomaly based on surface, multiplied by an "atenuation profile"?
tocn_pd = fitdist(tocn_anom(:,i_reg),'kernel');    
socn_pd = fitdist(socn_anom(:,i_reg),'kernel');
omeg_pd = makedist('Normal','mu',0.1429,'sigma',0.1167);

probs(end+1) = tocn_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(tocn_anom(:,i_reg),[min(tocn_anom(:,i_reg)), max(tocn_anom(:,i_reg))]);
iOpts.Marginals(end).Name = 't_anom';

probs(end+1) = socn_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(socn_anom(:,i_reg),[min(socn_anom(:,i_reg)), max(socn_anom(:,i_reg))]);
iOpts.Marginals(end).Name = 's_anom';

probs(end+1) = omeg_pd;
iOpts.Marginals(end+1).Type     = 'Normal';
iOpts.Marginals(end).Parameters = [0.1429 0.1167];
iOpts.Marginals(end).Name = 'omega';

%% Gets the glacier forcing

% Subglacial discharge (same procedure as for T and S)
% TODO: We might need something amplitude-related here, not the anomalies
% themselves, as small changes from non-summer months will dominate the signal
[q_clim,q_anom,~,~] = get_var_clim_by_region(fjords_processed,'Qsg');
q_forcing         = repmat(q_clim,size(q_anom,1)/size(q_clim,1),1,1);
q_pd              = makedist('Uniform','lower',0.5,'upper',1.5);

% It makes no sense to have Qsg in terms of anomalies, because in the
% frequency distribution, the very small variations in non-summer months would
% dominate the signal. So we prescribe a uniform distribution of
% "amplifying factors" instead
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [0.5 1.5];
iOpts.Marginals(end).Name       = 'q_amp';
probs(end+1)                    = q_pd;

%Plot to check parameter space
% figure('Name','Qsg parameter space'); histogram(random(q_pd,1000));

% Solid-ice discharge (same procedure as for T and S)
[d_forcing,d_anom,~,~] = get_var_clim_by_region(fjords_processed,'D');
d_pd              = fitdist(d_anom(:,i_reg),'kernel');

probs(end+1) = d_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(d_anom(:,i_reg),[min(d_anom(:,i_reg)), max(d_anom(:,i_reg))]);
iOpts.Marginals(end).Name = 'd_anom';

p_pd              = makedist('Uniform','lower',15,'upper',35);
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [15 35];
iOpts.Marginals(end).Name       = 'P0';
probs(end+1)                    = p_pd;

%% interpolates time series variables to the actual time steps used by the model
datasets.opts.dt            = 0.2; % time step in days
fjord_dummy = prepare_boxmodel_input(datasets,fjords_compilation(1));
time_axis = fjord_dummy.t;

tocn_forcing = interp1(fjords_processed(1).t,tocn_forcing,time_axis,'linear','extrap');
tocn_decay = interp1(fjords_processed(1).t,tocn_decay,time_axis,'linear','extrap');
socn_forcing = interp1(fjords_processed(1).t,socn_forcing,time_axis,'linear','extrap');
socn_decay = interp1(fjords_processed(1).t,socn_decay,time_axis,'linear','extrap');
q_forcing    = interp1(fjords_processed(1).t,q_forcing,time_axis,'linear','extrap');
d_forcing    = interp1(fjords_processed(1).t,d_forcing,time_axis,'linear','extrap');

%% Compiles the parameters into the structure to be used by the wrapper function

params.H    = max(fjord_stats.H.total);
params.zs   = -depths;
params.Tocn = tocn_forcing(:,:,i_reg);
params.Tdec = tocn_decay(:,:,i_reg);
params.Socn = socn_forcing(:,:,i_reg);
params.Sdec = socn_decay(:,:,i_reg);
params.Qglc = q_forcing(:,i_reg);
params.Dglc = d_forcing(:,i_reg);
params.t    = time_axis;
params.zi   = fjord_dummy.f.zi;
params.xi   = fjord_dummy.f.xi;
params.I0   = fjord_dummy.a.I0;

% Add the names to the input distributions we created
iOpts.Marginals(1).Name = 'length';
iOpts.Marginals(2).Name = 'width';
iOpts.Marginals(3).Name = 'z_sill';
iOpts.Marginals(4).Name = 'z_gl';
iOpts.Marginals(5).Name = 't_anom';
iOpts.Marginals(6).Name = 's_anom';
iOpts.Marginals(7).Name = 'omega';
iOpts.Marginals(8).Name = 'q_sg';
iOpts.Marginals(9).Name = 'q_ice';

if isfield(verbose,'plot') && verbose.plot
    run plot_distributions.m
end

end


