function [params,iOpts,probs,fjords_processed] = define_model_param_distrib(datasets,fjords_compilation,i_reg,experiments_taxis,experiments_time_dt)
%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: W, L, Zg, Zs, Ta, Sa, Qa, Da

%% Get compiled fjord data in pre-processed structure
% also selects the desired period in time
verbose.plot=1;  % change to 1 to produce plots of all parameter distributions
verbose.print=0; % change 1 to print basic fjord statistics

if nargin < 4 || isempty(experiments_taxis)
    datasets.opts.time_start = datetime(2010,01,15);
    datasets.opts.time_end   = datetime(2018,12,15);
else
    datasets.opts.time_start = experiments_taxis(1);
    datasets.opts.time_end   = experiments_taxis(end);
end
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
if nargin < 5
    experiments_time_dt       = 0.10; % time step (in days) for the actual experiments
end
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

%% Gets distribution of fjord geometry parameters
% since the geometry does not change with time, we can stick to the
% original fjords structure
fjord_stats = print_fjord_statistics(fjords_compilation,verbose);

probs(1) = fjord_stats.L.pd;
iOpts.Marginals(1) = uq_KernelMarginals(fjord_stats.L.total', [min(fjord_stats.L.total) max(fjord_stats.L.total)]);

probs(end+1) = fjord_stats.H.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.H.total',[min(fjord_stats.H.total), max(fjord_stats.H.total)]);

% probs(end+1) = fjord_stats.W.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.W.total',[min(fjord_stats.W.total), max(fjord_stats.W.total)]);
probs(end+1) = fjord_stats.a.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.a.total',[1, max(fjord_stats.a.total)]);

% probs(end+1) = fjord_stats.Zs.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zs.total', [min(fjord_stats.Zs.total) -50]);
probs(end+1) = fjord_stats.Zsr.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zsr.total', [min(abs(fjord_stats.Zs.total)/max(fjord_stats.H.total)) max(abs(fjord_stats.Zs.total)/max(fjord_stats.H.total))]);
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zsr.total', [0.05 0.95]);

% probs(end+1) = fjord_stats.Zg.pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zg.total', [min(fjord_stats.Zg.total) -50]);
probs(end+1) = fjord_stats.Zgr.pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zgr.total', [min(abs(fjord_stats.Zg.total)/max(fjord_stats.H.total)) max(abs(fjord_stats.Zg.total)/max(fjord_stats.H.total))]);
% iOpts.Marginals(end+1) = uq_KernelMarginals(fjord_stats.Zgr.total', [0.05 0.95]);

%% Gets the ocean forcing

% [tocn_forcing, tocn_anom, tocn_decay, depths] = get_var_clim_by_region(fjords_processed,'Ts');
% [socn_forcing, socn_anom, socn_decay, ~]      = get_var_clim_by_region(fjords_processed,'Ss');
[tocn_forcing, tocn_anom, tocn_decay, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
[socn_forcing, socn_anom, socn_decay, ~]      = get_var_forcing_by_region(fjords_processed,'Ss');

% reduce anomaly to a standardised/normalised profile and a "dT factor"
% Anomaly based on surface, multiplied by an "atenuation profile/decay function"
tocn_pd = fitdist(tocn_anom{i_reg},'kernel');    
socn_pd = fitdist(socn_anom{i_reg},'kernel');
omeg_pd = makedist('Normal','mu',0.1429,'sigma',0.1167);

probs(end+1) = tocn_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(tocn_anom{i_reg},[min(tocn_anom{i_reg}), max(tocn_anom{i_reg})]);

probs(end+1) = socn_pd;
iOpts.Marginals(end+1) = uq_KernelMarginals(socn_anom{i_reg},[min(socn_anom{i_reg}), max(socn_anom{i_reg})]);

% only applies omega to eastern Greenland
if ismember(Parameters.regID,[2,4,6]) 
    iOpts.Marginals(end+1).Type     = 'Gaussian';
    iOpts.Marginals(end).Parameters = [0.1429 0.1167]; % [1/7 days, 3.5/30 days] as per Harden et al. (2011; JClimate)
    iOpts.Marginals(end).Bounds     = [0 1];
else
    iOpts.Marginals(end+1).Type = 'Constant' ;
    iOpts.Marginals(end).Parameters = 0;
end

probs(end+1) = omeg_pd;
%% Gets the glacier forcing

% Subglacial discharge - we want the variation in amplitude rather than
% anomalies; otherwise the signal will be dominated by small variations in non-summer months
str_reg = {'SW','SE','CW','CE','NW','NE','NO'};
q_amp_max=NaN(size(fjords_processed));
q_amp_min=NaN(size(fjords_processed));
for k=1:length(fjords_processed)
    if strcmp(fjords_processed(k).m.region,str_reg{i_reg})
        [v_peaks,~] = findpeaks(fjords_processed(k).f.Qsg,'MinPeakProminence',10);
        q_amp_max(k) = max(v_peaks)./mean(v_peaks);
        q_amp_min(k) = min(v_peaks)./mean(v_peaks);
    end
end
q_amp_max = mean(q_amp_max,'omitnan');
q_amp_min = mean(q_amp_min,'omitnan');
% [q_clim,q_anom,~,~] = get_var_clim_by_region(fjords_processed,'Qsg');
[q_forcing,~,~,~] = get_var_forcing_by_region(fjords_processed,'Qsg');
q_pd              = makedist('Uniform','lower',q_amp_min,'upper',q_amp_max);

% It makes no sense to have Qsg in terms of anomalies, because in the
% frequency distribution, the very small variations in non-summer months would
% dominate the signal. So we prescribe a uniform distribution of
% "amplifying factors" instead
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [q_amp_min q_amp_max];
iOpts.Marginals(end).Bounds     = [0 1.5*q_amp_max];
probs(end+1)                    = q_pd;

%Plot to check parameter space
% figure('Name','Qsg parameter space'); histogram(random(q_pd,1000));

% considering we use an entrainment coefficient of 0.1, P0=[5,30] is
% equivalent to a plume width of 50-300 m
p_pd              = makedist('Uniform','lower',5,'upper',30); 
iOpts.Marginals(end+1).Type     = 'uniform';
iOpts.Marginals(end).Parameters = [5 30];
iOpts.Marginals(end).Bounds     = [5 30];
probs(end+1)                    = p_pd;

% Solid-ice discharge (same procedure as for T and S)
% [d_forcing,d_anom,~,~] = get_var_clim_by_region(fjords_processed,'D');
[d_forcing,d_anom,~,~] = get_var_forcing_by_region(fjords_processed,'D');
% d_pd                   = fitdist(d_anom{i_reg},'kernel');
% Ignoring D (for now?)
% probs(end+1) = d_pd;
% iOpts.Marginals(end+1) = uq_KernelMarginals(d_anom{i_reg},[min(d_anom{i_reg}), max(d_anom{i_reg})]);


%% interpolates time series variables to the actual time steps used by the model
datasets.opts.dt            = experiments_time_dt; % time step in days
fjord_dummy = prepare_boxmodel_input(datasets,fjords_compilation(1));
time_axis = fjord_dummy.t;

tocn_forcing = interp1(fjords_processed(1).t,tocn_forcing,time_axis,'linear','extrap');
tocn_decay = interp1(fjords_processed(1).t,tocn_decay,time_axis,'linear','extrap');
socn_forcing = interp1(fjords_processed(1).t,socn_forcing,time_axis,'linear','extrap');
socn_decay = interp1(fjords_processed(1).t,socn_decay,time_axis,'linear','extrap');
q_forcing    = interp1(fjords_processed(1).t,q_forcing,time_axis,'linear','extrap');
d_forcing    = interp1(fjords_processed(1).t,d_forcing,time_axis,'linear','extrap');

%% Compiles the parameters into the structure to be used by the wrapper function

params.regID = i_reg;
params.H     = max(fjord_stats.H.total);
params.zs    = -depths;
params.Tocn  = tocn_forcing(:,:,i_reg);
params.Tdec  = tocn_decay(:,:,i_reg);
params.Socn  = socn_forcing(:,:,i_reg);
params.Sdec  = socn_decay(:,:,i_reg);
params.Qglc  = q_forcing(:,i_reg);
params.Dglc  = d_forcing(:,i_reg) .*0; % setting to zero (for now?)
params.t     = time_axis;
params.zi    = fjord_dummy.f.zi;
params.xi    = fjord_dummy.f.xi;
params.I0    = fjord_dummy.a.I0;

% Add the names to the input distributions we created
iOpts.Marginals(1).Name = 'L';
% iOpts.Marginals(2).Name = 'W';
iOpts.Marginals(2).Name = 'H';
iOpts.Marginals(3).Name = 'L/W'; % iOpts.Marginals(2).Name = 'W';
iOpts.Marginals(4).Name = 'Zs/H';
iOpts.Marginals(5).Name = 'Zg/H';
% iOpts.Marginals(3).Name = 'Zs';
% iOpts.Marginals(4).Name = 'Zg';
iOpts.Marginals(6).Name = 'Ta';
iOpts.Marginals(7).Name = 'Sa';
iOpts.Marginals(8).Name = 'omega';
iOpts.Marginals(9).Name = 'Qa';
iOpts.Marginals(10).Name = 'P0';
% iOpts.Marginals(11).Name = 'Da';

end


