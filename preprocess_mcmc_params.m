%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: W, L, Zg, Zs, Ta, Sa, Da, Qa
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/mcmc/']; % where the model output files will be saved
figs_path = [project_path,'/figs/mcmc/'];                     % where the figures and animations will be saved

[datasets,fjords_compilation,~,~] = compile_datasets(data_path);

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
fjord_stats = print_fjord_statistics(fjords_compilation);%,verbose);

W_space  = fjord_stats.W.pd.random(10000,1);
L_space  = fjord_stats.L.pd.random(10000,1);
Zg_space = fjord_stats.Zg.pd.random(10000,1);
Zs_space = fjord_stats.Zs.pd.random(10000,1);

% figure for checking how the distribution looks
% figure; 
% subplot(2,2,1); histogram(L_space);
% subplot(2,2,2); histogram(W_space);
% subplot(2,2,3); histogram(Zg_space);
% subplot(2,2,4); histogram(Zs_space);

%% Sets up the ocean forcing

% [tocn_clim, socn_clim, tocn_anom, socn_anom, depths] = get_ocean_clim_by_region(fjords_processed);

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

% Compute PDF of anomalies
tanom_space = cell(1,7);
sanom_space = cell(1,7);
tocn_pd     = cell(1,7);
socn_pd     = cell(1,7);
for i_reg=1:7
    tocn_pd{i_reg} = fitdist(tocn_anom_upper(:,i_reg),'kernel');    
    socn_pd{i_reg} = fitdist(socn_anom_upper(:,i_reg),'kernel');
    tanom_space{i_reg} = tocn_pd{i_reg}.random(10000,1);
    sanom_space{i_reg} = socn_pd{i_reg}.random(10000,1);
end

% Plot to check resulting PDFs
% iplt=1;
% figure('T and S parameter space');
% for i_reg=1:7
%     subplot(7,2,iplt);
%     %ksdensity(tocn_anom_upper(:,i_reg));
%     histogram(tanom_space{i_reg});
%     subplot(7,2,iplt+1)
%     %ksdensity(socn_anom_upper(:,i_reg));
%     histogram(sanom_space{i_reg});
%     iplt=iplt+2;
% end

%% Sets up the glacier forcing

% Subglacial discharge (same procedure as for T and S)
[q_clim,q_anom,~] = get_var_clim_by_region(fjords_processed,'Qsg');
qanom_space = cell(1,7);
q_pd        = cell(1,7);
for i_reg=1:7
    q_pd{i_reg} = fitdist(q_anom(:,i_reg),'kernel');
    qanom_space{i_reg} = q_pd{i_reg}.random(10000,1);
end

%Plot to check parameter space
% figure('Name','Qsg parameter space'); for i_reg=1:7,subplot(7,1,i_reg); histogram(tanom_space{i_reg}); end

% Solid-ice discharge (same procedure as for T and S)
[d_clim,d_anom,~] = get_var_clim_by_region(fjords_processed,'D');
danom_space = cell(1,7);
d_pd        = cell(1,7);
for i_reg=1:7
    d_pd{i_reg} = fitdist(d_anom(:,i_reg),'kernel');
    danom_space{i_reg} = d_pd{i_reg}.random(10000,1);
end

%Plot to check parameter space
% figure('Name','Ice discharge parameter space'); for i_reg=1:7,subplot(7,1,i_reg); histogram(danom_space{i_reg}); end