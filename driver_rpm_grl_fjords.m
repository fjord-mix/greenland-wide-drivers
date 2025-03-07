%% Driver file for simulating the compiled fjords
clearvars
run setup_paths % Configuring paths 

plot_ensemble = 0;   % whether we want the ensemble to be plotted at the end of `run_model_loop_for_year`
n_runs        = 300; % number of runs per fjord
dt_in_h       = 1;   % model time step in hours
dt_plume_h    = 12;  % time step for updating the plume dynamics
n_years       = 10;  % how many years we want to run
tgt_days      = [n_years*365-180,n_years*365-105];  % which days of the run we want vertical profiles for

%% Compiling data 

file_fjords_compiled = [data_path,'/fjordmix/fjords_digitisation/fjords_gl_sill_depths_reduced_v2.xlsx'];
folder_ctd_casts     = [data_path,'/obs/OMG_all_casts'];
file_lengths         = [data_path,'/fjordmix/fjords_digitisation/fjords_centreline.shp'];
file_fjords          = [data_path,'/fjordmix/fjords_digitisation/fjords_grl.shp'];

fjords_digitised  = shaperead(file_fjords);
fjords_centreline = shaperead(file_lengths);

fjord_matrix = readtable(file_fjords_compiled);
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

%% Define parameter space
param_names = {'A0','Wp','C0'};
param_units = {'m^2','m','s'};

sermilik_max_bergs = 3e8;  % maximum submerged iceberg area within Sermilik fjord
% sermilik_area      = 1.1850e09; % area of Sermilik fjord (W*L)
% sermilik_vagl      = 7.7028e11; % volume above grounding line of Sermilik fjord (W*L*Hgl)
% iceberg_congestion = sermilik_max_bergs/sermilik_area;

range_params = {[0,1.5*sermilik_max_bergs],...  % A0 (if log scale, starts at 1)
                [10,700],... % wp 
                log10([5e1,5e5])};  % C0

% create probability functions using UQLab
rng('default') % set the seed for reproducibility
uqlab
for i_param=1:length(param_names)
    iOpts.Marginals(i_param).Type       = 'uniform'; % if we want to sample from a specific distribution, this can be done here
    iOpts.Marginals(i_param).Parameters = range_params{i_param};
    iOpts.Marginals(i_param).Bounds     = range_params{i_param};
    iOpts.Marginals(i_param).Name       = param_names{i_param};
end
input = uq_createInput(iOpts);

X = uq_getSample(input,n_runs,'LHS'); % create Latin Hypercube
X(:,3) = 10.^X(:,3); % reverting from log quantities to the ones we actually need
range_params{3} = 10.^range_params{3};
disp('Parameter space created.')
% plot_lhs(X,param_names,param_units,1); % quick check of input params distribution

file_inputs = [outs_path,'inputs_GRL_fjords_n',num2str(n_runs),'_dt',num2str(dt_in_h),'h.mat'];
save(file_inputs,'-v7.3','X','param_names','param_units','range_params','fjord_matrix');

%% Run the model for every year we want
for which_year=2016:2020
    [path_fout,tgt_days] = run_model_loop_for_year(which_year,fjords_digitised,fjords_centreline,fjord_matrix,...
                                                   folder_ctd_casts,X,param_names,n_years,tgt_days,dt_in_h,dt_plume_h,...
                                                   plot_ensemble);
    close all
end
delete(gcp('nocreate')) % ensure there is no parallel pool running at the end
%% Batch processing all years together
if ~exist('X','var') || ~exist('param_names','var') || ~exist('range_params','var')
    file_inputs = [outs_path,'inputs_GRL_fjords_n',num2str(n_runs),'_dt',num2str(dt_in_h),'h'];
    load(file_inputs);
end
if ~exist('path_fout','var')
    which_year=2020;
    path_fout = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year),'_dtp',num2str(12),'h_dtm',num2str(1),'h.mat'];
end
load(path_fout)
run postprocess_plot_ensembles

%% Plotting figures
i_yr_plt = 5;
% Summary of simulated fjords (Fig. 1b)
hf_fig1 = plot_best_runs_map(data_path,ensemble_yr,res_box_yr);
exportgraphics(hf_fig1,[figs_path,'sketch_processes/fjords_simulated_n',num2str(n_runs),'.png'],'BackgroundColor','none','Resolution',300)
% 
% Proof of concept that the model works (Fig. 2)
hf_ts = plot_ensemble_tempsalt(fjord_model_yr,ensemble_yr,res_box_yr,res_obs_yr,n_runs,tgt_days(2),2,{'0','28','89'},i_yr_plt);
exportgraphics(hf_ts,[figs_path,'2_temp_salt_example_fjords',num2str(2020),'_n',num2str(n_runs),'_fitall.png'],'Resolution',300)
% close all
%
% Sensitivity plots for select fjords (Fig. 3)
[hf_t,~,~] = plot_sensitivity_ensemble(X,ensemble_yr{i_yr_plt},res_box_yr{i_yr_plt},res_obs_yr{i_yr_plt},param_names,{'0','28','89'},0,0);
exportgraphics(hf_t,[figs_path,'3_sensitivity_temp_',num2str(2020),'_n',num2str(n_runs),'_v2.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'supp/sensitivity_QVs_',num2str(2020),'_n',num2str(n_runs),'.png'],'Resolution',300)
% close all
% 
% Plotting best parameters (Fig. 4)
hf_hist = plot_best_params_time_hist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,2);
exportgraphics(hf_hist,[figs_path,'4_best_params_GRL_hist_n',num2str(n_runs),'_filtered2.png'],'Resolution',300)
% close all

%% Supplementary/unused

% hf_mis_par = plot_misfits_per_parameter(X,ensemble_yr,res_box_yr,param_names,[]);  % most generalised example, showing all runs
% exportgraphics(hf_mis_par,[figs_path,'supp/rmse_vs_param_example_fjords_',num2str(n_runs),'_all.png'],'Resolution',300)
% hf_mis_par = plot_misfits_per_parameter(X,ensemble_yr(end),res_box_yr(end),param_names,{'0','28','89'}); % just showing our example fjords
% exportgraphics(hf_mis_par,[figs_path,'supp/rmse_vs_param_example_fjords_',num2str(n_runs),'.png'],'Resolution',300)

% hf_lhs = plot_lhs(X,param_names,param_units);
% exportgraphics(hf_lhs,[figs_path,'supp/lhs_design_n',num2str(n_runs),'.png'],'Resolution',300)

% If we find a distribution that fits
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,2,'kernel');
% exportgraphics(gcf,[figs_path,'best_params_GRL_prob_density',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
% close all;

% [hf_dst,hf_loc] = plot_ocn_cast_pairs(folder_ctd_casts,fjord_matrix,res_box_yr);
% exportgraphics(hf_dst,[figs_path,'supp/casts_dst_comparison/dst_OMG_fjord_shelf_casts.png'],'Resolution',300)
% exportgraphics(hf_loc,[figs_path,'overview_fig/loc_OMG_fjord_shelf_casts.png'],'Resolution',300)