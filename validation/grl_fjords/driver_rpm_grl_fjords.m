%% Driver file for simulating the compiled fjords
clearvars
run setup_paths % Configuring paths

plot_ensemble = 0;
n_runs        = 400;         % number of runs per fjord
dt_in_h       = 3;
n_years       = 10;           % how many years we want to run
tgt_days      = [n_years*365-180,n_years*365-105];  % which days of the run we want vertical profiles for

%% Define parameter space
param_units = {'m^2','m','s^{-1}','-'};
param_names = {'A0','wp','C0','K0'};

range_params = {[0,3e8],...    % A0
                [10,700],...   % wp no crashes with [10,400]
                [1e2,5e4],...  % C0
                [1e-4,1e-3]};  % K0 or do we stick to [1e4,1e-3]?

rng('default') % set the seed for reproducibility
uqlab
for i_param=1:length(param_names)
    iOpts.Marginals(i_param).Type       = 'uniform';
    iOpts.Marginals(i_param).Parameters = range_params{i_param};
    iOpts.Marginals(i_param).Bounds     = range_params{i_param};
    iOpts.Marginals(i_param).Name       = param_names{i_param};
end

input = uq_createInput(iOpts); % create probability functions
X = uq_getSample(input,n_runs,'LHS');  % training dataset
disp('Parameter space created.')

%% Compiling data 

% Read compilation of all fjords around Greenland (requires data-compilation repository)
% only needed if not using fjord geometries from Cowton et al.
% run compile_process_fjords 
% fjord_ids = [4,9,17,20,22,23,24,25,28,29,30,31,24,37]; % These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
file_fjords_compiled = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_gl_sill_depths_reduced_v2.xlsx'];
folder_ctd_casts     = [data_path,'/greenland/obs/OMG_all_casts'];
file_lengths = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_centreline.shp'];
file_fjords = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_grl.shp'];


fjords_digitised  = shaperead(file_fjords);
fjords_centreline = shaperead(file_lengths);

fjord_matrix = readtable(file_fjords_compiled);
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

%% Run the model for every year we want
for which_year=2016:2020
    [path_fout,tgt_days] = run_model_loop_for_year(which_year,fjords_digitised,fjords_centreline,fjord_matrix,...
                                                   folder_ctd_casts,X,param_names,n_years,tgt_days,dt_in_h,...
                                                   plot_ensemble);
    close all
end

%% Batch processing all years together
% path_fout = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year),'_',num2str(60),'layers_dt',num2str(3),'h'];
% load(path_fout)
% run postprocess_plot_ensembles.m

% hf_fw = plot_hist_zfw_export(ensemble_yr,res_box_yr);
% exportgraphics(hf_fw,[figs_path,'2_hist_fw_export.png'],'Resolution',300)

%% Plotting best parameters
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,2);
% exportgraphics(gcf,[figs_path,'best_params_GRL_hist_n',num2str(n_runs),'_ign20.png'],'Resolution',300)
% close all;

% plot_best_params_time(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,2);
% exportgraphics(gcf,[figs_path,'supp/best_params_GRL_scatter_allyrs_n',num2str(n_runs),'.png'],'Resolution',300)
% plot_best_params_rmse(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
% exportgraphics(gcf,[figs_path,'best_params_GRL_rmse',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)

% If we find a distribution that fits
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,2,'kernel');
% exportgraphics(gcf,[figs_path,'best_params_GRL_prob_density',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
% close all;


% [hf_dst,hf_loc] = plot_ocn_cast_pairs(folder_ctd_casts,fjord_matrix,res_box_yr);
% exportgraphics(hf_dst,[figs_path,'supp/casts_dst_comparison/dst_OMG_fjord_shelf_casts.png'],'Resolution',300)
% exportgraphics(hf_loc,[figs_path,'overview_fig/loc_OMG_fjord_shelf_casts.png'],'Resolution',300)