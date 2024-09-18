%% Driver file for simulating the compiled fjords
clearvars
run setup_paths % Configuring paths

plot_ensemble = 1;
n_runs        = 100;         % number of runs per fjord
dt_in_h       = 3;
%% Define parameter space
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
file_fjords_compiled = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_gl_sill_depths_reduced.xlsx'];
folder_ctd_casts     = [data_path,'/greenland/obs/OMG_all_casts'];
file_lengths = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_centreline.shp'];
file_fjords = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_grl.shp'];


fjords_digitised  = shaperead(file_fjords);
fjords_centreline = shaperead(file_lengths);

fjord_matrix = readtable(file_fjords_compiled);
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

%% Run the model for every year we want
for year=2016:2020
    [path_fout,tgt_days] = run_model_loop_for_year(year,fjords_digitised,fjords_centreline,fjord_matrix,folder_ctd_casts,X,param_names,dt_in_h,plot_ensemble); % {2016,2017,2018,2019,2020}
    close all
end

% plot_ensemble_profiles(fjord_model_yr{end},ensemble_yr{end},res_box_yr{end},res_obs_yr{end},n_runs,param_names,tgt_days(2),[],1,[],1,1,0)
% exportgraphics(gcf,[figs_path,'supp/series_temp_GRL_',num2str(year),'_n',num2str(n_runs),'.png'],'Resolution',300)

% plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2,[],0,0,1);
% exportgraphics(gcf,[figs_path,'rmse_temp_rpm_shelf_GRL_',num2str(which_year),'_n',num2str(n_runs),'.png'],'Resolution',300)
%% Checking how results compare from year to year

% plot_best_params_dist_qq(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,'poisson');
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
% param_distributions = {'hn','hn','wbl','hn'};
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,param_distributions);

%% Batch processing all years together
% path_fout = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(year),'_',num2str(60),'layers_dt',num2str(3),'h'];
% load(path_fout)
% run postprocess_plot_ensembles.m