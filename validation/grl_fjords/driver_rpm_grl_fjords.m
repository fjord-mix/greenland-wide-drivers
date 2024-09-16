%% Driver file for simulating the compiled fjords
clearvars
run setup_paths % Configuring paths

for year=2016:2020
    run_model_loop_for_year(year); % {2016,2017,2018,2019,2020}
    close all
end

% plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(1),[],1,[],1,1,0)

% plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2,[],0,0,1);
% exportgraphics(gcf,[figs_path,'rmse_temp_rpm_shelf_GRL_',num2str(which_year),'_n',num2str(n_runs),'.png'],'Resolution',300)
%% Checking how results compare from year to year

% plot_best_params_dist_qq(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,'poisson');
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
% param_distributions = {'hn','hn','wbl','hn'};
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,param_distributions);

%% Batch processing all years together

% run postprocess_plot_ensembles.m