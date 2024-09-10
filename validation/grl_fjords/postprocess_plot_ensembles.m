n_years_to_plot = 5;
res_box_yr      = cell([n_years_to_plot,1]);
res_obs_yr     = cell(size(res_box_yr));
ensemble_yr     = cell(size(res_box_yr));
fjord_model_yr  = cell(size(res_box_yr));
fjord_IDs       = 0:height(fjord_matrix); %char(65:1:65+length(fjord_names)-1);

for i_yr_load=1:n_years_to_plot
    which_year_load = 2015+i_yr_load;
    file_in = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year_load),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h'];
    load(file_in);

    [res_obs_yr{i_yr_load},res_box_yr{i_yr_load}] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
    idx = arrayfun(fun,res_obs_yr{i_yr_load});
    res_obs_yr{i_yr_load}(idx)=[]; % remove the empty elements
    idx = arrayfun(fun,res_box_yr{i_yr_load});
    res_box_yr{i_yr_load}(idx)=[]; % remove the empty elements

    ensemble_yr{i_yr_load}    = ensemble;
    fjord_model_yr{i_yr_load} = fjord_model;
    fprintf('Postprocessing %d done.\n',which_year_load)
    % plot_ensemble_profiles(fjord_model,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},res_obs_yr{i_yr_load},n_runs,param_names,tgt_days(2),name_days,2);
    plot_ensemble_profiles(fjord_model,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},res_obs_yr{i_yr_load},n_runs,param_names,tgt_days(2),name_days,2,[],0,0,0,[],1);
    % exportgraphics(gcf,[figs_path,'profiles_GRL_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(gcf,[figs_path,'rmse_temp_rpm_shelf_GRL_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % plot_sensitivity_profiles_v3(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},res_obs_yr{i_yr_load},param_names,2,0,[],i_yr_load, which_fj_sens{i_yr_load});
    % exportgraphics(gcf,[figs_path,'sensitivity_profiles_temp_simple_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    close all
end
plot_best_params_time(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
exportgraphics(gcf,[figs_path,'best_params_GRL_scatter',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
plot_best_params_rmse(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
exportgraphics(gcf,[figs_path,'best_params_GRL_rmse',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
% exportgraphics(gcf,[figs_path,'best_params_GRL_hist_n',num2str(n_runs),'.png'],'Resolution',300)
% close all;
% plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,'kernel');
% exportgraphics(gcf,[figs_path,'best_params_GRL_prob_density',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
% close all;