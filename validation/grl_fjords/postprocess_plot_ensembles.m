n_years_to_plot = 5;
res_box_yr      = cell([n_years_to_plot,1]);
ensemble_yr     = cell(size(res_box_yr));
fjord_model_yr  = cell(size(res_box_yr));
fjord_IDs       = 0:height(fjord_matrix); %char(65:1:65+length(fjord_names)-1);

for i_yr_load=1:n_years_to_plot
    which_year_load = 2015+i_yr_load;
    file_in = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year_load),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h'];
    load(file_in);

    [~,res_box_yr{i_yr_load}] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
    ensemble_yr{i_yr_load}    = ensemble;
    fjord_model_yr{i_yr_load} = fjord_model;
    fprintf('Postprocessing %d done.\n',which_year_load)
    plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2);
    exportgraphics(gcf,[figs_path,'GRL_wide/profiles_GRL_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    close all
end
% plot_best_params_time(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2);
plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,'kernel');
exportgraphics(gcf,[figs_path,'GRL_wide/best_params_GRL_ksdensity',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
close all;
plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,2,[]);
exportgraphics(gcf,[figs_path,'GRL_wide/best_params_GRL_hist',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
close all;