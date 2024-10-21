clear_empty     = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array
n_years_to_plot = 5;
if ~exist('ensemble_yr','var')
    res_box_yr      = cell([n_years_to_plot,1]);
    res_obs_yr      = cell(size(res_box_yr));
    ensemble_yr     = cell(size(res_box_yr));
    fjord_model_yr  = cell(size(res_box_yr));
end
fjord_IDs       = 0:height(fjord_matrix); %char(65:1:65+length(fjord_names)-1);

which_fj_sens = {{'0','17','30','108'}; % 2016
                {'0','14','28','108'};  % 2017
                {'17','24','89','108'};  % 2018 
                {'8','12','30','108'};  % 2019 
                {'0','28','89','108'}};  % 2020 

for i_yr_load=1:n_years_to_plot
    %% Loading and postprocessing
    which_year_load = 2015+i_yr_load;
    if isempty(ensemble_yr{i_yr_load})
        file_in = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year_load),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h'];
        load(file_in);
    
        [res_obs_yr{i_yr_load},res_box_yr{i_yr_load}] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
        idx = arrayfun(clear_empty,res_obs_yr{i_yr_load});
        res_obs_yr{i_yr_load}(idx)=[]; % remove the empty elements
        idx = arrayfun(clear_empty,res_box_yr{i_yr_load});
        res_box_yr{i_yr_load}(idx)=[]; % remove the empty elements
        ensemble_yr{i_yr_load}    = ensemble;
        fjord_model_yr{i_yr_load} = fjord_model;
    end
    fprintf('Postprocessing %d done.\n',which_year_load)

    %% Plotting and saving the main figures

    % Profiles and series for all fjords
    % [hf_profiles,hfs_profiles,hf_series,hf_rmse] = plot_ensemble_profiles(fjord_model,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},res_obs_yr{i_yr_load},n_runs,param_names,tgt_days(2),[],2,[],1,1,0,[],1);
    % exportgraphics(hf_profiles,[figs_path,'supp/profiles_temp_salt/profiles_GRL_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hfs_profiles,[figs_path,'supp/profiles_temp_salt/profiles_GRL_salt_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_series,[figs_path,'supp/series/series_GRL_temp_',num2str(which_year_load),'_n',num2str(n_runs),'_ign20.png'],'Resolution',300)
    % exportgraphics(hf_rmse,[figs_path,'supp/rmse/rmse_temp_rpm_shelf_GRL_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % close all

    % % Sensitivity plots for select fjords
    % [hf_t,hf_e] = plot_sensitivity_ensemble(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},param_names,which_fj_sens{i_yr_load});
    % exportgraphics(hf_t,[figs_path,'supp/sensitivity_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_e,[figs_path,'supp/sensitivity_fwex_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % close all;

    %% Other, older, plots
    % Only plotting and saving temperature for select fjords
    % plot_ensemble_profiles(fjord_model,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},res_obs_yr{i_yr_load},n_runs,param_names,tgt_days(2),[],2,[],0,0,0,which_fj_sens{i_yr_load},0);
    % exportgraphics(gcf,[figs_path,'supp/profiles_ts_selected/profiles_GRLsel_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % close all

    % Plotting sensitivity profiles - fjord cast time for T and S, peak fw export for freshwater export
    % [hf_t,hf_s,~] = plot_sensitivity_profiles_v3(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},param_names,2,1,0,which_fj_sens{i_yr_load});
    % [~,~,hf_e]    = plot_sensitivity_profiles_v3(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},param_names,1,0,1,which_fj_sens{i_yr_load});
    % exportgraphics(hf_t,[figs_path,'supp/sens_profiles/sensitivity_profiles_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_s,[figs_path,'supp/sens_profiles/sensitivity_profiles_salt_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_e,[figs_path,'supp/sens_profiles/sensitivity_profiles_fwfl_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % close all

    % Plotting sensitivity scatter
    % [hf_t,hf_s,~,~] = plot_sensitivity_scatter(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},param_names,2,1,0,0,which_fj_sens{i_yr_load});
    % [~,~,hf_f,hf_z] = plot_sensitivity_scatter(X,ensemble_yr{i_yr_load},res_box_yr{i_yr_load},param_names,1,0,1,1,which_fj_sens{i_yr_load});
    % exportgraphics(hf_t,[figs_path,'supp/sens_scatter/temp/sensitivity_scatter_temp_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_s,[figs_path,'supp/sens_scatter/salt/sensitivity_scatter_salt_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_f,[figs_path,'supp/sens_scatter/fw_export/sensitivity_scatter_fwexport_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % exportgraphics(hf_z,[figs_path,'supp/sens_scatter/fw_z/sensitivity_scatter_zfwexport_',num2str(which_year_load),'_n',num2str(n_runs),'.png'],'Resolution',300)
    % close all

end
% close all