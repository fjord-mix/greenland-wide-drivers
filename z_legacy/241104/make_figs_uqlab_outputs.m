%% Showing all input parameters first

% close all
% plot_reg_ocn_forcings(datasets,fjords_compilation)
% figure(1); exportgraphics(gcf,[figs_path,'hovmoller_Ts.png'],'Resolution',300)
% figure(2); exportgraphics(gcf,[figs_path,'hovmoller_Ss.png'],'Resolution',300)

% hs = plot_fjords_sectors(datasets,fjords_map,fjords_compilation,glaciers_compilation);
% exportgraphics(hs.hf,[figs_path,'grl_fjords_compilation_allglaciers.png'],'Resolution',300)

% hs = plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
% hf = plot_distributions(datasets,fjords_compilation,glaciers_compilation);
% exportgraphics(hf,[figs_path,'summary_input_probs2010-2018_ratios.png'],'Resolution',300)


%% Plotting the results of our analyses

% we want to know how many runs were successful
ok_runs  = zeros([1, n_regions]);
% ok_vruns = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(1,:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];

    % total_runs = ohc_vld(:,i_reg);
    % ok_vruns(i_reg) = sum(~isnan(total_runs));
end

%% Plot the ocean forcing (and its variability) for each region
plot_reg_ocn_profiles(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg_temp','.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg_salt','.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg_dens','.eps'],'Resolution',300)

plot_reg_ocn_forcings(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'sections_ocn_forcing_reg_temp','.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'sections_ocn_forcing_reg_salt','.png'],'Resolution',300)

plot_reg_glacier_forcings(datasets,fjords_compilation,glaciers_compilation)
% exportgraphics(gcf,[figs_path,'series_glacier_discharge.png'],'Resolution',300)

plot_reg_ts(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'ts_diag_ocn_forcing_reg','.png'],'Resolution',300)

%% Addressing which "temperature metric" best links fjord-modified tempts in contact w/ glacier and fjord processes

plot_heat_dt_correlations(time_axis,ensemble,regions);
% exportgraphics(gcf,[figs_path,'xcorr_heat_transp_layers_gline_v2','.png'],'Resolution',300)

%% Plotting hovmoller of ensemble-averaged Tf and Ts per region to compare
for i_reg=3:7
    plot_hov_fjord_shelf(ensemble,time_axis,i_reg,regions);
    % figure(i_reg);
    exportgraphics(gcf,[figs_path,'hovmoller_temp_',regions{i_reg},'.png'],'Resolution',300)
end

%% Plot the time series to see how they all behave
% will also receive a formatted timetable for easier operations with dT and dS
% although not a good practice to mix processing and figure plotting, this minimises redundant code/computations
% [~,tt_ensemble] = plot_ensemble_dt_ds(ensemble,time_axis,regions_lbl);
% exportgraphics(gcf,[figs_path,'ensemble_series_mp_n',num2str(n_runs),'.png'],'Resolution',300)

[~,tmp_ensemble,tmp_spread_ensemble] = plot_ensemble_dt_ds_layers(ensemble,time_axis,2,regions);
% exportgraphics(gcf,[figs_path,'ensemble_series_mp_midlayer_n',num2str(n_runs),'.png'],'Resolution',300)

%% Show how different dT and dS are for summer and non-summer months

[~] = plot_seasonal_cycle(datasets,fjords_compilation,glaciers_compilation,tmp_ensemble,tmp_spread_ensemble,regions,1); % requires running "plot_ensemble_dt_ds" first, for 'tt_ensemble'
% exportgraphics(gcf,[figs_path,'seasonal_cycles_ensemble_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plots the lag between fjord and shelf
% plot_ensemble_ts_lags(ensemble,360); % very convoluted...
% exportgraphics(gcf,[figs_path,'lags_ts_mp_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting surrogate vs numerical model

plot_model_fits(sur_model_ohc,Ynum_ohc,Ysur_ohc,[],[],ok_runs,'temperature','^oC');
% exportgraphics(gcf,[figs_path,'pce_fit_dt_n',num2str(n_runs),'.png'],'Resolution',300)

plot_model_fits(sur_model_osc,Ynum_osc,Ysur_osc,[],[],ok_runs,'salinity','');
% exportgraphics(gcf,[figs_path,'pce_fit_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model kernel density

[hf,hp,hl] = plot_surrogate_model_results(ohc_out,osc_out,ohc_ks_eval,osc_ks_eval);
% exportgraphics(gcf,[figs_path,'prob_dist_sur_dt_ds','.png'],'Resolution',300)

%% construct the numerical model kernel density plot (just for comparison)

% [hf,hp,hl] = plot_numerical_model_distributions(ohc_out,osc_out,regions_lbl);
% exportgraphics(gcf,[figs_path,'prob_dist_boxmodel_dt_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting the Borgonovo Indices (quantifying role of inputs in each region's fjord-shelf differences for each region)
% Original publication: https://doi.org/10.1016/j.ress.2006.04.015

hf = plot_index_matrix([],BorgonovoA_ohc,'Delta',IOpts,regions,'dT Borgonovo Index \delta_i');
% exportgraphics(gcf,[figs_path,'borgonovo_delta_dT.png'],'Resolution',300)
hf = plot_index_matrix([],BorgonovoA_osc,'Delta',IOpts,regions,'dS Borgonovo Index \delta_i');
% exportgraphics(gcf,[figs_path,'borgonovo_delta_dS.png'],'Resolution',300)



%% Plotting the Sobol indices (quantifying role of inputs in the variability between fjords in the same region)

hf = plot_index_matrix([],sobolA_ohc,'FirstOrder',IOpts,regions,'dT First-order Sobol index S_i');
% exportgraphics(gcf,[figs_path,'sobol_first_dT.png'],'Resolution',300)
hf = plot_index_matrix([],sobolA_ohc,'Total',IOpts,regions,'dT Total Sobol index S');
% exportgraphics(gcf,[figs_path,'sobol_total_dT.png'],'Resolution',300)


%% Convergence test to see if our choice of n_runs was enough

hf = plot_convergence_test(x_subsample,Yconv_ohc,Yconv_osc,ok_runs,n_runs);
% exportgraphics(gcf,[figs_path,'nruns_convergence_dt_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model boxplot
% Yeval_mat_ohc = NaN([n_regions,n_runs]);
% Yeval_mat_osc = NaN([n_regions,n_runs]);
% for i_reg=1:n_regions
%     Yeval_mat_ohc(i_reg,:) = mean(Yboo_ohc{i_reg},2,'omitnan');
%     Yeval_mat_osc(i_reg,:) = mean(Yboo_osc{i_reg},2,'omitnan');
% end
% figure; 
% subplot(2,1,1); hold on; box on; grid on;
% boxplot(Yeval_mat_ohc','boxstyle','filled','symbol','.','labels',regions,'colors',lines(n_regions))
% ylim([-1 1])
% xline(0,'--k');
% ylabel('Temperature difference (^oC)'); text(0.02,0.95,'(a)','Units','normalized','fontsize',14)
% set(gca,'fontsize',14)
% subplot(2,1,2); hold on; box on; grid on;
% boxplot(Yeval_mat_osc','boxstyle','filled','symbol','.','labels',regions,'colors',lines(n_regions))
% ylim([-1 1])
% xline(0,'--k')
% ylabel('Salinity difference'); text(0.02,0.95,'(b)','Units','normalized','fontsize',14)
% set(gca,'fontsize',14)
%% Surrogate model moments

for i_reg=1:7
    fprintf('Surrogate model results for %s:\n',regions{i_reg})
    fprintf('Mean dT: %.2f +- %.2f\n',sur_model_ohc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_ohc{i_reg}.PCE.Moments.Var))
    fprintf('Mean dS: %.2f +- %.2f\n',sur_model_osc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_osc{i_reg}.PCE.Moments.Var))
    disp('===========================================')
end


%% Median and prctiles of surrogate and numerical models
for i_reg=1:7
    pct_num_ohc = prctile(Ynum_ohc{i_reg},[50, 5, 95]);
    pct_num_osc = prctile(Ynum_osc{i_reg},[50, 5, 95]);
    pct_sur_ohc = prctile(Yeval_ohc{i_reg},[50, 5, 95]);
    pct_sur_osc = prctile(Yeval_osc{i_reg},[50, 5, 95]);
    fprintf('Results for %s:\n',regions{i_reg})
    fprintf('Emulator dT: %.2f [%.2f - %.2f]\n',pct_sur_ohc)
    fprintf('Simulator dT: %.2f [%.2f - %.2f]\n',pct_num_ohc)
    disp(' ')
    fprintf('Emulator dS: %.2f [%.2f - %.2f]\n',pct_sur_osc)
    fprintf('Simulator dS: %.2f [%.2f - %.2f]\n',pct_num_osc)
    disp('===========================================')
end


%% Model accuracy from the bootstraped model runs

% for i_reg=1:7
%     mean_ohc = mean(mean(Yboo_ohc{i_reg},2),'omitnan');
%     std_ohc  = prctile(mean(Yboo_ohc{i_reg},2,'omitnan'),[5 95]);
%     mean_osc = mean(mean(Yboo_osc{i_reg},2),'omitnan');
%     std_osc  = prctile(mean(Yboo_osc{i_reg},2,'omitnan'),[5 95]);
%     fprintf('Bootstrapped model results for %s:\n',regions{i_reg})
%     fprintf('Mean dT: %.2f (%.2f - %.2f)\n',mean_ohc,std_ohc(1),std_ohc(2))
%     fprintf('Mean dS: %.2f (%.2f - %.2f)\n',mean_osc,std_osc(1),std_osc(2))
%     disp('')
%     fprintf('Surrogate model results for %s:\n',regions{i_reg})
%     fprintf('Mean dT: %.2f +- %.2f\n',sur_model_ohc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_ohc{i_reg}.PCE.Moments.Var))
%     fprintf('Mean dS: %.2f +- %.2f\n',sur_model_osc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_osc{i_reg}.PCE.Moments.Var))
%     disp('===========================================')
% end