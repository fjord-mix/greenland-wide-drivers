%% Showing all input parameters first

% close all
% plot_reg_ocn_forcings(datasets,fjords_compilation)
% figure(1); exportgraphics(gcf,[figs_path,'hovmoller_Ts.png'],'Resolution',300)
% figure(2); exportgraphics(gcf,[figs_path,'hovmoller_Ss.png'],'Resolution',300)
% figure(3); exportgraphics(gcf,[figs_path,'series_discharge_hc_sc.png'],'Resolution',300)

% hs = plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
% hf = plot_distributions(datasets,fjords_compilation);
% exportgraphics(hf,[figs_path,'summary_input_probs2010-2018_ratios.png'],'Resolution',300)


%% Plotting the results of our analyses

% we want to know how many runs were successful
ok_runs  = zeros([1, n_regions]);
ok_vruns = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];

    % total_runs = ohc_vld(:,i_reg);
    % ok_vruns(i_reg) = sum(~isnan(total_runs));
end

% time axis for plotting the results, excluding t0
time_axis_plt = time_axis(2:end); % datetime(2010,01,15)+1:1:datetime(2018,12,15); 

% get the range of results for showing the probability distributions
ohc_x = linspace(1.2*min(ohc_out(:)),1.2*max(ohc_out(:)),1000);
osc_x = linspace(1.2*min(osc_out(:)),1.2*max(osc_out(:)),1000);

% get the heat/salt content trends from the (undisturbed) shelf
% datasets.opts.time_start = time_axis(1);
% datasets.opts.time_end   = time_axis(end);
% datasets.opts.dt         = 30.;

%% Plot the ocean forcing (and its variability) for each region
plot_reg_ocn_profiles(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg_temp','.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg_salt','.png'],'Resolution',300)

plot_reg_ocn_forcings(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'sections_ocn_forcing_reg_temp','.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'sections_ocn_forcing_reg_salt','.png'],'Resolution',300)

%% Plot the time series to see how they all behave
plot_ensemble_dt_ds(ensemble,time_axis_plt,regions_lbl);
% exportgraphics(gcf,[figs_path,'ensemble_series_mp_bootstrapped_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plots the lag between fjord and shelf
plot_ensemble_ts_lags(ensemble,360);
% exportgraphics(gcf,[figs_path,'lags_ts_mp_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting surrogate vs numerical model

% plot_model_fits(sur_model_ohc,Ynum_ohc,Ysur_ohc,Yind_ohc,Yvld_ohc,ok_runs,'temperature','^oC');
plot_model_fits(sur_model_ohc,Ynum_ohc,Ysur_ohc,[],[],ok_runs,'temperature','^oC');
% exportgraphics(gcf,[figs_path,'pce_fit_dt_n',num2str(n_runs),'.png'],'Resolution',300)

% plot_model_fits(sur_model_osc,Ynum_osc,Ysur_osc,Yind_osc,Yvld_osc,ok_runs,'salinity','');
plot_model_fits(sur_model_osc,Ynum_osc,Ysur_osc,[],[],ok_runs,'salinity','');
% exportgraphics(gcf,[figs_path,'pce_fit_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model kernel density

[hf,hp,hl] = plot_surrogate_model_results(ohc_x,osc_x,ohc_ks_eval,osc_ks_eval);
% exportgraphics(gcf,[figs_path,'ks_surrogate_dt_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% construct the numerical model kernel density plot (just for comparison)

[hf,hp,hl] = plot_numerical_model_distributions(ohc_x,osc_x,ohc_ks,osc_ks,regions_lbl);
% exportgraphics(gcf,[figs_path,'ks_boxmodel_dt_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting the Borgonovo Indices (quantifying role of inputs in each region's fjord-shelf differences for each region)
% Original publication: https://doi.org/10.1016/j.ress.2006.04.015
figure('Name','Borgonovo Indices','position',[40 40 1000 400])
for i_reg=1:7
    borgResults_ohc  = BorgonovoA_ohc{i_reg}.Results;
    borgResults_osc  = BorgonovoA_osc{i_reg}.Results;
    borgIndices = [borgResults_ohc.Delta borgResults_osc.Delta];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts{i_reg}.Marginals), borgIndices, 1.,'edgecolor','black');%,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts{i_reg}.Marginals),'XTickLabel', borgResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('Borgonovo indices'); end
    ylim([0 0.5])
end
hl = legend(hb,{'Temperature','Salinity'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
% exportgraphics(gcf,[figs_path,'borgonovo_delta_n',num2str(n_runs),'.png'],'Resolution',300)


% for i_reg=1:7
%     uq_display(BorgonovoA_ohc{1},1,'Joint PDF',3);
% end

%% Plotting the Sobol indices (quantifying role of inputs in the variability between fjords in the same region)

figure('Name','First-order Sobol Indices','position',[40 40 1000 400])
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    sobolResults_osc  = sobolA_osc{i_reg}.Results;
    sobolIndices = [sobolResults_ohc.FirstOrder sobolResults_osc.FirstOrder];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts{i_reg}.Marginals), sobolIndices, 1.,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts{i_reg}.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('First-order Sobol indices'); end
    ylim([0 1])
end
hl = legend(hb,{'Temperature','Salinity'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
% exportgraphics(gcf,[figs_path,'sobol_first_n',num2str(n_runs),'.png'],'Resolution',300)


figure('Name','Total Sobol Indices','position',[40 40 1000 400])
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    sobolResults_osc  = sobolA_osc{i_reg}.Results;
    sobolIndices = [sobolResults_ohc.Total sobolResults_osc.Total];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts{i_reg}.Marginals), sobolIndices, 1.,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts{i_reg}.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('Total Sobol indices'); end
    ylim([0 1])
end
hl = legend(hb,{'Temperature','Salinity'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
% exportgraphics(gcf,[figs_path,'sobol_total_n',num2str(n_runs),'.png'],'Resolution',300)


%% Convergence test to see if our choice of n_runs was enough

hf = plot_convergence_test(x_subsample,Yconv_ohc,Yconv_osc,ok_runs,n_runs);
% exportgraphics(gcf,[figs_path,'nruns_convergence_dt_ds_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model boxplot
Yeval_mat_ohc = NaN([n_regions,n_runs]);
Yeval_mat_osc = NaN([n_regions,n_runs]);
for i_reg=1:n_regions
    Yeval_mat_ohc(i_reg,:) = mean(Yboo_ohc{i_reg},2,'omitnan');
    Yeval_mat_osc(i_reg,:) = mean(Yboo_osc{i_reg},2,'omitnan');
end
figure; 
subplot(2,1,1); hold on; box on; grid on;
boxplot(Yeval_mat_ohc','boxstyle','filled','symbol','.','labels',regions,'colors',lines(n_regions))
xline(0,'--k');
ylabel('Temperature difference (^oC)'); text(0.02,0.95,'(a)','Units','normalized','fontsize',14)
set(gca,'fontsize',14)
subplot(2,1,2); hold on; box on; grid on;
boxplot(Yeval_mat_osc','boxstyle','filled','symbol','.','labels',regions,'colors',lines(n_regions))
xline(0,'--k')
ylabel('Salinity difference'); text(0.02,0.95,'(b)','Units','normalized','fontsize',14)
set(gca,'fontsize',14)
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
%     fprintf('Botstrapped model results for %s:\n',regions{i_reg})
%     fprintf('Mean dT: %.2f (%.2f - %.2f)\n',mean_ohc,std_ohc(1),std_ohc(2))
%     fprintf('Mean dS: %.2f (%.2f - %.2f)\n',mean_osc,std_osc(1),std_osc(2))
%     disp('')
%     fprintf('Surrogate model results for %s:\n',regions{i_reg})
%     fprintf('Mean dT: %.2f +- %.2f\n',sur_model_ohc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_ohc{i_reg}.PCE.Moments.Var))
%     fprintf('Mean dS: %.2f +- %.2f\n',sur_model_osc{i_reg}.PCE.Moments.Mean,sqrt(sur_model_osc{i_reg}.PCE.Moments.Var))
%     disp('===========================================')
% end