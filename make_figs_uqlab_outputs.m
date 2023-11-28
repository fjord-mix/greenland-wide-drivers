%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs  = zeros([1, n_regions]);
ok_vruns = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];

    total_runs = ohc_vld(:,i_reg);
    ok_vruns(i_reg) = sum(~isnan(total_runs));
end

% get the range of results for computing the probability distributions
ohc_x = linspace(1.2*min(ohc_out(:)),1.2*max(ohc_out(:)),1000);
osc_x = linspace(1.2*min(osc_out(:)),1.2*max(osc_out(:)),1000);

% get the heat/salt content trends from the (undisturbed) shelf
% datasets.opts.time_start = time_axis(1);
% datasets.opts.time_end   = time_axis(end);
% datasets.opts.dt         = 30.;
% fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
% for i=1:length(fjords_compilation),fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));end



region_line_color = lines(7);
time_axis_plt = datetime(2010,01,15)+1:1:datetime(2018,12,15);

%% Plot the ocean forcing (and its variability) for each region
plot_reg_ocn_profiles(datasets,fjords_compilation)
% exportgraphics(gcf,[figs_path,'profiles_ocn_forcing_reg,','.png'],'Resolution',300)

%% Plot the time series to see how they all behave
plot_ensemble_dt_ds(ensemble,time_axis_plt,regions_lbl);
% exportgraphics(gcf,[figs_path,'ensemble_series_bootstrapped_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plots the lag between fjord and shelf
plot_ensemble_ts_lags(ensemble,180);

%% Plotting surrogate vs numerical model

plot_model_fits(Ynum_ohc,Ysur_ohc,Yind_ohc,Yvld_ohc,ok_runs,ok_vruns,'temperature','^oC');
% exportgraphics(gcf,[figs_path,'pce_fit_ohc_n',num2str(n_runs),'.png'],'Resolution',300)

plot_model_fits(Ynum_osc,Ysur_osc,Yind_osc,Yvld_osc,ok_runs,ok_vruns,'salinity','');
% exportgraphics(gcf,[figs_path,'pce_fit_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model kernel density

[hf,hp,hl] = plot_surrogate_model_results(ohc_x,osc_x,ohc_ks_eval,osc_ks_eval);
% exportgraphics(gcf,[figs_path,'kssur_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% construct the numerical model kernel density plot (just for comparison)

[hf,hp,hl] = plot_numerical_model_distributions(ohc_x,osc_x,ohc_ks,osc_ks,regions_lbl);
% exportgraphics(gcf,[figs_path,'ksnum_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting the Sobol indices

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

%% Plotting the Borgonovo Indices
figure('Name','Borgonovo Indices','position',[40 40 1000 400])
for i_reg=1:7
    borgResults_ohc  = BorgonovoA_ohc{i_reg}.Results;
    borgResults_osc  = BorgonovoA_osc{i_reg}.Results;
    borgIndices = [borgResults_ohc.Delta borgResults_osc.Delta];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts{i_reg}.Marginals), borgIndices, 1.,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts{i_reg}.Marginals),'XTickLabel', borgResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('Borgonovo indices'); end
    ylim([0 0.5])
end
hl = legend(hb,{'Temperature','Salinity'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;

%% Convergence test to see if our choice of n_runs was enough

hf = plot_convergence_test(x_subsample,Yconv_ohc,Yconv_osc,ok_runs);
% exportgraphics(gcf,[figs_path,'nruns_convergence_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)
