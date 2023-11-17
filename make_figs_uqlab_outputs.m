%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs  = zeros([1, n_regions]);
% ok_vruns = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];

    % total_runs = ohc_vld(:,i_reg);
    % ok_vruns(i_reg) = sum(~isnan(total_runs));
end

% get the range of results for computing the probability distributions
ohc_x = linspace(1.2*min(ohc_out(:)),1.2*max(ohc_out(:)),1000);
osc_x = linspace(1.2*min(osc_out(:)),1.2*max(osc_out(:)),1000);

% get the heat/salt content trends from the (undisturbed) shelf
datasets.opts.time_start = time_axis(1);
datasets.opts.time_end   = time_axis(end);
datasets.opts.dt         = 30.;
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation),fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));end
[temp_forcing, ~, ~, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
[salt_forcing, ~, ~, ~] = get_var_forcing_by_region(fjords_processed,'Ss');
% fjord_rho = (fjords_processed(1).p.betaS*salt_forcing - fjords_processed(1).p.betaT*temp_forcing);
% sc_reg = squeeze(trapz(depths,salt_forcing.* fjord_rho,2)./max(abs(depths)));
% hc_reg = squeeze(trapz(depths,(temp_forcing+273.15).* fjord_rho,2)./max(abs(depths)).* fjords_processed(1).p.cw);
sc_reg = squeeze(trapz(depths,salt_forcing,2)./max(abs(depths)));
hc_reg = squeeze(trapz(depths,(temp_forcing+273.15),2)./max(abs(depths)));
taxis_shelf = 1:1:size(sc_reg,1);
tr_ohc_shelf = NaN(size(regions));
tr_osc_shelf = NaN(size(regions));
for i_reg=1:length(regions)
    p = polyfit(taxis_shelf,hc_reg(:,i_reg),1);
    tr_ohc_shelf(i_reg) = p(1)*12;
    % tr_ohc_shelf(i_reg) = mean(hc_reg(end-12:end,i_reg))-mean(hc_reg(1:12,i_reg));
    % p = polyfit(taxis_shelf,sc_reg(:,i_reg),1);
    tr_osc_shelf(i_reg) = p(1)*12;
    % tr_osc_shelf(i_reg) = mean(sc_reg(end-12:end,i_reg))-mean(sc_reg(1:12,i_reg));
end
region_line_color = lines(7);
time_axis_plt = datetime(2010,01,15)+1:1:datetime(2018,12,15);

%% Plot the time series to see how they all behave
figure('Name','time series model outputs','Position',[20 20 1000 500])
region_handles = [];
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis_plt)]);
    osc_reg=NaN([n_runs,length(time_axis_plt)]);
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp)
            [ohc_reg(k_run,:),osc_reg(k_run,:)] = get_active_fjord_contents(ensemble(k_run,i_reg));
        else
            ohc_reg(k_run,:) = NaN;
            osc_reg(k_run,:) = NaN;
        end
    end
    subplot(1,2,1); hold on; box on;
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],ohc_reg-273.15));
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],ohc_reg-273.15));
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;

    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 
    

    subplot(1,2,2); hold on; box on;
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],osc_reg));
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],osc_reg));
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;

    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    hp = plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 
    
    region_handles=[region_handles hp];
end
subplot(1,2,1)
text(0.03,1.03,'(a)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Temperature (^oC m^{-3})');
set(gca,'fontsize',14)
subplot(1,2,2)
text(0.03,1.03,'(b)','fontsize',14,'units','normalized')
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','southeast');
hl.NumColumns=3;
xlabel('Time'); ylabel('Salinity (m^{-3})');
set(gca,'fontsize',14)
% exportgraphics(gcf,[figs_path,'ensemble_series_bootstrapped_n',num2str(n_runs),'.png'],'Resolution',300)


%% Plotting surrogate vs numerical model
figure('Name','Model fit for HC','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on    
    uq_plot(gca,Ynum_ohc{i_reg},Ysur_ohc{i_reg},'+')
    uq_plot(gca,Yind_ohc{i_reg},Yvld_ohc{i_reg},'o','color',[0.8500 0.3250 0.0980])
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    box on; grid on;

    rms = rmse(Yvld_ohc{i_reg},Yind_ohc{i_reg},'omitnan');
    % mdl = fitlm(Ynum_ohc{i_reg},Ysur_ohc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s (n=%d)',letters{i_reg},regions{i_reg},ok_runs(i_reg)),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('$\\epsilon_{LOO}=%0.2f$',sur_model_ohc{i_reg}.Error.LOO),'interpreter','latex','units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% hl = legend('Training dataset','Validation dataset','fontsize',12);
% hl.Position(1)=hl.Position(1)+0.175;
% exportgraphics(gcf,[figs_path,'pce_fit_ohc_n',num2str(n_runs),'.png'],'Resolution',300)

figure('Name','Model fit for SC','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on;
    uq_plot(gca,Ynum_osc{i_reg},Ysur_osc{i_reg},'+')    
    uq_plot(gca,Yind_osc{i_reg},Yvld_osc{i_reg},'o','color',[0.8500 0.3250 0.0980])
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    box on; grid on;

    rms = rmse(Yvld_osc{i_reg},Yind_osc{i_reg},'omitnan');
    % mdl = fitlm(Ynum_ohc{i_reg},Ysur_ohc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s (n=%d)',letters{i_reg},regions{i_reg},ok_runs(i_reg)),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('$\\epsilon_{LOO}=%0.2f$',sur_model_osc{i_reg}.Error.LOO),'interpreter','latex','units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% hl = legend('Training dataset','Validation dataset','fontsize',12);
% hl.Position(1)=hl.Position(1)+0.175;
% exportgraphics(gcf,[figs_path,'pce_fit_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Surrogate model kernel density
figure('Name','Surrogate model kernel density','Position',[40 40 850 500]); hold on;
subplot(2,2,1), hold on; box on; grid on
for i_reg=1:n_regions
    plot(1e-3.*ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    % xline(1e-3.*tr_ohc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
    scatter(1e-3.*tr_ohc_shelf(i_reg),pdf(ohc_ks_eval{i_reg},tr_ohc_shelf(i_reg)),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
ylabel('Probability density');
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
xlim([-1.5e-4 1.5e-4]);
set(gca,'fontsize',14)
subplot(2,2,3), hold on; box on; grid on
for i_reg=1:n_regions
    plot(1e-3.*ohc_x,cdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    % xline(1e-3.*tr_ohc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
    scatter(1e-3.*tr_ohc_shelf(i_reg),cdf(ohc_ks_eval{i_reg},tr_ohc_shelf(i_reg)),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Temperature trend (^oC m^{-3}yr^{-1})',fontsize=14);
ylabel('Cumulative probability density');
text(0.05,0.95,'(c)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-1.5e-4 1.5e-4]);
% handle_plots = [];
subplot(2,2,2), hold on; box on; grid on
for i_reg=1:n_regions
    hp = plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    % xline(tr_osc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
    scatter(tr_osc_shelf(i_reg),pdf(osc_ks_eval{i_reg},tr_osc_shelf(i_reg)),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
    % handle_plots = [handle_plots hp];
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlim([-0.1 0.1])
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
ylabel('Probability density');
set(gca,'fontsize',14)
subplot(2,2,4), hold on; box on; grid on
handle_plots = [];
for i_reg=1:n_regions
    hp = plot(osc_x,cdf(osc_ks_eval{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    % xline(tr_osc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
    scatter(tr_osc_shelf(i_reg),cdf(osc_ks_eval{i_reg},tr_osc_shelf(i_reg)),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
    handle_plots = [handle_plots hp];
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Salinity trend (m^{-3}yr^{-1})',fontsize=14);
ylabel('Cumulative probability density');
text(0.05,0.95,'(d)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-0.1 0.1])
hl = legend(handle_plots,regions,'fontsize',14,'Location','east');
% exportgraphics(gcf,[figs_path,'kssur_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% construct the numerical model kernel density plot (just for comparison)
% figure('Name','Numerical model kernel density','Position',[40 40 850 300]); 
% handle_plots = [];
% subplot(1,2,1), hold on; 
% for i_reg=1:n_regions
%     hp = plot(1e-3.*ohc_x,pdf(ohc_ks{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
%     xline(1e-3.*tr_ohc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
% end
% xlabel('Temperature trend (^oC m^{-3}yr^{-1})',fontsize=14); ylabel('Probability',fontsize=14);  box on
% text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
% set(gca,'fontsize',14)
% xlim([-2 4])
% subplot(1,2,2), hold on; 
% for i_reg=1:n_regions
%     hp = plot(osc_x,pdf(osc_ks{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
%     xline(tr_osc_shelf(i_reg),'linewidth',1,'linestyle','--','color',region_line_color(i_reg,:)); 
%     handle_plots = [handle_plots hp];
% end
% xlabel('Salinity trend (m^{-3}yr^{-1})',fontsize=14); box on
% text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
% set(gca,'fontsize',14)
% xlim([-0.1 0.2])
% hl = legend(handle_plots,regions_lbl,'fontsize',14,'Location','northeast');
%exportgraphics(gcf,[figs_path,'ksnum_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%TODO: add histogram of differences?
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

%% Convergence test to see if our choice of n_runs was enough

figure('Name','Convergence test for n_runs','Position',[40 40 850 300]); hold on;
region_handles = [];
subplot(1,2,1), hold on; box on; grid on
for i_reg=1:n_regions
    plot(x_subsample,Yconv_ohc(:,i_reg),'linewidth',2,'Color',region_line_color(i_reg,:));
    % xline(ok_runs(i_reg),'linewidth',1,'linestyle','--','Color',region_line_color(i_reg,:));
    [d,i_conv] = min(abs(double(ok_runs(i_reg))-double(x_subsample)));
    scatter(ok_runs(i_reg),Yconv_ohc(i_conv,i_reg),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
end
ylabel('Avg. temperature trend (^oC m^{-3}yr^{-1})',fontsize=14); xlabel('experimental design size (n)','fontsize',14);  
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([0 n_runs])
subplot(1,2,2), hold on; box on; grid on
for i_reg=1:n_regions
    hp = plot(x_subsample,Yconv_osc(:,i_reg),'linewidth',2,'Color',region_line_color(i_reg,:)); 
    % xline(ok_runs(i_reg),'linewidth',1,'linestyle','--','Color',region_line_color(i_reg,:));
    [d,i_conv] = min(abs(double(ok_runs(i_reg))-double(x_subsample)));
    scatter(ok_runs(i_reg),Yconv_osc(i_conv,i_reg),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
    region_handles = [region_handles hp];
end
ylabel('Avg. salinity trend (m^{-3}yr^{-1})','fontsize',14); xlabel('experimental design size (n)','fontsize',14);   
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([0 n_runs])
hl = legend(region_handles,regions,'fontsize',14,'Location','southwest');
hl.NumColumns=2;
% exportgraphics(gcf,[figs_path,'nruns_convergence_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Box plots showing results of bootstrapping

% figure('Name','PCE accuracy test','Position',[40 40 850 300]); hold on;
% for i_reg=1:n_regions
%     ohc_boostrap = Yboo_ohc{i_reg};
%     osc_boostrap = Yboo_osc{i_reg};
%     subplot(1,2,1), hold on; box on
%     boxplot(i_reg,ohc_boostrap(:)')
% 
%     subplot(1,2,2), hold on; box on
%     boxplot(i_reg,osc_boostrap(:)')
% end