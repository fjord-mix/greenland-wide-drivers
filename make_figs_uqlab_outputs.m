%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];
end

% get the range of results for computing the probability distributions
ohc_x = linspace(0.99*min(ohc_out(:)),1.01*max(ohc_out(:)),1000);
osc_x = linspace(0.99*min(osc_out(:)),1.01*max(osc_out(:)),1000);

%% Plot the time series to see how they all behave
region_line_color = lines(7);
time_axis = datetime(2010,01,15)+1:1:datetime(2018,12,15);
region_handles = [];
figure('Name','time series model outputs','Position',[20 20 1000 500])
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis)]);
    osc_reg=NaN([n_runs,length(time_axis)]);
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).ohc)
            % ohc_reg(k_run,:) = ensemble(k_run,i_reg).ohc; 
            % compute quantities per unit volume
            fjord_run = ensemble(k_run,i_reg);
            ohc_reg(k_run,:) = sum(fjord_run.ohc,1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            osc_reg(k_run,:) = sum(fjord_run.osc,1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            % above sill only
            % ohc_reg(k_run,:) = sum(fjord_run.ohc(1:fjord_run.p.N,:),1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            % osc_reg(k_run,:) = sum(fjord_run.osc(1:fjord_run.p.N,:),1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            % below sill only
            % ohc_reg(k_run,:) = fjord_run.ohc(end,:)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            % osc_reg(k_run,:) = fjord_run.osc(end,:)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            
            
        else
            ohc_reg(k_run,:) = NaN;
            osc_reg(k_run,:) = NaN;
        end
    end
    subplot(1,2,1); hold on; box on;
    % multiplying by 1e-3 to change units to kJ
    lower_bnd = 1e-3.*min(ohc_reg,[],'omitnan');
    upper_bnd = 1e-3.*max(ohc_reg,[],'omitnan');
    % lower_bnd = 1e-3.*prctile(ohc_reg,25,1);
    % upper_bnd = 1e-3.*prctile(ohc_reg,75,1);
    
    median_ln = 1e-3.*median(ohc_reg,1,'omitnan');
    x2 = [time_axis, fliplr(time_axis)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    hp = plot(time_axis,median_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 

    subplot(1,2,2); hold on; box on;
    lower_bnd = min(osc_reg,[],'omitnan');
    upper_bnd = max(osc_reg,[],'omitnan');
    % lower_bnd = prctile(osc_reg,25,1);
    % upper_bnd = prctile(osc_reg,75,1);
    
    median_ln = median(osc_reg,1,'omitnan');
    x2 = [time_axis, fliplr(time_axis)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    hp = plot(time_axis,median_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 
    region_handles=[region_handles hp];
end
subplot(1,2,1)
text(0.03,1.03,'(a)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Heat content (kJ m^{-3})');
set(gca,'fontsize',14)
subplot(1,2,2)
text(0.03,1.03,'(b)','fontsize',14,'units','normalized')
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','northeast');
hl.NumColumns=3;
xlabel('Time'); ylabel('Salt content (g m^{-3})');
set(gca,'fontsize',14)
exportgraphics(gcf,[figs_path,'ensemble_series_totalspread_n',num2str(n_runs),'.png'],'Resolution',300)


%% Plotting surrogate vs numerical model
figure('Name','Model fit for HC','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on    
    uq_plot(gca,Ynum_ohc{i_reg},Ysur_ohc{i_reg},'+')
    uq_plot(gca,Yind_ohc{i_reg},Yvld_ohc{i_reg},'o','color',[0.8500 0.3250 0.0980])
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    box on; grid on;

    % rms = rmse(Ysur_ohc{i_reg},Ynum_ohc{i_reg},'omitnan');
    mdl = fitlm(Yind_ohc{i_reg},Yvld_ohc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('RMSE=%0.2f J m^{-3}',rms),'units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
hl = legend(hb,{'Training dataset','Validation dataset'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
exportgraphics(gcf,[figs_path,'pce_fit_ohc_n',num2str(n_runs),'.png'],'Resolution',300)

figure('Name','Model fit for SC','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on;
    uq_plot(gca,Ynum_osc{i_reg},Ysur_osc{i_reg},'+')    
    uq_plot(gca,Yind_osc{i_reg},Yvld_osc{i_reg},'o','color',[0.8500 0.3250 0.0980])
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    box on; grid on;

    % rms = rmse(Ysur_osc{i_reg},Ynum_osc{i_reg},'omitnan');
    mdl = fitlm(Yind_osc{i_reg},Yvld_osc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('RMSE=%0.2f J m^{-3}',rms),'units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
hl = legend(hb,{'Training dataset','Validation dataset'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
exportgraphics(gcf,[figs_path,'pce_fit_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting the resulting kernel density

figure('Name','Surrogate model kernel density','Position',[40 40 850 300]); hold on;
subplot(1,2,1), hold on; box on
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content change (J m^{-3})',fontsize=14); ylabel('Kernel density',fontsize=14);  
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-200 100])
subplot(1,2,2), hold on; box on
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content change (g m^{-3})',fontsize=14); 
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-7e-3 4e-3])
hl = legend(regions,'fontsize',14,'Location','west');
exportgraphics(gcf,[figs_path,'kssur_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

% construct the numerical model kernel density plot (just for comparison
% figure('Name','Numerical model kernel density','Position',[40 40 850 300]); 
% subplot(1,2,1), hold on; 
% for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks_as{i_reg},ohc_x),'linewidth',2); end
% xlabel('Heat content change (J m^{-3})',fontsize=14); ylabel('Probability',fontsize=14);  box on
% text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
% set(gca,'fontsize',14)
% subplot(1,2,2), hold on; 
% for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks_as{i_reg},osc_x),'linewidth',2); end
% xlabel('Salt content change (g m^{-3})',fontsize=14); box on
% text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
% set(gca,'fontsize',14)
% hl = legend(regions_lbl,'fontsize',14,'Location','west');
% %exportgraphics(gcf,[figs_path,'ksnum_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)


%% Plotting the Sobol indices

figure('Name','First-order Sobol Indices','position',[40 40 1000 400])
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    sobolResults_osc  = sobolA_osc{i_reg}.Results;
    sobolIndices = [sobolResults_ohc.FirstOrder sobolResults_osc.FirstOrder];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts.Marginals), sobolIndices, 1.,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('First-order Sobol indices'); end
    ylim([0 1])
end
hl = legend(hb,{'Heat content','Salt content'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
exportgraphics(gcf,[figs_path,'sobol_first_n',num2str(n_runs),'.png'],'Resolution',300)


figure('Name','Total Sobol Indices','position',[40 40 1000 400])
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    sobolResults_osc  = sobolA_osc{i_reg}.Results;
    sobolIndices = [sobolResults_ohc.Total sobolResults_osc.Total];

    subplot(2,4,i_reg); hold on; box on; grid on
    hb=uq_bar(gca,1:length(IOpts.Marginals), sobolIndices, 1.,'grouped');
    text(0.01,1.075,sprintf('(%s) %s',letters{i_reg},regions{i_reg}),'units','normalized','fontsize',12)
    set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
    if i_reg < 4, xlabel(''); end
    if i_reg==1 || i_reg==5, ylabel('Total Sobol indices'); end
    ylim([0 1])
end
hl = legend(hb,{'Heat content','Salt content'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;
exportgraphics(gcf,[figs_path,'sobol_total_n',num2str(n_runs),'.png'],'Resolution',300)

%% Convergence test to see if our choice of n_runs was enough

figure('Name','Convergence test for n_runs','Position',[40 40 850 300]); hold on;
subplot(1,2,1), hold on; box on
for i_reg=1:n_regions,plot(x_subsample,Yconv_ohc(:,i_reg),'linewidth',2); end
ylabel('Avg. heat content change (J m^{-3})',fontsize=14); xlabel('experimental design size (n)','fontsize',14);  
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
% xlim([-200 100])
subplot(1,2,2), hold on; box on
for i_reg=1:n_regions,plot(x_subsample,Yconv_osc(:,i_reg),'linewidth',2); end
ylabel('Avg. salt content change (g m^{-3})','fontsize',14); xlabel('experimental design size (n)','fontsize',14);   
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
% xlim([-7e-3 4e-3])
hl = legend(regions,'fontsize',14,'Location','southeast');
hl.NumColumns=2;
exportgraphics(gcf,[figs_path,'nruns_convergence_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Box plots showing results of bootstrapping

figure('Name','PCE accuracy test','Position',[40 40 850 300]); hold on;
for i_reg=1:n_regions
    ohc_boostrap = Yboo_ohc{i_reg};
    osc_boostrap = Yboo_osc{i_reg};
    subplot(1,2,1), hold on; box on
    boxplot(i_reg,ohc_boostrap(:)')

    subplot(1,2,2), hold on; box on
    boxplot(i_reg,osc_boostrap(:)')
end