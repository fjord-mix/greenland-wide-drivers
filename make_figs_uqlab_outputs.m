%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];
end

%% Plot the time series to see how they all behave
region_line_color = lines(7);
time_axis = datetime(2010,01,15)+1:1:datetime(2018,12,15);
region_handles = [];
figure('Name','time series model outputs','Position',[20 20 600 400])
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis)]);
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).ohc)
            ohc_reg(k_run,:) = ensemble(k_run,i_reg).ohc; 
        else
            ohc_reg(k_run,:) = NaN;
        end
    end
    % multiplying by 1e-3 to change units to kJ
    % lower_bnd = 1e-3.*min(ohc_reg,[],'omitnan');
    % upper_bnd = 1e-3.*max(ohc_reg,[],'omitnan');
    lower_bnd = 1e-3.*prctile(ohc_reg,25,1);
    upper_bnd = 1e-3.*prctile(ohc_reg,75,1);
    
    median_ln = 1e-3.*median(ohc_reg,1,'omitnan');
    x2 = [time_axis, fliplr(time_axis)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(time_axis,median_ln,'Color',region_line_color(i_reg,:)); 
    region_handles=[region_handles hp];
end
hl = legend(region_handles,regions_lbl,'fontsize',14,'Location','northeastoutside');
% hl.NumColumns=3;
xlabel('Time'); ylabel('Heat content (kJ m^{-3})');
set(gca,'fontsize',14)
exportgraphics(gcf,[figs_path,'ensemble_series_ohc_n',num2str(n_runs),'.png'],'Resolution',300)

%% construct the numerical model kernel density plot
ohc_x = linspace(0.99*min(ohc_out(:)),1.01*max(ohc_out(:)),1000);
osc_x = linspace(0.99*min(osc_out(:)),1.01*max(osc_out(:)),1000);

figure('Name','Numerical model kernel density','Position',[40 40 850 300]); 
subplot(1,2,1), hold on; 
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content change (J m^{-3})',fontsize=14); ylabel('Probability',fontsize=14);  box on
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
subplot(1,2,2), hold on; 
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content change (g m^{-3})',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
hl = legend(regions_lbl,'fontsize',14,'Location','west');
exportgraphics(gcf,[figs_path,'ksnum_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)



%% Plotting surrogate vs numerical model
figure('Name','Model fit','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    uq_plot(gca,Ynum_ohc{i_reg},Ysur_ohc{i_reg},'+')
    box on; grid on;

    rms = rmse(Ysur_ohc{i_reg},Ynum_ohc{i_reg},'omitnan');
    mdl = fitlm(Ynum_ohc{i_reg},Ysur_ohc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('RMSE=%0.2f J m^{-3}',rms),'units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% savefig([figs_path,'pce_fit_indices_n30.png'],'Resolution',300)
exportgraphics(gcf,[figs_path,'pce_fit_ohc_n',num2str(n_runs),'.png'],'Resolution',300)

figure('Name','Model fit','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on;
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    uq_plot(gca,Ynum_osc{i_reg},Ysur_osc{i_reg},'+')    
    box on; grid on;

    rms = rmse(Ysur_osc{i_reg},Ynum_osc{i_reg},'omitnan');
    mdl = fitlm(Ynum_osc{i_reg},Ysur_osc{i_reg});

    text(0.05,0.95,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('RMSE=%0.2f J m^{-3}',rms),'units','normalized','horizontalAlignment','right','fontsize',14)
    text(0.98,0.09,sprintf('R^2=%0.2f',mdl.Rsquared.Adjusted),'units','normalized','horizontalAlignment','right','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% savefig([figs_path,'pce_fit_indices_n30.png'],'Resolution',300)
exportgraphics(gcf,[figs_path,'pce_fit_osc_n',num2str(n_runs),'.png'],'Resolution',300)

%% Plotting the surrogate model kernel density
figure('Name','Surrogate model kernel density','Position',[40 40 850 300]); hold on;
subplot(1,2,1), hold on; box on
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content change (J m^{-3})',fontsize=14); ylabel('Kernel density',fontsize=14);  
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-180 100])
subplot(1,2,2), hold on; box on
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content change (g m^{-3})',fontsize=14); 
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-6e-3 5e-3])
hl = legend(regions,'fontsize',14,'Location','west');
exportgraphics(gcf,[figs_path,'kssur_ohc_osc_n',num2str(n_runs),'.png'],'Resolution',300)


%% Plotting the Sobol indices

figure('Name','Sobol Indices','position',[40 40 1000 400])
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    sobolResults_osc  = sobolA_osc{i_reg}.Results;
    % sobolIndices_ohc = [sobolResults_ohc.Total-sobolResults_ohc.FirstOrder sobolResults_ohc.FirstOrder];
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

% % Plotting indices grouped per region
% sobolTotal_ohc = [];
% sobolFirstOrder_ohc = [];
% sobolTotal_osc = [];
% sobolFirstOrder_osc = [];
% 
% % Gather up all regions here
% for i_reg=1:7
%     sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
%     % uq_print(sobolAnalysis)
%     sobolTotal_ohc      = [sobolTotal_ohc sobolResults_ohc.Total];
%     sobolFirstOrder_ohc = [sobolFirstOrder_ohc sobolResults_ohc.FirstOrder];
% 
%     sobolResults_osc  = sobolA_osc{i_reg}.Results;    
%     sobolTotal_osc      = [sobolTotal_osc sobolResults_osc.Total];
%     sobolFirstOrder_osc = [sobolFirstOrder_osc sobolResults_osc.FirstOrder];
% end

% uq_figure('Name', 'Total Sobol'' Indices','Position',[40 40 1000 800])
% barWidth = 1.;
% subplot(2,1,1)
% uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_ohc, barWidth)
% xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
% ylabel('Total Sobol indices'); % xlabel('Variable'); 
% text(0.01,1.075,'(a) Heat content change','units','normalized','fontsize',20)
% set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
% 
% subplot(2,1,2)
% uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_osc, barWidth)
% xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
% xlabel('Variable'); ylabel('Total Sobol indices')
% text(0.01,1.075,'(b) Salt content change','units','normalized','fontsize',20)
% set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
% legend(regions,'Location', 'northeastoutside')
% % exportgraphics(gcf,[figs_path,'test_total_sobol_indices_n30.png'],'Resolution',300)
% 
% uq_figure('Name', 'First-Order Sobol'' Indices','Position',[40 40 1000 500])
% subplot(1,2,1)
% uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_ohc, barWidth)
% xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
% xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OHC)')
% set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
% subplot(1,2,2)
% uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_osc, barWidth)
% xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
% xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OSC)')
% set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
% uq_legend(regions,'Location', 'northeastoutside')
% exportgraphics(gcf,[figs_path,'test_1st_sobol_indices_n50.png'],'Resolution',300)