%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];
end

% construct the plot itself
ohc_x = linspace(0.99*min(ohc_out(:)),1.01*max(ohc_out(:)),1000);
osc_x = linspace(0.99*min(osc_out(:)),1.01*max(osc_out(:)),1000);
figure('Position',[40 40 850 300]); 
subplot(1,2,1), hold on; 
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content (J m^{-3})',fontsize=14); ylabel('Probability',fontsize=14);  box on
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
subplot(1,2,2), hold on; 
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content (g m^{-3})',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
hl = legend(regions_lbl,'fontsize',14,'Location','north');
% exportgraphics(gcf,[figs_path,'test_output_ohc_osc_ks_n30.png'],'Resolution',300)


% figure; hold on; for k=1:n_runs, plot(ensemble(k).time,ensemble(k).ohc); end

%% Plotting surrogate vs numerical model
figure('Name','Model fit','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg)
    uq_plot(gca,Ynum_ohc{i_reg},Ysur_ohc{i_reg},'+')
    text(0.05,0.95,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% savefig([figs_path,'pce_fit_indices_n30.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'pce_fit_indices_n30.png'],'Resolution',300)

%% Plotting the surrogate model kernel density
figure('Position',[40 40 1000 400]); hold on;
subplot(1,2,1), hold on; box on
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content (J m^{-3})',fontsize=14); ylabel('Kernel density',fontsize=14);  
set(gca,'fontsize',14)
subplot(1,2,2), hold on; box on
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content (g m^{-3})',fontsize=14); 
set(gca,'fontsize',14)
hl = legend(regions,'fontsize',14,'Location','northeast');
% exportgraphics(gcf,[figs_path,'test_pce_ks_n1e6_n30.png'],'Resolution',300)


%% Plotting the Sobol indices
sobolTotal_ohc = [];
sobolFirstOrder_ohc = [];
sobolTotal_osc = [];
sobolFirstOrder_osc = [];

% Gather up all regions here
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    % uq_print(sobolAnalysis)
    sobolTotal_ohc      = [sobolTotal_ohc sobolResults_ohc.Total];
    sobolFirstOrder_ohc = [sobolFirstOrder_ohc sobolResults_ohc.FirstOrder];

    sobolResults_osc  = sobolA_osc{i_reg}.Results;    
    sobolTotal_osc      = [sobolTotal_osc sobolResults_osc.Total];
    sobolFirstOrder_osc = [sobolFirstOrder_osc sobolResults_osc.FirstOrder];
end

% Plotting indices
uq_figure('Name', 'Total Sobol'' Indices','Position',[40 40 1000 500])
barWidth = 1.;
subplot(1,2,1)
uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_ohc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('Total Sobol indices (OHC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)

subplot(1,2,2)
uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_osc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('Total Sobol indices (OSC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
uq_legend(regions,'Location', 'northeastoutside')
% exportgraphics(gcf,[figs_path,'test_total_sobol_indices_n30.png'],'Resolution',300)

uq_figure('Name', 'First-Order Sobol'' Indices','Position',[40 40 1000 500])
subplot(1,2,1)
uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_ohc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OHC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
subplot(1,2,2)
uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_osc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OSC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
uq_legend(regions,'Location', 'northeastoutside')
% exportgraphics(gcf,[figs_path,'test_1st_sobol_indices_n50.png'],'Resolution',300)