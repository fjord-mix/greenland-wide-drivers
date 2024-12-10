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