function hf = plot_borgonovo_bars(BorgonovoA_ohc,BorgonovoA_osc,IOpts,letters,regions)
hf = figure('Name','Borgonovo Indices','position',[40 40 1000 400]);
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
    % ylim([0 0.35])
end
hl = legend(hb,{'Temperature','Salinity'},'fontsize',12);
hl.Position(1)=hl.Position(1)+0.175;

end