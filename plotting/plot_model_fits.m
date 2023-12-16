function hf = plot_model_fits(Ynum,Ysur,Yind,Yvld,ok_runs,ok_vruns,var,units)

letters = {'a','b','c','d','e','f','g','h'};
regions = {'SW','SE','CW','CE','NW','NE','NO'};
n_regions = length(regions);

hf = figure('Name',['Model fit for ',var],'position',[40 40 1000 400]);
for i_reg=1:n_regions
    subplot(2,4,i_reg); hold on    
    uq_plot(gca,Ynum{i_reg},Ysur{i_reg},'+')
    if ~isempty(Yvld) && ~isempty(Yind)
        uq_plot(gca,Yind{i_reg},Yvld{i_reg},'o','color',[0.8500 0.3250 0.0980])
        rms = rmse(Yvld{i_reg},Yind{i_reg},'omitnan');
        text(0.98,0.09,sprintf('RMSE=%0.2f %s',rms,units),'units','normalized','horizontalAlignment','right','fontsize',14)
    end
    hl = refline(1,0); hl.LineStyle='--'; hl.Color='r';
    box on; grid on;

    % mdl = fitlm(Yind_ohc{i_reg},Yvld_ohc{i_reg});
    text(0.05,0.95,sprintf('(%s) %s (n=%d)',letters{i_reg},regions{i_reg},ok_runs(i_reg)),'units','normalized','fontsize',14)
    % text(0.98,0.07,sprintf('$\\epsilon_{LOO}=%0.2f$',sur_model_ohc{i_reg}.Error.LOO),'interpreter','latex','units','normalized','horizontalAlignment','right','fontsize',14)
    % set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% hl = legend('Training dataset','Validation dataset','fontsize',12);
% hl.Position(1)=hl.Position(1)+0.175;
end