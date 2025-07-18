function [hf,hp,hl] = plot_surrogate_model_results(ohc_out,osc_out,ohc_ks_eval,osc_ks_eval)

% get the range of results for showing the probability distributions
ohc_x = linspace(1.2*min(ohc_out(:)),1.2*max(ohc_out(:)),1000);
osc_x = linspace(1.2*min(osc_out(:)),1.2*max(osc_out(:)),1000);


regions = {'SW','SE','CW','CE','NW','NE','NO'};
n_regions=length(regions);
region_line_color = lines(n_regions);

hf = figure('Name','Surrogate model kernel density','Position',[40 40 850 500]); hold on;
subplot(2,2,1), hold on; box on; grid on
for i_reg=1:n_regions
    plot(ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
ylabel('Probability density');
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
xlim([-1. 1]);
set(gca,'fontsize',14)
subplot(2,2,3), hold on; box on; grid on
for i_reg=1:n_regions
    plot(ohc_x,cdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Mean temperature difference (^oC)',fontsize=14);
ylabel('Cumulative probability density');
text(0.05,0.95,'(c)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-1. 1]);
subplot(2,2,2), hold on; box on; grid on
for i_reg=1:n_regions
    hp = plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlim([-0.35 0.35])
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
ylabel('Probability density');
set(gca,'fontsize',14)
subplot(2,2,4), hold on; box on; grid on
handle_plots = [];
for i_reg=1:n_regions
    hp = plot(osc_x,cdf(osc_ks_eval{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    handle_plots = [handle_plots hp];
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Mean salinity difference',fontsize=14);
ylabel('Cumulative probability density');
text(0.05,0.95,'(d)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-0.35 0.35])
hl = legend(handle_plots,regions,'fontsize',14,'Location','west');
end