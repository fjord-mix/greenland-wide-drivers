function [hf,hp,hl] = plot_surrogate_model_results(ohc_x,osc_x,ohc_ks_eval,osc_ks_eval)

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
xlim([-2 1]);
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
xlim([-2 1]);
subplot(2,2,2), hold on; box on; grid on
for i_reg=1:n_regions
    hp = plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlim([-2.5 0.5])
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
xlim([-2.5 0.5])
hl = legend(handle_plots,regions,'fontsize',14,'Location','west');
end