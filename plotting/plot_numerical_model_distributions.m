function [hf,hp,hl] = plot_numerical_model_distributions(ohc_x,osc_x,ohc_ks,osc_ks,regions_lbl)

n_regions = length(regions_lbl);
region_line_color = lines(n_regions);

hf = figure('Name','Numerical model kernel density','Position',[40 40 850 300]); 
handle_plots = [];
subplot(1,2,1), hold on; box on; grid on;
for i_reg=1:n_regions
    hp = plot(ohc_x,pdf(ohc_ks{i_reg},ohc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Temperature difference (^oC)',fontsize=14); ylabel('Probability',fontsize=14);  box on
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-2 1])
subplot(1,2,2), hold on; box on; grid on;
for i_reg=1:n_regions
    hp = plot(osc_x,pdf(osc_ks{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    handle_plots = [handle_plots hp];
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Salinity difference',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-2.5 0.5])
hl = legend(handle_plots,regions_lbl,'fontsize',14,'Location','west');
end