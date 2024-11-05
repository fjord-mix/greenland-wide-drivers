function [hf,hp,hl] = plot_numerical_model_distributions(ohc_out,osc_out,regions_lbl)

n_regions = length(regions_lbl);
region_line_color = lines(n_regions);

%% Computing the distributions from the model outputs
% get the range of results for showing the probability distributions
ohc_x = linspace(1.2*min(ohc_out(:)),1.2*max(ohc_out(:)),1000);
osc_x = linspace(1.2*min(osc_out(:)),1.2*max(osc_out(:)),1000);

ohc_ks  = cell([1, n_regions]);
osc_ks  = cell([1, n_regions]);
for i_reg=1:n_regions
    % ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
    % osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
    ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
    osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
end

%% Plotting
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
xlim([-0.5 1])
subplot(1,2,2), hold on; box on; grid on;
for i_reg=1:n_regions
    hp = plot(osc_x,pdf(osc_ks{i_reg},osc_x),'linewidth',2,'color',region_line_color(i_reg,:)); 
    handle_plots = [handle_plots hp];
end
xline(0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
xlabel('Salinity difference',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([-0.3 0.2])
hl = legend(handle_plots,regions_lbl,'fontsize',14,'Location','west');
end