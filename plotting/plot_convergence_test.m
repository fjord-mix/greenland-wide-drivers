function hf = plot_convergence_test(x_subsample,Yconv_ohc,Yconv_osc,ok_runs)
n_regions = length(Yconv_ohc);
region_line_color = lines(n_regions);

hf = figure('Name','Convergence test for n_runs','Position',[40 40 850 300]); hold on;
region_handles = [];
subplot(1,2,1), hold on; box on; grid on
for i_reg=1:n_regions
    plot(x_subsample,Yconv_ohc(:,i_reg),'linewidth',2,'Color',region_line_color(i_reg,:));
    % xline(ok_runs(i_reg),'linewidth',1,'linestyle','--','Color',region_line_color(i_reg,:));
    [d,i_conv] = min(abs(double(ok_runs(i_reg))-double(x_subsample)));
    scatter(ok_runs(i_reg),Yconv_ohc(i_conv,i_reg),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
end
ylabel('Avg. temperature difference (^oC)',fontsize=14); xlabel('experimental design size (n)','fontsize',14);  
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([0 n_runs])
subplot(1,2,2), hold on; box on; grid on
for i_reg=1:n_regions
    hp = plot(x_subsample,Yconv_osc(:,i_reg),'linewidth',2,'Color',region_line_color(i_reg,:)); 
    % xline(ok_runs(i_reg),'linewidth',1,'linestyle','--','Color',region_line_color(i_reg,:));
    [d,i_conv] = min(abs(double(ok_runs(i_reg))-double(x_subsample)));
    scatter(ok_runs(i_reg),Yconv_osc(i_conv,i_reg),40,'filled','o','MarkerFaceColor',region_line_color(i_reg,:));
    region_handles = [region_handles hp];
end
ylabel('Avg. salinity difference','fontsize',14); xlabel('experimental design size (n)','fontsize',14);   
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([0 n_runs])
hl = legend(region_handles,regions,'fontsize',14,'Location','southwest');
hl.NumColumns=2;

end