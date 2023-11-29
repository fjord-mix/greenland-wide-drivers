function hf = plot_ensemble_dt_ds(ensemble,time_axis_plt,regions_lbl)

n_regions=size(ensemble,2);
n_runs=size(ensemble,1);
region_handles = [];
region_line_color = lines(n_regions);

hf = figure('Name','time series model outputs','Position',[20 20 1000 500]);
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis_plt)]);
    osc_reg=NaN([n_runs,length(time_axis_plt)]);
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp)
            [ohc_fjord,osc_fjord] = get_active_fjord_contents(ensemble(k_run,i_reg));
            ohc_shelf = squeeze(trapz(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts)./max(abs(ensemble(k_run,i_reg).zs)));
            osc_shelf = squeeze(trapz(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss)./max(abs(ensemble(k_run,i_reg).zs)));

            ohc_reg(k_run,:) = ohc_fjord - ohc_shelf;
            osc_reg(k_run,:) = osc_fjord - osc_shelf;
        else
            ohc_reg(k_run,:) = NaN;
            osc_reg(k_run,:) = NaN;
        end
    end
    subplot(1,2,1); hold on; box on;
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],ohc_reg));
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],ohc_reg));
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 
    

    subplot(1,2,2); hold on; box on;
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],osc_reg));
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],osc_reg));
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    hp = plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',2); 
    
    region_handles=[region_handles hp];
end
subplot(1,2,1)
text(0.03,1.03,'(a)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Temperature difference (^oC)');
set(gca,'fontsize',14)
subplot(1,2,2)
text(0.03,1.03,'(b)','fontsize',14,'units','normalized')
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','southeast');
hl.NumColumns=3;
xlabel('Time'); ylabel('Salinity difference');
set(gca,'fontsize',14)

end