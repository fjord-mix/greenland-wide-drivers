function hf = plot_ensemble_dt_ds(ensemble,time_axis_plt,regions_lbl)

n_regions=size(ensemble,2);
n_runs=size(ensemble,1);
region_handles = [];
region_line_color = lines(n_regions);

avg_silldepth   = NaN([1,size(ensemble,2)]);
avg_activedepth = NaN([1,size(ensemble,2)]);

hf = figure('Name','time series model outputs','Position',[20 20 1200 800]);
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis_plt)]);
    osc_reg=NaN([n_runs,length(time_axis_plt)]);
    ohc_shf=NaN([n_runs,length(time_axis_plt)]);
    osc_shf=NaN([n_runs,length(time_axis_plt)]);
    ohc_fjd=NaN([n_runs,length(time_axis_plt)]);
    osc_fjd=NaN([n_runs,length(time_axis_plt)]);
    hact_fjd=NaN(size(time_axis_plt));
    zsill_fjd=NaN([1,n_runs]);

    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp)
            [ohc_fjd(k_run,:),osc_fjd(k_run,:),hact_fjd(k_run)] = get_active_fjord_contents(ensemble(k_run,i_reg));

            zs0 = unique(sort([0,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.silldepth]));
            i_sill = find(zs0 == ensemble(k_run,i_reg).p.silldepth);
            zs0 = zs0(i_sill:end);

            Ss0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs0,'pchip','extrap');
            Ts0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs0,'pchip','extrap');

            ohc_shf(k_run,:) = squeeze(trapz(zs0,Ts0)./abs(ensemble(k_run,i_reg).p.silldepth));
            osc_shf(k_run,:) = squeeze(trapz(zs0,Ss0)./abs(ensemble(k_run,i_reg).p.silldepth));

            ohc_reg(k_run,:) = ohc_fjd(k_run,:) - ohc_shf(k_run,:);
            osc_reg(k_run,:) = osc_fjd(k_run,:) - osc_shf(k_run,:);

            zsill_fjd(k_run) = abs(ensemble(k_run,i_reg).p.silldepth);
        else
            ohc_reg(k_run,:) = NaN;
            osc_reg(k_run,:) = NaN;
        end
        avg_activedepth(i_reg) = mean(hact_fjd,'omitnan');
        avg_silldepth(i_reg)   = mean(zsill_fjd,'omitnan');
    end
    subplot(2,2,1); hold on; box on;
    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],ohc_reg));
    std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_reg);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1.5); 
    

    subplot(2,2,2); hold on; box on;
    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],osc_reg));
    std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],osc_reg);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp = plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1.5); 
    
    region_handles=[region_handles hp];

    subplot(2,2,3); hold on; box on;
    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],ohc_fjd));
    std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_fjd);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    % fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp_fjd = plot(time_axis_plt,mean_ln,'Color','k','linewidth',1.); 
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1.); 

    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],ohc_shf));
    hp_shf = plot(time_axis_plt,mean_ln,'Color','k','linewidth',1,'linestyle','--'); 
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 
    

    subplot(2,2,4); hold on; box on;
    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],osc_fjd));
    std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],osc_fjd);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    % fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp = plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1); 

    mean_ln = mean(bootstrp(1e3,@(x)[mean(x,1,'omitnan')],osc_shf));
    plot(time_axis_plt,mean_ln,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 
end
subplot(2,2,1)
text(0.03,1.07,'(a)','fontsize',14,'units','normalized')
ylabel('Temperature difference (^oC)'); % xlabel('Time');
set(gca,'fontsize',14)
subplot(2,2,2)
text(0.03,1.07,'(b)','fontsize',14,'units','normalized')
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','southeast');
hl.NumColumns=3;
ylabel('Salinity difference'); % xlabel('Time'); 
set(gca,'fontsize',14)

subplot(2,2,3)
text(0.03,1.07,'(c)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Temperature (^oC)');
hl2 = legend([hp_fjd hp_shf],{'fjord','shelf'},'fontsize',10,'Location','northeast');
set(gca,'fontsize',14)
subplot(2,2,4)
text(0.03,1.07,'(d)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Salinity');
set(gca,'fontsize',14)


end