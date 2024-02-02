function hf = plot_seasonal_cycle(tt_ensemble,regions_lbl,verbose)

if nargin < 3, verbose=0; end

n_regions=length(regions_lbl);
region_handles = [];
sf_handles = [];
region_line_color = lines(n_regions);

hf = figure('Name','Seasonal cycles','Position',[20 20 1000 600]);

for i_reg=1:n_regions
    clim_mon = groupsummary(tt_ensemble{i_reg},"Time","monthofyear","mean");

    subplot(2,2,1); hold on; grid on; box on;
    hp_r = plot(clim_mon.mean_Tf,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    plot(clim_mon.mean_Ts,'LineStyle','--','Color',region_line_color(i_reg,:),'linewidth',1.5);

    subplot(2,2,2); hold on; grid on; box on;
    hp_f = plot(clim_mon.mean_Sf,'-k','linewidth',1.);
    hp_s = plot(clim_mon.mean_Ss,'--k','linewidth',1.);
    plot(clim_mon.mean_Sf,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    plot(clim_mon.mean_Ss,'LineStyle','--','Color',region_line_color(i_reg,:),'linewidth',1.5);

    subplot(2,2,3); hold on; grid on; box on;
    plot(clim_mon.mean_dT,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    
    subplot(2,2,4); hold on; grid on; box on;
    plot(clim_mon.mean_dS,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);

    region_handles = [region_handles hp_r];
    sf_handles = [hp_f hp_s];

    if verbose
        [~,imax_tf] = max(clim_mon.mean_Tf);
        [~,imin_tf] = min(clim_mon.mean_Tf);
        [~,imax_ts] = max(clim_mon.mean_Ts);
        [~,imin_ts] = min(clim_mon.mean_Ts);
        [~,imax_sf] = max(clim_mon.mean_Sf);
        [~,imin_sf] = min(clim_mon.mean_Sf);
        [~,imax_ss] = max(clim_mon.mean_Ss);
        [~,imin_ss] = min(clim_mon.mean_Ss);
        disp('======================================')
        fprintf('region: %s\n',regions_lbl{i_reg})
        fprintf('Tf max in %d and min in %d\n',imax_tf,imin_tf)
        fprintf('Ts max in %d and min in %d\n',imax_ts,imin_ts)
        disp('')
        fprintf('Sf max in %d and min in %d\n',imax_sf,imin_sf)
        fprintf('Ss max in %d and min in %d\n',imax_ss,imin_ss)
    end
end

subplot(2,2,1)
text(0.02,1.05,'(a)','fontsize',14,'Units','normalized')
ylabel('Temperature (^oC)'); 
hl = legend(region_handles,regions_lbl,'location','northwest','fontsize',12);
hl.NumColumns = 2;
xlim([1 12])
set(gca,'fontsize',14)

subplot(2,2,2)
text(0.02,1.05,'(b)','fontsize',14,'Units','normalized')
ylabel('Salinity');
hl = legend(sf_handles,{'Fjord','Shelf'},'fontsize',12);
hl.NumColumns = 2;
xlim([1 12])
set(gca,'fontsize',14)

subplot(2,2,3)
text(0.02,1.05,'(c)','fontsize',14,'Units','normalized')
ylabel('Temperature difference (^oC)'); 
xlabel('Month'); xlim([1 12])
set(gca,'fontsize',14)

subplot(2,2,4)
text(0.02,1.05,'(d)','fontsize',14,'Units','normalized')
ylabel('Salinity difference');
xlabel('Month'); xlim([1 12])
set(gca,'fontsize',14)

end