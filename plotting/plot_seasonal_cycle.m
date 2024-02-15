function hf = plot_seasonal_cycle(datasets,fjords_compilation,tt_ensemble,regions_lbl,verbose)

if nargin < 5, verbose=0; end

%% Plotting auxiliary variables
n_regions=length(regions_lbl);
region_handles = [];
sf_handles = [];
region_line_color = lines(n_regions);

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings

%% Getting subglacial discharge
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end
% time_axis = fjords_processed(1).t;
time_axis = datasets.opts.time_start:datasets.opts.dt:datasets.opts.time_end;
% time_axis_plt = time_axis(2:end)';

qsg_reg = NaN([length(fjords_processed),length(time_axis),length(regions_lbl)]);
d_reg   = NaN(size(qsg_reg));
for i_fjord=1:length(fjords_processed)
    fjord = fjords_processed(i_fjord);
    qsg_reg(i_fjord,:,fjord.m.regionID) = fjord.f.Qsg;
    d_reg(i_fjord,:,fjord.m.regionID)   = fjord.f.D;
end


tt_q = cell(size(tt_ensemble));
for i_reg=1:n_regions
    % mean_ln_qsg = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],qsg_reg(:,:,i_reg)),'omitnan');
    mean_ln_qsg = mean(qsg_reg(:,:,i_reg),1,'omitnan');
    mean_ln_d   = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],d_reg(:,:,i_reg)),'omitnan');
    tt_q{i_reg} = timetable(time_axis',mean_ln_qsg',mean_ln_d','VariableNames',{'Qsg','D'});
end


%% Plotting
hf = figure('Name','Seasonal cycles','Position',[20 20 1000 600]);
for i_reg=1:n_regions
    clim_mon = groupsummary(tt_ensemble{i_reg},"Time","monthofyear","mean");

    subplot(2,2,1); hold on; grid on; box on;
    hp_r = plot(clim_mon.mean_Tf,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    plot(clim_mon.mean_Ts,'LineStyle','--','Color',region_line_color(i_reg,:),'linewidth',1.5);

    subplot(2,2,2); hold on; grid on; box on;
    if i_reg==1
        hp_f = plot(clim_mon.mean_Sf,'-k','linewidth',1.);
        hp_s = plot(clim_mon.mean_Sf,'--k','linewidth',1.);
    end
    plot(clim_mon.mean_Sf,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    plot(clim_mon.mean_Ss,'LineStyle','--','Color',region_line_color(i_reg,:),'linewidth',1.5);

    subplot(2,2,3); hold on; grid on; box on;
    plot(clim_mon.mean_dT,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);
    
    subplot(2,2,4); hold on; grid on; box on;
    plot(clim_mon.mean_dS,'LineStyle','-','Color',region_line_color(i_reg,:),'linewidth',1.5);

    region_handles = [region_handles hp_r];
    sf_handles = [hp_f hp_s];

    if verbose
        [max_dt,imax_dt] = max(clim_mon.mean_dT);
        [min_dt,imin_dt] = min(clim_mon.mean_dT);
        [max_tf,imax_tf] = max(clim_mon.mean_Tf);
        [min_tf,imin_tf] = min(clim_mon.mean_Tf);
        [max_ts,imax_ts] = max(clim_mon.mean_Ts);
        [min_ts,imin_ts] = min(clim_mon.mean_Ts);
        [max_sf,imax_sf] = max(clim_mon.mean_Sf);
        [min_sf,imin_sf] = min(clim_mon.mean_Sf);
        [max_ss,imax_ss] = max(clim_mon.mean_Ss);
        [min_ss,imin_ss] = min(clim_mon.mean_Ss);
        disp('======================================')
        fprintf('region: %s\n',regions_lbl{i_reg})
        fprintf('Tf max in %d and min in %d\n',imax_tf,imin_tf)
        fprintf('Ts max in %d and min in %d\n',imax_ts,imin_ts)
        fprintf('Amplitude of dT:%.2f\n',max_dt-min_dt)
        fprintf('Amplitude of Tf:%.2f\n',max_tf-min_tf)
        fprintf('Amplitude of Ts:%.2f\n',max_ts-min_ts)
        disp('')
        fprintf('Sf max in %d and min in %d\n',imax_sf,imin_sf)
        fprintf('Ss max in %d and min in %d\n',imax_ss,imin_ss)
    end
end
for i_reg=1:n_regions % making its own loop to get around Matlab adding extra markers to the lineplots
    qcli_mon = groupsummary(tt_q{i_reg},"Time","monthofyear","mean");
    subplot(2,2,3); hold on; grid on; box on
    yyaxis right
    x2 = [double(qcli_mon.monthofyear_Time), fliplr(double(qcli_mon.monthofyear_Time))];
    inBetween = [zeros(size(qcli_mon.mean_Qsg)), fliplr(qcli_mon.mean_Qsg)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    % plot(qcli_mon.mean_Qsg,'LineStyle','--','Color',region_line_color(i_reg,:),'linewidth',0.5);
    ax=gca;
    ax.YColor=[0 0 0];
    yyaxis left
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
yyaxis right
ylabel('Subglacial discharge (m^3s^{-1})','fontsize',14)
subplot(2,2,4)
text(0.02,1.05,'(d)','fontsize',14,'Units','normalized')
ylabel('Salinity difference');
xlabel('Month'); xlim([1 12])
set(gca,'fontsize',14)

end