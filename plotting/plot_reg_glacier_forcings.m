function plot_reg_glacier_forcings(datasets,fjords_compilation,glaciers_compilation)

letters = {'a','b','c','d','e','f','g','h'};
regions_lbl = {'SW','SE','CW','CE','NW','NE','NO'};
region_line_color = lines(length(regions_lbl));

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end
time_axis     = fjords_processed(1).t;
%% Plotting subglacial and solid-ice discharge

if datasets.opts.restrict_to_fjords
    qsg_reg=NaN([length(fjords_processed),length(time_axis),length(regions_lbl)]);
    d_reg=NaN(size(qsg_reg));
    for i_fjord=1:length(fjords_processed)
        fjord = fjords_processed(i_fjord);
        qsg_reg(i_fjord,:,fjord.m.regionID) = fjord.f.Qsg;
        d_reg(i_fjord,:,fjord.m.regionID)   = fjord.f.D;
    end
else
    runtime_axis = [datetime(datasets.opts.time_start):datasets.opts.dt:datetime(datasets.opts.time_end)]';
    qsg_reg=NaN([length(glaciers_compilation),length(time_axis),length(regions_lbl)]);
    d_reg=NaN(size(qsg_reg));
    for i_glacier=1:length(glaciers_compilation)
        glacier = glaciers_compilation(i_glacier);

        time_discharge = glacier.iceberg.time_axis;
        seconds_in_year = eomday(time_discharge.Year,time_discharge.Month)*86400*12; % we need to convert from total discharge in the month to discharge per second
        Gt_to_m3_ice = 1e9*1e3/916.7; % Giga * kg / (kg/m^3) assuming glacier ice density of 916.7 kg/m^3
        D_s = glacier.iceberg.discharge./seconds_in_year .* Gt_to_m3_ice;
        [~] = data_overlap_check(runtime_axis,time_discharge,'iceberg'); % simple check that the data is available for the simulation period
        
        d_reg(i_glacier,:,glacier.regionID)   = interp1(time_discharge,D_s,runtime_axis,'linear','extrap'); 
        qsg_reg(i_glacier,:,glacier.regionID) = get_total_glacier_runoff(glacier,runtime_axis);
    end

end

time_axis_dq = runtime_axis';%fjords_processed(1).m.time_axis';
region_handles = [];
figure('Name','Subglacial and solid-ice discharge','Position',[20 20 1200 400])
for i_reg=1:length(regions_lbl)
    
    subplot(1,2,1); hold on; box on
    ylabel('Subglacial discharge (m^{3}s^{-1})'); 
    % lower_bnd = prctile(qsg_reg(:,:,i_reg),25,1);
    % upper_bnd = prctile(qsg_reg(:,:,i_reg),75,1);
    % median_ln = median(qsg_reg(:,:,i_reg),1,'omitnan');
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],qsg_reg(:,:,i_reg)),'omitnan');
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],qsg_reg(:,:,i_reg)),'omitnan');
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;
    x2 = [time_axis_dq, fliplr(time_axis_dq)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(time_axis_dq,mean_ln,'Color',region_line_color(i_reg,:)); 
    region_handles=[region_handles hp];
    % ylim([0 1e3])

    subplot(1,2,2); hold on; box on
    ylabel('Solid-ice discharge (m^{3}s^{-1})');
    % lower_bnd = prctile(d_reg(:,:,i_reg),25,1);
    % upper_bnd = prctile(d_reg(:,:,i_reg),75,1);
    % median_ln = median(d_reg(:,:,i_reg),1,'omitnan');
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],d_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],d_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;
    x2 = [time_axis_dq, fliplr(time_axis_dq)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.1);
    plot(time_axis_dq,mean_ln,'Color',region_line_color(i_reg,:)); 

    % subplot(2,2,3); hold on; box on
    % ylabel('Shelf heat content (kJ m^{-3})',fontsize=14);
    % % lower_bnd = 1e-3.*prctile(hc_reg(:,:,i_reg),25,1);
    % % upper_bnd = 1e-3.*prctile(hc_reg(:,:,i_reg),75,1);
    % % median_ln = 1e-3.*median(hc_reg(:,:,i_reg),1,'omitnan');
    % mean_ln = 1e-3.*mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],hc_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    % std_ln  = 1e-3.*std(bootstrp(100,@(x)[mean(x,1,'omitnan')],hc_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    % upper_bnd = mean_ln+std_ln;
    % lower_bnd = mean_ln-std_ln;
    % x2 = [time_axis_dq, fliplr(time_axis_dq)];
    % inBetween = [lower_bnd, fliplr(upper_bnd)];
    % 
    % % hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.1);
    % plot(time_axis_dq,mean_ln,'Color',region_line_color(i_reg,:)); 
    % 
    % subplot(2,2,4); hold on; box on
    % ylabel('Shelf salt content (kg m^{-3})',fontsize=14);
    % % lower_bnd = 1e-3.*prctile(sc_reg(:,:,i_reg),25,1);
    % % upper_bnd = 1e-3.*prctile(sc_reg(:,:,i_reg),75,1);
    % % median_ln = 1e-3.*median(sc_reg(:,:,i_reg),1,'omitnan');
    % mean_ln = 1e-3.*mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],hc_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    % std_ln  = 1e-3.*std(bootstrp(100,@(x)[mean(x,1,'omitnan')],hc_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    % upper_bnd = mean_ln+std_ln;
    % lower_bnd = mean_ln-std_ln;
    % x2 = [time_axis_dq, fliplr(time_axis_dq)];
    % inBetween = [lower_bnd, fliplr(upper_bnd)];
    % 
    % % hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.1);
    % plot(time_axis_dq,mean_ln,'Color',region_line_color(i_reg,:)); 
end
hl = legend(region_handles,regions_lbl,'fontsize',14,'Location','northeast');
hl.NumColumns=3;
subplot(1,2,1); 
xlabel('Time'); text(0.02,1.07,'(a)','units','normalized','fontsize',14); set(gca,'fontsize',14)
subplot(1,2,2); 
xlabel('Time'); text(0.02,1.07,'(b)','units','normalized','fontsize',14); set(gca,'fontsize',14)
% subplot(2,2,3); 
% xlabel('Time');
% text(0.02,1.07,'(c)','units','normalized','fontsize',14); set(gca,'fontsize',14)
% subplot(2,2,4); 
% xlabel('Time');
% text(0.02,1.07,'(d)','units','normalized','fontsize',14); set(gca,'fontsize',14)

end