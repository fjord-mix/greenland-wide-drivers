function plot_reg_ocn_forcings(datasets,fjords_compilation,temp_amp_seasonal)
letters = {'a','b','c','d','e','f','g','h'};
regions_lbl = {'SW','SE','CW','CE','NW','NE','NO'};

if nargin < 3 || isempty(temp_amp_seasonal)
    temp_amp_seasonal=0.5; % temperature amplitude threshold to be considered as the "seasonal layer"
end

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

[temp_forcing, ~, ~, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
[salt_forcing, ~, ~, ~] = get_var_forcing_by_region(fjords_processed,'Ss');

% fjord_rho = (fjords_processed(1).p.betaS*salt_forcing - fjords_processed(1).p.betaT*temp_forcing);
% sc_reg = squeeze(trapz(depths,salt_forcing.* fjord_rho,2)./max(abs(depths)));
% squeeze(trapz(depths,(temp_forcing+273.15).* fjord_rho,2)./max(abs(depths)).* fjords_processed(1).p.cw);

% time_axis = datasets.opts.time_start:datasets.opts.dt:datasets.opts.time_end;
region_line_color = lines(7);
time_axis     = fjords_processed(1).t;

depth_seasonal=NaN([1,size(temp_forcing,3)]);
for i_reg=1:length(regions_lbl)    
    temp_amp_seasonal = std(detrend(temp_forcing(:,end,i_reg)));
    for i_z=length(depths):-1:1 % will go down the profile until it finds the deepest amplitude
        temp_detrend=detrend(temp_forcing(:,i_z,i_reg));
        % if (max(temp_detrend)-min(temp_detrend)) > temp_amp_seasonal
        if std(temp_detrend) > 0.15*temp_amp_seasonal
            depth_seasonal(i_reg) = depths(i_z);
        end
    end
end

%% Plotting Temperature and salinity
figure('Name','Temperature undisturbed','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg)
    imagesc(time_axis,-depths,temp_forcing(:,:,i_reg)');
    % yline(-depth_seasonal(i_reg),'--r','linewidth',2)
    clim([-2.5 7]);
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('time (days)'); end
    set(gca,'fontsize',14)
    % p = polyfit(time_axis,hc_reg(:,i_reg),1);
    % trend_ohc = p(1)*1e3;
    text(0.05,1.07,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
colormap(cmocean('thermal'))
hc = colorbar('fontsize',14);
ylabel(hc,'Temperature (^oC)');
hc.Position(1)=hc.Position(1)+0.15;


figure('Name','Salinity undisturbed','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg)
    imagesc(time_axis,-depths,salt_forcing(:,:,i_reg)');
    clim([30 36]);
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('time (days)'); end
    set(gca,'fontsize',14)
    % p = polyfit(time_axis,sc_reg(:,i_reg),1);
    % trend_osc = p(1)*1e3;

    % text(0.05,1.07,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
    text(0.05,1.07,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
colormap(cmocean('haline'))
hc = colorbar('fontsize',14);
ylabel(hc,'Salinity');
hc.Position(1)=hc.Position(1)+0.15;


%% Plotting subglacial and solid-ice discharge
qsg_reg=NaN([length(fjords_processed),length(time_axis),length(regions_lbl)]);
d_reg=NaN(size(qsg_reg));
hc_reg=NaN(size(qsg_reg));
sc_reg=NaN(size(qsg_reg));
for i_fjord=1:length(fjords_processed)
    fjord = fjords_processed(i_fjord);
    qsg_reg(i_fjord,:,fjord.m.regionID) = fjord.f.Qsg;
    d_reg(i_fjord,:,fjord.m.regionID)   = fjord.f.D;

    % compute heat and salt contents per unit area
    % fjord_rho = (fjord.p.betaS*fjord.f.Ss - fjord.p.betaT*fjord.f.Ts);
    % sc_reg(i_fjord,:,fjord.m.regionID) = trapz(fjord.f.zs,fjord.f.Ss.* fjord_rho)./fjord.p.H;
    % hc_reg(i_fjord,:,fjord.m.regionID) = trapz(fjord.f.zs,(fjord.f.Ts+273.15).* fjord.p.cw .* fjord_rho)./fjord.p.H;
end


time_axis_dq = fjords_processed(1).m.time_axis';
region_handles = [];
figure('Name','Subglacial and solid-ice discharge','Position',[20 20 1200 400])
for i_reg=1:length(regions_lbl)
    
    subplot(1,2,1); hold on; box on
    ylabel('Subglacial discharge (m^{3}s^{-1})'); 
    % lower_bnd = prctile(qsg_reg(:,:,i_reg),25,1);
    % upper_bnd = prctile(qsg_reg(:,:,i_reg),75,1);
    % median_ln = median(qsg_reg(:,:,i_reg),1,'omitnan');
    mean_ln = mean(bootstrp(100,@(x)[mean(x,1,'omitnan')],qsg_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    std_ln  = std(bootstrp(100,@(x)[mean(x,1,'omitnan')],qsg_reg(:,:,i_reg)),'omitnan'); % multiplying by 1e-3 to change units to kJ
    upper_bnd = mean_ln+std_ln;
    lower_bnd = mean_ln-std_ln;
    x2 = [time_axis_dq, fliplr(time_axis_dq)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];

    hp = fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(time_axis_dq,mean_ln,'Color',region_line_color(i_reg,:)); 
    region_handles=[region_handles hp];
    ylim([0 2e3])

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