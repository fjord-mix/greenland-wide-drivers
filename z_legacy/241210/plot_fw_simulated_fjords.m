function hf_zfw = plot_fw_simulated_fjords(data_path,ensemble_yr,res_box_yr)

rmse_threshold = 0.5;
n_years   = length(ensemble_yr);
n_layers  = ensemble_yr{1}(1,1).p.N;
lcolor    = lines(n_years);
hb_colors = [];
lbl_yrs   = {};
fw_discharge_stacked = [];
fw_export_stacked    = [];
% zfw_stacked         = [];
% znb_stacked         = [];
% tfw_stacked         = [];
% t_qsg_stacked       = [];

proj      = projcrs(3413,"Authority",'EPSG');
z_general = 5:10:1205;

%% Getting the fjord outline

path_bedmachine = [data_path,'/greenland/BedMachineGreenland-v5.nc']; % DOI: https://doi.org/10.5067/GMEVBWFLWA7X
x       = single(ncread(path_bedmachine,'x'));
y       = single(ncread(path_bedmachine,'y'));
mask    = single(ncread(path_bedmachine,'mask'));
surface = ncread(path_bedmachine,'surface');
x       = x(1:10:end);
y       = y(1:10:end);
topo    = surface(1:10:end,1:10:end);
mask    = mask(1:10:end,1:10:end);

mask(mask == 4) = 0;
mask(mask==2 | mask==3) = 1;
topo(mask==0) = NaN;

%% Summary figure
hf_zfw = figure('Position',[50 50, 1200 800]); 
tiledlayout(2,4,'TileSpacing','Compact');

% Plotting Greenland
nexttile(1,[1 1]); hold on;
title('(a) Simulated fjords','fontsize',14)
imagesc2(x,y,topo'); box on; axis xy equal %tight; 
clim([0 2.5e3])
colormap('gray');

% Adding the bubbles
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    res_box  = res_box_yr{i_year};
    n_fjords = size(ensemble,1);
    x_fjd=NaN([1,n_fjords]);
    y_fjd=NaN(size(x_fjd));
    rmse_fjd = NaN(size(x_fjd));
    for i_fjord=1:n_fjords
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s) %&& res_box(i_fjord).rmse_tf(i_run,2) < rmse_threshold
                [x_fjd(i_fjord),y_fjd(i_fjord)] = projfwd(proj,ensemble(i_fjord,i_run).m.lat,ensemble(i_fjord,i_run).m.lon);
                rmse_fjd(i_fjord) = min(res_box(i_fjord).rmse_tf(:,2),[],'omitnan');
                break
            end
        end
    end % i_fjord
    hb = bubblechart(x_fjd,y_fjd,rmse_fjd,lcolor(i_year,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5);
    hb_colors = [hb_colors hb];
    lbl_yrs{end+1} = sprintf("%d (n=%d)",2015+i_year,n_fjords);
end % i_year
bubblesize([5 15])
bubblelim([0 2.1])
blgd = bubblelegend('RMSE (^oC)','Location','southwest','fontsize',8);
hl = legend(hb_colors,lbl_yrs,'location','layout','fontsize',16);
hl.Layout.Tile = 5;
blgd.Title.FontWeight='normal';
set(blgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

%% Profiles
% for i_year=n_years:-1:1
for i_year=1:n_years
    ensemble = ensemble_yr{i_year};
    res_box  = res_box_yr{i_year};

    % zfw = NaN(size(ensemble));
    % znb = NaN(size(ensemble));
    % t_fw_max = NaN(size(ensemble));
    % t_qsg_max = NaN([size(ensemble,1),1]);
    fw_discharge_norm = NaN([length(z_general),size(ensemble,1)]);
    fw_export_norm = NaN(size(fw_discharge_norm));
    for i_fjord=1:size(ensemble,1)
        sum_fw_export = NaN([n_layers,size(ensemble,2)]);
        sum_fw_discharge = NaN(size(sum_fw_export));
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s) && res_box(i_fjord).rmse_tf(i_run,1) < rmse_threshold
                z_fjord = ensemble(i_fjord,i_run).s.z;
                sum_fw_export(:,i_run)    = ensemble(i_fjord,i_run).s.fw_profile_export/trapz(ensemble(i_fjord,i_run).s.fw_profile_export,ensemble(i_fjord,i_run).s.z);
                sum_fw_discharge(:,i_run) = ensemble(i_fjord,i_run).s.fw_profile_discharge/trapz(ensemble(i_fjord,i_run).s.fw_profile_discharge,ensemble(i_fjord,i_run).s.z);
                % sum_fw_export_t(:,i_run) = ensemble(i_fjord,i_run).s.fw_export_t(end-364:end);
                % sum_qsg(:,i_run)         = ensemble(i_fjord,i_run).s.Qsg(end-364:end);
                

                nexttile(2,[1 1]); hold on; box on; grid on;
                plot(1:365,-ensemble(i_fjord,i_run).s.fw_export_t(end-364:end),'Color',[lcolor(i_year,:),0.1],'linewidth',0.5);

                nexttile(6,[1 1]); hold on; box on; grid on;
                plot(1:365,ensemble(i_fjord,i_run).s.Qsg(end-364:end),'color',[lcolor(i_year,:),0.1],'linewidth',0.5)

                nexttile(3,[2 1]); hold on; box on; grid on;
                plot(-sum_fw_discharge(:,i_run),-ensemble(i_fjord,i_run).s.z,'color',[lcolor(i_year,:) 0.01],'linewidth',2.5)
                set(gca,'fontsize',16)

                nexttile(4,[2 1]); hold on; box on; grid on;
                plot(-sum_fw_export(:,i_run),-ensemble(i_fjord,i_run).s.z,'color',[lcolor(i_year,:) 0.01],'linewidth',2.5)
                set(gca,'fontsize',16)
                

            end
        end

        % fw_export_t_allruns = -median(sum_fw_export_t,2,'omitnan');
        % qsg_allruns         = median(sum_qsg,2,'omitnan');

        fw_export_allruns            = -mean(sum_fw_export,2,'omitnan');
        fw_discharge_allruns         = -mean(sum_fw_discharge,2,'omitnan');
        fw_discharge_norm(:,i_fjord) = interp1(z_fjord,fw_discharge_allruns,z_general,'linear');%,'extrap');
        fw_export_norm(:,i_fjord)    = interp1(z_fjord,fw_export_allruns,z_general,'linear');%,'extrap');

        % nexttile(2,[1 1]); hold on; box on; grid on;
        % plot(1:365,fw_export_t_allruns,'Color',[lcolor(i_year,:),0.5],'linewidth',0.5);
        % 
        % nexttile(6,[1 1]); hold on; box on; grid on;
        % plot(1:365,qsg_allruns,'color',[lcolor(i_year,:),0.5],'linewidth',0.5)

        % nexttile(3,[2 1]); hold on; box on; grid on;
        % plot(fw_discharge_allruns,-z_fjord,'color',[lcolor(i_year,:) 0.5],'linewidth',1.5)
        % 
        % nexttile(4,[2 1]); hold on; box on; grid on;
        % plot(fw_export_allruns,-z_fjord,'color',[lcolor(i_year,:) 0.5],'linewidth',1.5)

    end
    % t_qsg_stacked = [t_qsg_stacked; t_qsg_max(:)];
    fw_discharge_stacked = [fw_discharge_stacked, mean(fw_discharge_norm,2,'omitnan')];
    fw_export_stacked = [fw_export_stacked, mean(fw_export_norm,2,'omitnan')];

    % nexttile(3,[2 1]); hold on; box on; grid on;
    % plot(mean(fw_discharge_norm,2,'omitnan'),-z_general,'color',[lcolor(i_year,:) 1],'linewidth',2.5)
    % 
    % nexttile(4,[2 1]); hold on; box on; grid on;
    % plot(mean(fw_export_norm,2,'omitnan'),-z_general,'color',[lcolor(i_year,:) 1],'linewidth',2.5)

    % h_sit = [h_hist h1];
    % lbl_hist{end+1} = num2str(2015+i_year);
end
nexttile(2,[1 1]); hold on;
title('(b) Peak FW export','fontsize',14)
xlim([100 300]);
% ylim([0 7e3]);
ylabel('Freshwater export (m^3 s^{-1})','fontsize',16)
set(gca,'fontsize',16)

nexttile(6,[1 1]); hold on;
title('(c) Peak subglacial discharge','fontsize',14)
xlim([100 300])
% ylim([0 7e3]);
ylabel('Discharge (m^3 s^{-1})','fontsize',16)
xlabel('Day of year','fontsize',16)
set(gca,'fontsize',16)

nexttile(3,[2 1]); hold on;
title('(d) Plume to fjord','fontsize',14)
plot(mean(fw_discharge_stacked,2,'omitnan'),-z_general,'color',[0 0 0],'linewidth',2.5)
xlabel('Norm. mean freshwater flux','fontsize',16)
ylabel('Depth (m)','fontsize',16)
ylim([-450 0])
% xlim([-0.1 0.05])
set(gca,'fontsize',16)

nexttile(4,[2 1]); hold on;
title('(e) Fjord to ocean','fontsize',14)
plot(mean(fw_export_stacked,2,'omitnan'),-z_general,'color',[0 0 0],'linewidth',2.5)
xlabel('Norm. mean freshwater export','fontsize',16)
ylim([-450 0])
xlim([-0.5 0.5])
set(gca,'fontsize',16)
% legend(h_hist,lbl_hist,'location','southeast')

end