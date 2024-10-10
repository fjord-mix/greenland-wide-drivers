function hf_zfw = plot_fw_simulated_fjords(data_path,ensemble_yr,res_box_yr)

rmse_threshold = 0.5;
n_years = length(ensemble_yr);
lcolor=lines(n_years);
h_hist = [];
lbl_hist = {};
zfw_stacked = [];
znb_stacked = [];
tfw_stacked = [];
t_qsg_stacked = [];
proj = projcrs(3413,"Authority",'EPSG');

%% Getting the fjord outline

path_bedmachine = [data_path,'/greenland/BedMachineGreenland-v5.nc']; % DOI: https://doi.org/10.5067/GMEVBWFLWA7X
x       = single(ncread(path_bedmachine,'x'));
y       = single(ncread(path_bedmachine,'y'));
mask    = single(ncread(path_bedmachine,'mask'));
surface = ncread(path_bedmachine,'surface');

x = x(1:10:end);
y = y(1:10:end);
topo    = surface(1:10:end,1:10:end);
mask    = mask(1:10:end,1:10:end);
mask(mask == 4) = 0;
mask(mask==2 | mask==3) = 1;
topo(mask==0) = NaN;



%% Summary figure
hf_zfw = figure('Position',[50 50, 1200 800]); 
ht_main = tiledlayout(2,4,'TileSpacing','Compact');

nexttile(1,[1 1]); hold on;
title('(a) Simulated fjords','fontsize',14)
hi_sfc = imagesc2(x,y,topo'); box on; axis xy equal %tight; 
% hi_sfc = contour(xplt,yplt,outline,'-k','linewidth',1.5); box on; axis xy equal %tight; 
clim([0 2.5e3])
colormap('gray');
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
    bubblechart(x_fjd,y_fjd,rmse_fjd,lcolor(i_year,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5);
end % i_year
bubblesize([5 15])
bubblelim([0 2.1])
% blgd = bubblelegend('RMSE (^oC)','Style','telescopic','Location','southeast','fontsize',12);
blgd = bubblelegend('RMSE (^oC)','Location','southwest','fontsize',8);
blgd.Title.FontWeight='normal';
blgd.Title.FontSize=8;
set(blgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

%% Histograms
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    
    zfw = NaN(size(ensemble));
    znb = NaN(size(ensemble));
    t_fw_max = NaN(size(ensemble));
    t_qsg_max = NaN([size(ensemble,1),1]);
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s) %&& res_box(i_fjord).rmse_tf(i_run,2) < rmse_threshold
                % depths at target day
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_tgt_day);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb(i_tgt_day);

                % % average depth over the last year
                % zfw(i_fjord,i_run) =
                % mean(ensemble(i_fjord,i_run).s.z_max_export_t(end-365:end)); % average of peak FW exp depth
                zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export;   % depth of the average peak FW exp
                znb(i_fjord,i_run) = mean(ensemble(i_fjord,i_run).s.znb_t(end-365:end));
                % 
                % % depth of export maximum over the last year
                [~,i_fz] = min(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export_t(end-365+i_fz);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb_t(end-365+i_fz);

                t_fw_max(i_fjord,i_run) = i_fz;
                [~,t_qsg_max(i_fjord)] = max(ensemble(i_fjord,i_run).s.Qsg(end-365:end));

                zfw_stacked = [zfw_stacked; zfw(:)];
                znb_stacked = [znb_stacked; znb(:)];
                tfw_stacked = [tfw_stacked; t_fw_max(:)];
            end
        end
    end
    t_qsg_stacked = [t_qsg_stacked; t_qsg_max(:)];

    nexttile(2,[1 1]); hold on; box on; grid on;
    histogram(t_fw_max(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','vertical','FaceAlpha',0.5);
    set(gca,'fontsize',16)
    
    nexttile(6,[1 1]); hold on; box on; grid on;
    histogram(t_qsg_max(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','vertical','FaceAlpha',0.5);
    set(gca,'fontsize',16,'ydir','reverse')

    nexttile(3,[2 1]); hold on; box on; grid on;
    histogram(znb(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca, 'xdir', 'reverse','fontsize',16);
    ylim([-450 0])

    nexttile(4,[2 1]); hold on; box on; grid on;
    h1 = histogram(zfw(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca,'fontsize',16)
    ylim([-450 0])

    h_sit = [h_hist h1];
    lbl_hist{end+1} = num2str(2015+i_year);
end
nexttile(2,[1 1]); hold on;
title('(b) Maximum FW export','fontsize',14)
histogram(tfw_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','vertical','FaceAlpha',1.);
xlim([100 300])
ylabel('Probability','fontsize',16)

nexttile(6,[1 1]); hold on;
title('(c) Maximum subglacial discharge','fontsize',14)
histogram(t_qsg_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','vertical','FaceAlpha',1.);
xlim([100 300])
ylabel('Probability','fontsize',16)
xlabel('Day of year','fontsize',16)

nexttile(3,[2 1]); hold on;
title('(d) Neutral buoyancy','fontsize',14)
histogram(znb_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);
xlabel('Probability','fontsize',16)
ylabel('Depth (m)','fontsize',16)

nexttile(4,[2 1]); hold on;
title('(e) Maximum FW export','fontsize',14)
histogram(zfw_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);
xlabel('Probability','fontsize',16)
legend(h_hist,lbl_hist,'location','southeast')

end