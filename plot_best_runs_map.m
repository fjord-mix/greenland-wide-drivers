function hf_map = plot_best_runs_map(data_path,ensemble_yr,res_box_yr)

fsize=16;
% rmse_threshold = 0.5;
n_years   = length(ensemble_yr);
% n_layers  = ensemble_yr{1}(1,1).p.N;
lcolor    = lines(n_years);
hb_colors = [];
lbl_yrs   = {};
% fw_discharge_stacked = [];
% fw_export_stacked    = [];
% zfw_stacked         = [];
% znb_stacked         = [];
% tfw_stacked         = [];
% t_qsg_stacked       = [];

proj      = projcrs(3413,"Authority",'EPSG');
% z_general = 5:10:1205;

%% Getting the fjord outline

path_bedmachine = [data_path,'/greenland_common/BedMachineGreenland-v5.nc']; % DOI: https://doi.org/10.5067/GMEVBWFLWA7X
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
hf_map = figure('Position',[50 50, 500 700]); 
% tiledlayout(2,4,'TileSpacing','Compact');

% Plotting Greenland
% nexttile(1,[1 1]); hold on;
hold on
% title('(a) Simulated fjords','fontsize',14)
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
bubblesize([10 30])
bubblelim([0 2.])
set(gca,'fontsize',fsize)
blgd = bubblelegend('RMSE (^oC)','Location','northeast','fontsize',fsize);
hl = legend(hb_colors,lbl_yrs,'location','southwest','fontsize',fsize);
% hl.Layout.Tile = 5;
blgd.Title.FontWeight='normal';
set(blgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
set(hl.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
set(gcf, 'color', 'none');    
set(gca, 'color', 'none');
end