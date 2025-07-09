function hf_grid = plot_best_runs_grid(ensemble_yr,res_box_yr,fjord_matrix)

fsize=18;
% rmse_threshold = 0.5;
n_years   = length(ensemble_yr);
% n_layers  = ensemble_yr{1}(1,1).p.N;
lcolor    = lines(n_years);
hb_colors = [];
lbl_yrs   = {};

grid_rmse_t = NaN([max(fjord_matrix.ID)+1,n_years]);
grid_rmse_s = NaN(size(grid_rmse_t));
grid_rmse_2 = NaN(size(grid_rmse_t));

%% Create colormap (red-yellow-green)
A = imread('/Users/mabrag/scripts/import/matlab/gyrcb.png');
basect = im2double(permute(A(40,24:517,:),[2 3 1]));
N = 256;
nb = size(basect,1);
CT1 = interp1(1:nb,basect,linspace(1,nb,N)); % interpolate to get a specified-length version

% y_grid =1;
%% Summary figure
% hf_grid = figure('Position',[50 50, 500 700]); 
% hold on; box on;
% Adding the bubbles
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    res_box  = res_box_yr{i_year};
    n_fjords = size(ensemble,1);
    x_fjd=NaN([1,n_fjords]);
    y_fjd=NaN(size(x_fjd));
    rmse_fjd = NaN(size(x_fjd));
    size_marker = NaN(size(rmse_fjd));
    for i_fjord=1:n_fjords
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s) %&& res_box(i_fjord).rmse_tf(i_run,2) < rmse_threshold

                w_rmse_t  = 0.5;
                z_rmse_t  = normalize(res_box(i_fjord).rmse_tf(:,2),"range");
                z_rmse_s  = normalize(res_box(i_fjord).rmse_sf(:,2),"range");
                rmse_both = (z_rmse_t + z_rmse_s)/2;
                rmse_fjd(i_fjord) = min(rmse_both,[],'omitnan');

                % Creating a proper grid with T and S RMSEs
                x_grid = int32(str2double(res_box(i_fjord).id));% - 64; % converting the fjord ID ('A' to 'N') to an integer                
                grid_rmse_t(x_grid+1,i_year) = min(res_box(i_fjord).rmse_tf(:,2),[],'omitnan');
                grid_rmse_s(x_grid+1,i_year) = min(res_box(i_fjord).rmse_sf(:,2),[],'omitnan');
                grid_rmse_2(x_grid+1,i_year) = min(rmse_both,[],'omitnan');

                
                
                % This is for the bubble chart if we use it
                

                % rmse_fjd(i_fjord) = min(res_box(i_fjord).rmse_tf(:,2),[],'omitnan');
                y_fjd(i_fjord) = int32(str2double(res_box(i_fjord).id));% - 64; % converting the fjord ID ('A' to 'N') to an integer
                x_fjd(i_fjord) = 2015+i_year;
                size_marker(i_fjord) = 0.5-rmse_fjd(i_fjord);
                
                % if best_fjord_params(i_fjord).rmse_t > max_rmse, max_rmse = best_fjord_params(i_fjord).rmse_t; end
                % if best_fjord_params(i_fjord).rmse_t < min_rmse, min_rmse = best_fjord_params(i_fjord).rmse_t; end

                break
            end
        end
    end % i_fjord
    % hb = bubblechart(x_fjd,y_fjd,rmse_fjd,lcolor(i_year,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5);
    % hb_colors = [hb_colors hb];
    % lbl_yrs{end+1} = sprintf("%d (n=%d)",2015+i_year,n_fjords);
    lbl_yrs{end+1} = sprintf("%d",2015+i_year);
    % y_grid = y_grid+2;

end % i_year

cell_rmse_t_filt = {};
cell_rmse_s_filt = {};
cell_rmse_2_filt = {};
x_fjords_filt = {};
for i_fjord=1:size(grid_rmse_t,1)
    if sum(isnan(grid_rmse_t(i_fjord,:))) == size(grid_rmse_t,2)
        continue
    end
    cell_rmse_t_filt{end+1} = grid_rmse_t(i_fjord,:);
    cell_rmse_s_filt{end+1} = grid_rmse_s(i_fjord,:);
    cell_rmse_2_filt{end+1} = grid_rmse_2(i_fjord,:);
    x_fjords_filt{end+1} = i_fjord;
end

grid_rmse_t_filt = NaN([length(cell_rmse_t_filt),n_years]);
grid_rmse_s_filt = NaN(size(grid_rmse_t_filt));
grid_rmse_2_filt = NaN(size(grid_rmse_t_filt));
for i_fjord=1:length(cell_rmse_t_filt)
    grid_rmse_t_filt(i_fjord,:) = cell_rmse_t_filt{i_fjord};
    grid_rmse_s_filt(i_fjord,:) = cell_rmse_s_filt{i_fjord};
    grid_rmse_2_filt(i_fjord,:) = cell_rmse_2_filt{i_fjord};
end


vert_lines_ids = [0 32 68 94 105 111];
vert_lines_x = zeros(size(vert_lines_ids));
rmse_x = cell2mat(x_fjords_filt)-1;%+0.5; %0.5:30:150+0.5;
rmse_y = 1.5:1:n_years+0.5;
fjord_lbl = {};
for i=1:length(rmse_x)
    fjord_lbl{i} = num2str(rmse_x(i));
end
x_plot_axis = 1:length(fjord_lbl);
for i=1:length(vert_lines_ids)
    vert_lines_x(i) = find(rmse_x == vert_lines_ids(i));
end
hf_grid = figure('Position',[50 50, 1200 900]); 
ht = tiledlayout(3,1);%,'TileSpacing','tight');%,'Padding','compact');
nexttile(1)
% figure;
hold on;
h = imagesc(x_plot_axis,rmse_y,grid_rmse_t_filt');
set(h,'AlphaData', 1-isnan(grid_rmse_t_filt'))
clim([0 2.5])
xlim([x_plot_axis(1)-0.5 x_plot_axis(end)+0.5]); ylim([1 n_years+1])
box on;
set(gca,'XTick',1:length(fjord_lbl),'XtickLabel',{});
set(gca,'Ytick',rmse_y,'YTickLabel',flip(lbl_yrs))
hcb = colorbar;
ylabel(hcb,'RMSE (^oC)');
set(gca,'FontSize',fsize)
text(0.01,1.01,"(a) Temperature",'units','normalized','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',fsize)
% ylabel('Year')
text(vert_lines_x(1)+0.5,1.5,"NW Greenland",'fontsize',fsize-4)
plot([vert_lines_x(2)+0.5 vert_lines_x(2)+0.5],[0 n_years+1],'-k')
text(vert_lines_x(3),1.5,"SE Greenland",'fontsize',fsize-4)
plot([vert_lines_x(3)-0.5 vert_lines_x(3)-0.5],[0 n_years+1],'-k')
text(vert_lines_x(5),1.5,"Scoresby Sund",'fontsize',fsize-4)
plot([vert_lines_x(5)-0.5 vert_lines_x(5)-0.5],[0 n_years+1],'-k')
plot([vert_lines_x(6)+0.5 vert_lines_x(6)+0.5],[0 n_years+1],'-k')


nexttile(2)
hold on;
h = imagesc(x_plot_axis,rmse_y,grid_rmse_s_filt');
set(h,'AlphaData', 1-isnan(grid_rmse_s_filt'))
clim([0 2.5])
xlim([x_plot_axis(1)-0.5 x_plot_axis(end)+0.5]); ylim([1 n_years+1])
box on;
set(gca,'XTick',1:length(fjord_lbl),'XtickLabel',{});
set(gca,'Ytick',rmse_y,'YTickLabel',flip(lbl_yrs))
% ylabel('Year')
set(gca,'FontSize',fsize)
text(0.01,1.01,"(b) Salinity",'units','normalized','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',fsize)
plot([vert_lines_x(2)+0.5 vert_lines_x(2)+0.5],[0 n_years+1],'-k')
plot([vert_lines_x(3)-0.5 vert_lines_x(3)-0.5],[0 n_years+1],'-k')
plot([vert_lines_x(5)-0.5 vert_lines_x(5)-0.5],[0 n_years+1],'-k')
plot([vert_lines_x(6)+0.5 vert_lines_x(6)+0.5],[0 n_years+1],'-k')
hcb = colorbar;
ylabel(hcb,'RMSE (g/Kg)');

nexttile(3)
hold on;
h = imagesc(x_plot_axis,rmse_y,grid_rmse_2_filt');
set(h,'AlphaData', 1-isnan(grid_rmse_2_filt'))
clim([0 0.12])
xlim([x_plot_axis(1)-0.5 x_plot_axis(end)+0.5]); ylim([1 n_years+1])
box on;
set(gca,'XTick',1:length(fjord_lbl),'XtickLabel',rmse_x);
set(gca,'Ytick',rmse_y,'YTickLabel',flip(lbl_yrs))

colormap(flip(colormap('hot')))
% xlabel('Fjord')
% ylabel('Year')
set(gca,'FontSize',fsize)
text(0.01,1.01,"(c) Averaged normalised",'units','normalized','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',fsize)
plot([vert_lines_x(2)+0.5 vert_lines_x(2)+0.5],[0 n_years+1],'-k')
plot([vert_lines_x(3)-0.5 vert_lines_x(3)-0.5],[0 n_years+1],'-k')
plot([vert_lines_x(5)-0.5 vert_lines_x(5)-0.5],[0 n_years+1],'-k')
plot([vert_lines_x(6)+0.5 vert_lines_x(6)+0.5],[0 n_years+1],'-k')
colormap(CT1)
hcb = colorbar;
ylabel(hcb,'RMSE (-)');
ylabel(ht,'Year','fontsize',fsize)
xlabel(ht,'Fjord','fontsize',fsize)
end