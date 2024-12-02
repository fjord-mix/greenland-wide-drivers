function [hf_dst,hf_loc] = plot_ocn_cast_pairs(folder_ctd_casts,fjord_matrix,res_box_yr)

% mcolor=lines(5);
mcolor   = cmocean('thermal',5);
% lcolor = lines(2);
leg_handles = [];
leg_labels  = {};
flon=NaN([5,height(fjord_matrix)]);
flat=NaN(size(flon));
slon=NaN(size(flon));
slat=NaN(size(flon));
rmse_fjord=NaN([5,length(flon)]);

hf_dst = figure; box on; hold on;
for which_year=2016:2020
    first_plot=1;
    i_year = which_year-2015;
    res_box = res_box_yr{i_year};
    res_box_ids=NaN(size(res_box));
    for i_res_box=1:length(res_box)
        res_box_ids(i_res_box) = int32(str2double(res_box(i_res_box).id));
    end
    
    for i_fjord=1:height(fjord_matrix)
        eval(sprintf('id_cast_shelf = num2str(fjord_matrix.shelf_%d(i_fjord));',which_year));
        eval(sprintf('id_cast_fjord = num2str(fjord_matrix.fjord_%d(i_fjord));',which_year));
        if ~strcmp(id_cast_shelf,'NaN') && ~strcmp(id_cast_fjord,'NaN')
            % digitised_id = find([fjords_digitised.id] == fjord_matrix.ID(i_fjord));
            id_simulated_fjord = find(res_box_ids == fjord_matrix.ID(i_fjord));
            rmse_fjord(i_year,i_fjord) = min(res_box(id_simulated_fjord).rmse_tf(:,2));

            omg_data_shelf = dir([folder_ctd_casts,'/*',id_cast_shelf,'*.nc']);
            omg_data_fjord = dir([folder_ctd_casts,'/*',id_cast_fjord,'*.nc']);

            flon(i_year,i_fjord)=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lon');
            flat(i_year,i_fjord)=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lat');

            slon(i_year,i_fjord)=ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'lon');
            slat(i_year,i_fjord)=ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'lat');

            dist_between_casts = m_lldist([flon(i_year,i_fjord),slon(i_year,i_fjord)],[flat(i_year,i_fjord),slat(i_year,i_fjord)]);
            % dist_between_casts = dist_between_casts-(res_box(id_simulated_fjord).L.*1e-3); % subtract fjord distance
            if first_plot
                hs = scatter(fjord_matrix.ID(i_fjord),dist_between_casts,150,'o','MarkerFaceColor',mcolor(which_year-2015,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                first_plot=0;
                leg_handles = [leg_handles hs];
                leg_labels{end+1} = num2str(which_year);
            else
                scatter(fjord_matrix.ID(i_fjord),dist_between_casts,100,'o','MarkerFaceColor',mcolor(which_year-2015,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
            end
        end
    end
end
% ylabel(hc,'Year')
legend(leg_handles,leg_labels,'location','northwest')
xlabel('Fjord','fontsize',16)
ylabel('Distance between casts (km)','fontsize',16)
set(gca,'fontsize',16)
xlim([-2 125])


% flon(isnan(flon)) = [];
% flat(isnan(flat)) = [];
% slon(isnan(slon)) = [];
% slat(isnan(slat)) = [];

hf_loc = figure('Position',[50 50 1000 1000]); box on; hold on;

tiledlayout('flow','TileSpacing','Compact');
nexttile; hold on
m_proj('UTM','lon',[-75 -10],'lat',[59.5 84]); % GRL wide
m_gshhs_h('color','k','linewidth',0.5);
m_grid;
plot_pts(slon,slat,flon,flat,rmse_fjord);%,lcolor)
hc = colorbar('southoutside');
ylabel(hc,'RMSE (^oC)','fontsize',14)
clim([0 2])

nexttile; hold on
m_proj('UTM','lon',[-75 -46],'lat',[67 84]); % NO/NW/CW
m_gshhs_h('color','k','linewidth',0.5);
m_grid;
plot_pts(slon,slat,flon,flat,rmse_fjord);%,lcolor)
clim([0 2])

nexttile; hold on
m_proj('UTM','lon',[-46 -30],'lat',[59.5 67]); % SE
m_gshhs_h('color','k','linewidth',0.5);
m_grid;
plot_pts(slon,slat,flon,flat,rmse_fjord);%,lcolor)
clim([0 2])

nexttile; hold on
m_proj('UTM','lon',[-46 -12],'lat',[67 84]); % CE/NO
m_gshhs_h('color','k','linewidth',0.5);
m_grid;
plot_pts(slon,slat,flon,flat,rmse_fjord);%,lcolor)
clim([0 2])
end

function plot_pts(slon,slat,flon,flat,rmse_fjord,~)
    m_scatter(slon,slat,40,'^','filled','MarkerFaceColor',[1 0 0])
    % for i=1:size(flon,1)
    for i=size(flon,1):-1:1
        for j=1:size(flon,2)
            m_line([flon(i,j),slon(i,j)],[flat(i,j),slat(i,j)],'linewidth',1.5,'linestyle','--','color','k')
            m_scatter(flon(i,j),flat(i,j),80,rmse_fjord(i,j),'o','filled','MarkerFaceAlpha',0.75)
        end
    end
end