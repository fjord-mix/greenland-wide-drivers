function fig_handles = plot_fjords_sectors(datasets,fjords_map,fjords,glaciers)
% set(0,'DefaultFigureWindowStyle','docked') % will dock the figure by default
set(0,'DefaultFigureWindowStyle','normal') % will have the figure in a separate window
ss = 10; % for subsampling and faster plot browsing
glacier_line_color = lines(7); % get glacier color based on their region
handles_regions    = [];

x       = datasets.grl.x;
y       = datasets.grl.y;
mask    = datasets.grl.mask;
bed     = datasets.grl.bed;
surface = datasets.grl.surface;

ocean_depth = datasets.full_ocn.depth;
depth_range = ocean_depth > 25 & ocean_depth < 200;
ocean_surf = squeeze(mean(datasets.full_ocn.thetao(:,:,depth_range,:),3,'omitnan'));
ocean_trend = NaN([size(ocean_surf,1),size(ocean_surf,2)]);
for i=1:size(ocean_surf,1)
    for j=1:size(ocean_surf,2)
        p = polyfit(datasets.full_ocn.time,ocean_surf(i,j,:),1);
        ocean_trend(i,j) = p(1)*12; % from month^-1 to year^-1
    end
end



topo = surface;
topo(mask~=2 & mask ~=3) = NaN;

x = x(1:ss:end);
y = y(1:ss:end);
bed = bed(1:ss:end,1:ss:end);
topo = topo(1:ss:end,1:ss:end);


%% Creating figure
fig_handles.hf = figure('Position',[10 10 730 580]);
ax0 = gca;
ax0.Visible = 'off';
ax0.XTick = []; ax0.YTick = [];

ax1 = axes;  hold on;
hi_bed = imagesc(x,y,bed'); box on; axis xy equal %tight;
demcmap([min(bed) max(bed)])
ax2 = axes; hold on;
hi_sfc = imagesc2(x,y,topo'); box on; axis xy equal %tight; 
clim(ax2,[0 3e3])
colormap(ax2,'gray');
ax3 = axes; hold on;

%% Plotting each fjord
for i=1:length(fjords)
fjord = fjords(i);
[~] = scatter(ax3,fjord.x,fjord.y,50,glacier_line_color(fjord.glaciers(1).regionID,:),'filled','MarkerEdgeColor','black','MarkerFaceAlpha',0.9,'Marker','square');
% text(fjord.x,fjord.y,num2str(i),'fontsize',12,'color',[1 1 1]);

% show ocean point linked to fjord
% hp_fo = plot([fjord.ocean.x,fjord.x],[fjord.ocean.y,fjord.y],':','linewidth',1.,'Color','black');
% hs_o = scatter(ax3,fjord.ocean.x,fjord.ocean.y,25,mean(fjord.ocean.Tz(1,:)),'filled','MarkerEdgeColor','black','MarkerFaceAlpha',0.9,'Marker','diamond');
% hp_fo = plot([fjord.ocean.x,fjord.x],[fjord.ocean.y,fjord.y],':','linewidth',2.5,'Color',fjord_color(i,:));
% hs_o = scatter(gca,fjord.ocean.x,fjord.ocean.y,25,fjord_marker_color,'filled','MarkerEdgeColor','black','MarkerFaceColor',fjord_color(i,:),'MarkerFaceAlpha',0.9,'Marker','diamond');

if datasets.opts.restrict_to_fjords
% Plotting each linked glacier to fjord
if isfield(fjord,'glaciers')
glaciers = fjord.glaciers;
for j=1:length(glaciers)
    glacier_color = glacier_line_color(glaciers(j).regionID,:);
    if isfield(glaciers(j),'runoff') && datasets.opts.plot_runoff
        if sum(glaciers(j).runoff.q) > 0 % glaciers without associated runoff have all zeros in their discharge time series            
            [~]  = scatter(ax3,glaciers(j).runoff.x,glaciers(j).runoff.y,25,glacier_color,'filled','MarkerEdgeColor','black','MarkerFaceAlpha',1.0,'Marker','o');
            hp_rg = plot([glaciers(j).x,glaciers(j).runoff.x],[glaciers(j).y,glaciers(j).runoff.y],':','linewidth',1.5,'Color','black');
        end
        if datasets.opts.plot_glaciers
            hp_gf = plot([glaciers(j).x,fjord.x],[glaciers(j).y,fjord.y],':','linewidth',1.,'Color','black');
            [~]  = scatter(ax3,glaciers(j).x,glaciers(j).y,30,glacier_color,'filled','MarkerEdgeColor','black','MarkerFaceAlpha',0.9,'Marker','^');
            % text(glaciers(j).x,glaciers(j).y,num2str(glaciers(j).glacierID),'fontsize',8,'color','white');
        end            
        
    end
end % for glaciers
end % isfield glaciers
end % restrict_to_tfjords
end % for fjords
if datasets.opts.restrict_to_fjords==0
    for i_glaciers=1:length(glaciers)
        glacier=glaciers(i_glaciers);
        glacier_color = glacier_line_color(glacier.regionID,:);

        if isfield(glacier,'runoff') && datasets.opts.plot_runoff
            if sum(glacier.runoff.q) > 0 % glaciers without associated runoff have all zeros in their discharge time series            
                [~]  = scatter(ax3,glacier.runoff.x,glacier.runoff.y,25,glacier_color,'filled','MarkerEdgeColor','black','MarkerFaceAlpha',1.0,'Marker','o');
                hp_rg = plot([glacier.x_ref,glacier.runoff.x],[glacier.y_ref,glacier.runoff.y],':','linewidth',1.5,'Color','black');
            end
        end
        if datasets.opts.plot_glaciers
            [~]  = scatter(ax3,glacier.x_ref,glacier.y_ref,40,glacier_color,'filled','MarkerEdgeColor','black','MarkerFaceAlpha',0.9,'Marker','^');
        end
    end

end
axis xy equal

linkaxes([ax1,ax2,ax3]);
ax1.Visible = 'off';
ax1.XTick = []; ax1.YTick = [];
set(ax1,'color','none','visible','off');
ax2.Visible = 'off';
ax2.XTick = []; ax2.YTick = [];
set(ax2,'color','none','visible','off');
ax3.Visible = 'off';
ax3.XTick = []; ax3.YTick = [];
set(ax3,'color','none','visible','off');


% TODO: add box around plot and a scale bar
box on;

% Dummy plots for the legend entry
hs_f = scatter(ax3,0,0,50,'filled','MarkerEdgeColor','black','MarkerFaceColor','black','Marker','square');
hs_g = scatter(ax3,0,0,25,'filled','MarkerEdgeColor','black','MarkerFaceColor','black','Marker','^');
hs_r = scatter(ax3,0,0,25,'filled','MarkerEdgeColor','black','MarkerFaceColor','black','Marker','o');
for i_reg=1:length(glacier_line_color)
    [h_reg] = scatter(ax3,0,0,50,glacier_line_color(i_reg,:),'filled','MarkerEdgeColor','black','MarkerFaceAlpha',0.9,'Marker','square');
    handles_regions = [handles_regions; h_reg];
end


set([ax1,ax2,ax3],'Position',[.05 .11 .685 .815]);
cb1 = colorbar(ax1);
ylabel(cb1,'Bedrock (m a.s.l.)','fontsize',14)
cb2 = colorbar(ax2);
ylabel(cb2,'Ice surface (m a.s.l.)','fontsize',14)


legend_handles = [hs_f; hs_g; hs_r; handles_regions];
legend_text    = {'Fjords','Glaciers','Runoff','SW','SE','CW','CE','NW','NE','NO'};
% hl = legend(legend_handles,legend_text,'Location','northwest','fontsize', 14);
hl = legend(legend_handles,legend_text,'fontsize', 14);
set(hl,'Box','on','Color',[1.,1.,1.], 'EdgeColor',[0.,0.,0.]) ;

cb1.Position = [.56 .53 .0675 .35];
cb2.Position = [.56 .165 .0675 .35];
hl.Position = [.12 .683 .075 .1];

%% return the handles to all plots
fig_handles.hi_bed = hi_bed;
fig_handles.hi_sfc = hi_sfc;
fig_handles.hs_r = hs_r;
fig_handles.hp_rg = hp_rg;
fig_handles.hs_g = hs_g;
if datasets.opts.restrict_to_fjords, fig_handles.hp_gf = hp_gf; end
fig_handles.hs_f = hs_f;
fig_handles.hl   = hl;
fig_handles.cb1 = cb1;
fig_handles.cb2 = cb2;
% fig_handles.cb3 = cb3;
end