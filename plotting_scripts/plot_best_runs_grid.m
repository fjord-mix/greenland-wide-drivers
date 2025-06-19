function hf_grid = plot_best_runs_grid(ensemble_yr,res_box_yr)

fsize=16;
% rmse_threshold = 0.5;
n_years   = length(ensemble_yr);
% n_layers  = ensemble_yr{1}(1,1).p.N;
lcolor    = lines(n_years);
hb_colors = [];
lbl_yrs   = {};


%% Summary figure
hf_grid = figure('Position',[50 50, 500 500]); 
hold on; box on;
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
                rmse_fjd(i_fjord) = min(res_box(i_fjord).rmse_tf(:,2),[],'omitnan');
                x_fjd(i_fjord) = int32(str2double(res_box(i_fjord).id));% - 64; % converting the fjord ID ('A' to 'N') to an integer
                y_fjd(i_fjord) = 2015+i_year;

                size_marker(i_fjord) = 0.5-rmse_fjd(i_fjord);
                
                % if best_fjord_params(i_fjord).rmse_t > max_rmse, max_rmse = best_fjord_params(i_fjord).rmse_t; end
                % if best_fjord_params(i_fjord).rmse_t < min_rmse, min_rmse = best_fjord_params(i_fjord).rmse_t; end

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
blgd = bubblelegend('RMSE (^oC)','Location','southeast','fontsize',fsize);
blgd.Title.FontWeight='normal';
blgd.LimitLabels = {'\geq 2','0'};
xlabel('Fjord')
ylabel('Year')
set(gca,'YTick',2015+[1:n_years])

% hl = legend(hb_colors,lbl_yrs,'location','southwest','fontsize',fsize);
% hl.Layout.Tile = 5;
blgd.Title.FontWeight='normal';
set(blgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
% set(hl.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
% set(gcf, 'color', 'none');    
% set(gca, 'color', 'none');
end