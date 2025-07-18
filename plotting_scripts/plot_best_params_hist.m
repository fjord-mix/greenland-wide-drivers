function hf = plot_best_params_hist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,res_obs_yr,param_names,param_units,range_params,i_tgt_day)

if nargin < 7, i_tgt_day=1; end
fsize    = 16;
n_fjords = length(fjord_IDs);
n_years  = length(fjord_model_yr);
n_params = length(param_names);
lcolor   = lines(n_years);
letters=lower(char(65:65+20));
i_letter = 1;

fjord_names = cell([1,length(fjord_IDs)]);
for i=1:length(fjord_names), fjord_names{i} = fjord_IDs(i); end

hf = figure('Name','Best parameters','Position',[40 40 1000 300]);
ht = tiledlayout(1,n_params,'TileSpacing','tight','Padding','compact');

param_entries_all = cell([1,n_params]);
h_yr = [];
lbl_years = cell([n_years,1]);
% for i_yr=1:n_years

%% find the max/min T, S of the entire ensemble so we can normalise the data
max_t = 0;
min_t = 99;
max_s = 0;
min_s = 99;
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    res_obs  = res_obs_yr{i_year};
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            % if max(res_box(i_fjord).rmse_tf(:,2),[],'omitnan') > max_rmse_t
            %     max_rmse_t = max(res_box(i_fjord).rmse_tf(:,2),[],'omitnan');
            % end
            % if max(res_box(i_fjord).rmse_sf(:,2),[],'omitnan') > max_rmse_s
            %     max_rmse_s = max(res_box(i_fjord).rmse_sf(:,2),[],'omitnan');
            % end
            if max(res_obs(i_fjord).tf,[],'omitnan') > max_t
                max_t = max(res_obs(i_fjord).tf,[],'omitnan');
            end
            if min(res_obs(i_fjord).tf,[],'omitnan') < min_t
                min_t = min(res_obs(i_fjord).tf,[],'omitnan');
            end
            if max(res_obs(i_fjord).sf,[],'omitnan') > max_s
                max_s = max(res_obs(i_fjord).sf,[],'omitnan');
            end
            if min(res_obs(i_fjord).sf,[],'omitnan') < min_s
                min_s = min(res_obs(i_fjord).sf,[],'omitnan');
            end
        end
    end
end

%% Run loop
for i_yr=n_years:-1:1
    fjord_model = fjord_model_yr{i_yr};
    ensemble    = ensemble_yr{i_yr};
    res_box     = res_box_yr{i_yr};
    n_fjord_runs = length(fjord_model);

    %% Finding the best combination of model parameters
    name_fields = {'best_t','best_s','best_2'};
    name_fields{2,1} = cell(size(fjord_model));
    best_fjord_params = struct(name_fields{:});
    param_entry = NaN([n_fjords,n_params]);
    
    for i_fjord=1:n_fjord_runs

        % Filter out runs that have a high RMSE(t), i.e., we want to focus on the fjords we can simulate well
        % rmse_tf_filtered = res_box(i_fjord).rmse_tf;
        w_rmse_t  = 0.5;
        % z_rmse_t  = normalize(res_box(i_fjord).rmse_tf(:,2),"zscore");
        % z_rmse_s  = normalize(res_box(i_fjord).rmse_sf(:,2),"zscore");
        z_rmse_t  = res_box(i_fjord).rmse_tf(:,2)./(max_t-min_t);
        z_rmse_s  = res_box(i_fjord).rmse_sf(:,2)./(max_s-min_s);
        rmse_both = w_rmse_t * z_rmse_t + (1-w_rmse_t) * z_rmse_s;
        rmse_ts_filtered = rmse_both;

        rmse_ts_threshold = 1.0; % 1.0 won't filter anything
        rmse_ts_filtered(rmse_ts_filtered>rmse_ts_threshold) = NaN;
        [best_rmse_t,inds_best_tf] = min(squeeze(rmse_ts_filtered),[],'all','omitnan');

        % % find run with the smallest RMSE(sigma)
        % rmse_df_filtered = res_box(i_fjord).rmse_df;
        % rmse_df_threshold = 0.5; %prctile(rmse_tf_filtered(:,i_tgt_day),10,1);
        % rmse_df_filtered(rmse_df_filtered>rmse_df_threshold) = NaN;

        % find run with the smallest RMSE
        % rmse_tf_threshold = 0.5;
        % rmse_tf_filtered(rmse_ts_filtered>rmse_tf_threshold) = NaN;
        % [best_rmse_t,inds_best_tf] = min(squeeze(rmse_tf_filtered(:,i_tgt_day)),[],'all','omitnan');

        if ~isnan(best_rmse_t) % we only compute it if that fjord has any runs below RMSE = 0.5
            for i_param=1:n_params
                best_fjord_params(i_fjord).best_t.(param_names{i_param}) = ensemble(i_fjord,inds_best_tf).p.(param_names{i_param});
                best_fjord_params(i_fjord).rmse_t = best_rmse_t;
                param_entry(i_fjord,i_param) = best_fjord_params(i_fjord).best_t.(param_names{i_param});
            end
        end
    end
    
    %% Plotting
    
    for i_param=1:n_params
        % i_scat = i_param+(i_param-1)*2;

        ax1 = nexttile(i_param); hold on; box on;

        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            set(gca,'XScale','log')
            if min(range_params{i_param}) == 0
                lower_bnd=1;
            else
                lower_bnd=min(range_params{i_param});
            end
            xlim([lower_bnd max(range_params{i_param})])
            x_ticks = get(gca,'xTick');
        else
            xlim([min(range_params{i_param}) max(range_params{i_param})])
            x_ticks = get(gca,'xTick');
        end
        x_lims = get(gca,'XLim');
        

        % Histogram
        
        % Matlab's histogram function might act funny when plotted on a log scale if we do not treat the bin edges
        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            binedges = logspace(log10(x_lims(1)),log10(x_lims(end)),25);
        else
            binedges = linspace(x_lims(1),x_lims(end),25);
        end
        param_entry_filtered = param_entry(:,i_param);
        param_entries_all{i_param} = [param_entries_all{i_param}; param_entry_filtered];

        % [n, edges] = histcounts(param_entry_filtered,1e2,'Normalization','probability');
        % edge_centre = 0.5*(edges(1:end-1)+edges(2:end));
        % edge_centre(end) = edge_centre(end).*1.1;
        
        % h = histogram(param_entry_filtered,binedges,'Normalization','probability','FaceColor',lcolor(i_yr,:),'FaceAlpha',0.5,'Orientation','horizontal');
        % bar(edge_centre, n, 'barwidth', 1,'FaceColor',lcolor(i_yr,:),'FaceAlpha',0.5,'Horizontal','on');
        % histogram(param_entry_filtered,30,'FaceColor',lcolor(i_yr,:),'FaceAlpha',0.5,'Orientation','horizontal');
        
        % set(gca,'XTick',x_ticks)
        % set(gca,'YTickLabel',[],'XTickLabel',[])
        % xlim(x_lims)
        % if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
        %     set(gca,'XScale','log')
        % end
        if i_yr==1
            % colororder({'k','k'})
            % [n, edges] = histcounts(param_entries_all{i_param},1e2,'Normalization','probability');
            % edge_centre = 0.5*(edges(1:end-1)+edges(2:end));
            % edge_centre(end) = edge_centre(end).*1.1;

            yyaxis(ax1,'left')
            hhist = histogram(param_entries_all{i_param},binedges,'Normalization','cumcount','FaceAlpha',0.25,'Orientation','vertical');
            % histogram(param_entries_all{i_param},binedges,'Normalization','cdf','FaceAlpha',0.25,'Orientation','vertical');
            set(gca,'YColor','k');
            if i_param == 1
                ylabel('Cumulative','Color','k')
                % set(gca,'YColor','k');
            else
                set(gca,'YTickLabel',{})
            end
            
            yyaxis(ax1,'right')
            hhist = histogram(param_entries_all{i_param},binedges,'Normalization','count','FaceAlpha',0.5,'Orientation','vertical');
            % histogram(param_entries_all{i_param},binedges,'Normalization','probability','FaceAlpha',0.5,'Orientation','vertical');
            set(gca,'YColor','k');
            if i_param == n_params
                ylabel('Count','Color','k')
                set(gca,'YColor','k');
            else
                set(gca,'YTickLabel',{})
            end

            if i_param==1
                hl = legend('Cumulative','Count','fontsize',fsize,'Location','west');
            end

            text(0.02,1.0,['(',letters(i_letter),')'],'HorizontalAlignment','left','VerticalAlignment','top','Units','normalized','fontsize',fsize)
            i_letter=i_letter+1;
            set(gca,'fontsize',fsize)
        end
        
        if strcmp(param_names{i_param},'A0')
            x_ticks = [1e0, 1e2, 1e4, 1e6, 1e8];
        elseif strcmp(param_names{i_param},'C0')
            x_ticks = [1e2, 1e3, 1e4, 1e5];
        end
        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            xt_labels = cellstr(num2str(round(log10(x_ticks(:))), '10^%d'));
            set(gca,'XTick',x_ticks)
            set(gca,'XTickLabel',xt_labels)
        end
        xlabel([param_names{i_param},' (',param_units{i_param},')'])
        
    end
    % h_yr = [h_yr h1];
    % lbl_years{i_yr} = num2str(2015+i_yr);
end
% hl = legend(flip(h_yr),lbl_years,'fontsize',fsize,'Location','southeast');
% hl.NumColumns = 2;

% nexttile(1,[1,2]); hold on; box on;
% blgd = bubblelegend('RMSE (^oC)','Location','southwest','fontsize',fsize-4);
% blgd.Title.FontWeight='normal';
% blgd.LimitLabels = {'\geq 2','0'};
end