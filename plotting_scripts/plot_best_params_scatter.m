function hf = plot_best_params_scatter(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,param_units,range_params,i_tgt_day)

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

hf = figure('Name','Best parameters','Position',[40 40 1000 500]);
ht = tiledlayout(1,n_params,'TileSpacing','tight','Padding','compact');

param_entries_all = cell([1,n_params]);
h_yr = [];
lbl_years = cell([n_years,1]);
% for i_yr=1:n_years
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

        % Filter out runs that have a high RMSE, i.e., we want to focus on the fjords we can simulate well
        rmse_tf_filtered = res_box(i_fjord).rmse_tf;
        rmse_tf_threshold = 0.5; %prctile(rmse_tf_filtered(:,i_tgt_day),10,1);
        rmse_tf_filtered(rmse_tf_filtered>rmse_tf_threshold) = NaN;
    
        % find run with the smallest RMSE
        [best_rmse_t,inds_best_tf] = min(squeeze(rmse_tf_filtered(:,i_tgt_day)),[],'all','omitnan');
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
        min_rmse = 10;
        max_rmse = 0;

        % scatter
        nexttile(i_param); hold on; box on;
        for i_fjord=1:n_fjord_runs
            x_var = int32(str2double(res_box(i_fjord).id));% - 64; % converting the fjord ID ('A' to 'N') to an integer
            % size_marker = 250 - abs(best_fjord_params(i_fjord).rmse_t).*20;
            if rmse_tf_threshold > 3
                size_marker = 2-best_fjord_params(i_fjord).rmse_t;
            else
                size_marker = rmse_tf_threshold-best_fjord_params(i_fjord).rmse_t;
            end
            % if best_fjord_params(i_fjord).rmse_t > max_rmse, max_rmse = best_fjord_params(i_fjord).rmse_t; end
            % if best_fjord_params(i_fjord).rmse_t < min_rmse, min_rmse = best_fjord_params(i_fjord).rmse_t; end

            % if size_marker < 0, size_marker = 10; end
            if ~isempty(best_fjord_params(i_fjord).best_t)
                if i_param==1 && i_fjord==1
                    % h1 = scatter(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,lcolor(i_yr,:),'filled','o','MarkerFaceAlpha',.5);
                    h1 = bubblechart(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,lcolor(i_yr,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                else
                    % scatter(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,lcolor(i_yr,:),'filled','o','MarkerFaceAlpha',.5);
                    bubblechart(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,lcolor(i_yr,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                end
            end
        end
        xlim([0 150])
        ylabel([param_names{i_param},' (',param_units{i_param},')'])
        ylim([0.5*min(range_params{i_param}) 1.1*max(range_params{i_param})])
    
        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            set(gca,'YScale','log')
            if min(range_params{i_param}) == 0
                lower_bnd=1;
            else
                lower_bnd=0.5.*min(range_params{i_param});
            end
            ylim([lower_bnd 7*max(range_params{i_param})])
        end
        % hline(range_params{i_param},'--','color',[0.75 0.75 0.75])
        plot([0 150],[range_params{i_param}(1),range_params{i_param}(1)],'linestyle','--','color',[0.75 0.75 0.75])
        plot([0 150],[range_params{i_param}(2),range_params{i_param}(2)],'linestyle','--','color',[0.75 0.75 0.75])
        y_lims = get(gca,'YLim');
        y_ticks = get(gca,'YTick');
        if i_yr==1
            text(0.02,1.0,['(',letters(i_letter),')'],'HorizontalAlignment','left','VerticalAlignment','top','Units','normalized','fontsize',fsize)
            i_letter=i_letter+1;
        end
        set(gca,'fontsize',fsize)
        xlabel('Fjord')
        bubblesize([1 15])
        bubblelim([0 0.2])
    end
    h_yr = [h_yr h1];
    lbl_years{i_yr} = num2str(2015+i_yr);
end

hl = legend(flip(h_yr),lbl_years,'fontsize',fsize,'Location','southeast');
hl.NumColumns = 2;

% nexttile(1,[1,2]); hold on; box on;
blgd = bubblelegend('RMSE (Kg m^{-3})','Location','southwest','fontsize',fsize-4);
blgd.Title.FontWeight='normal';
blgd.LimitLabels = {'\geq 0.2','0'};
end