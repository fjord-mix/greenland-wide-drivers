function hf = plot_misfits_per_parameter(X,ensemble_yr,res_box_yr,param_names,which_fjords,i_tgt_day)

if nargin < 6, i_tgt_day=2; end
n_years  = length(res_box_yr);
n_params = length(param_names);
lcolor = lines(n_years);
fsize = 16;

hf = figure('Name','Mistfit per parameter','Position',[40 40 1200 500]);
ht = tiledlayout(1,length(param_names));

for i_yr = 1:n_years
    ensemble = ensemble_yr{i_yr};
    res_box  = res_box_yr{i_yr};

    %% Getting bins
    % 10 bins of percentiles
    bins       = 0:10:100;
    param_bnds_prctile = NaN([length(bins),n_params]);
    for i_param=1:n_params
        param_bnds_prctile(:,i_param) = prctile(X(:,i_param),bins);
    end
    
    % loops over all ensemble entries
    struct_fields = param_names;
    struct_fields{2,1} = cell(size(ensemble));
    mask_bnds_prctile = struct(struct_fields{:});
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).p)
                for i_param=1:n_params
                    % for each parameter, assigns a mask of 1:10 if the value is is in the 1:100th percentile
                    for i_bin=1:length(bins)-1
                        if ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds_prctile(i_bin,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds_prctile(i_bin+1,i_param)
                            mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) = i_bin;
                        end
                    end
                end
            end
        end
    end

    %% Plotting
    for i_fjord = 1:length(res_box)
        if exist('which_fjords','Var') && ~isempty(which_fjords)
            for i_tgt_fjords=1:length(which_fjords)
                if strcmp(which_fjords{i_tgt_fjords},res_box(i_fjord).id) == 1
                    plot_fjord=1;
                    break
                else
                    plot_fjord=0;
                end
            end
        end
        if exist('plot_fjord','var') && ~plot_fjord,continue; end

        % rmse_arr = res_box(i_fjord).rmse_tf(:,i_tgt_day);
        % rmse_arr(rmse_arr > 4) = NaN;

        for i_param=1:n_params
            nexttile(i_param); hold on; box on;
            % param_arr = NaN(1,size(ensemble,2));
            % for i_run=1:size(ensemble,2)
            %     if ~isempty(ensemble(i_fjord,i_run).p)
            %         param_arr(i_run) = ensemble(i_fjord,i_run).p.(param_names{i_param});
            %     end
            % end
            % scatter(param_arr,rmse_arr,80,lcolor(i_yr,:),'filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
            % plot(param_arr,rmse_arr,'color',lcolor(i_yr,:),'linewidth',1.5)

            pctl_plot = NaN([1,length(bins)-1]);
            rmse_plot = NaN(size(pctl_plot));
            for i_bin=1:length(bins)-1
                rmse_bin = NaN(1,size(ensemble,2));
                % find profiles that fit into that bin for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) == i_bin
                        rmse_bin(i_run) = res_box(i_fjord).rmse_tf(i_run,i_tgt_day);
                    end
                end
                rmse_plot(i_bin) = mean(rmse_bin,'omitnan');
                pctl_plot(i_bin) = i_bin*10;
            end
            rmse_plot(rmse_plot > 10) = NaN;
            % scatter(pctl_plot,rmse_plot,50,'filled','MarkerFaceColor',lcolor(i_yr,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            plot(pctl_plot,rmse_plot,'o-','linewidth',1.5)

            if i_yr==n_years
                xlabel(sprintf("%s quantile",param_names{i_param}))
                xlim([0 110])
                set(gca,'fontsize',fsize)
            end
        end
    end
end
ylabel(ht,'Mean RMSE (^oC)','fontsize',fsize);
end