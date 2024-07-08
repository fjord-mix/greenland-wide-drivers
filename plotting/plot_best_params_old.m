function hf = plot_best_params_old(fjord_model,res_box,dims_ensemble,param_names,range_params,tgt_days,name_days,i_tgt_day)

if nargin < 8, i_tgt_day=1; end
n_fjord_runs = length(fjord_model);
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize=12;
n_days = length(tgt_days);
lcolor = lines(3+n_days);


%% Finding the best combination of model parameters
name_fields = {'inds_best_tf','inds_best_sf','inds_best2'};
name_fields{2,1} = cell(size(fjord_model));
best_fjord_params = struct(name_fields{:});

for i=1:n_fjord_runs

    % find run with the smallest RMSE
    if n_days > 1
        rmse_both = w_rmse_t.*res_box(i).rmse_tf + (1-w_rmse_t).*res_box(i).rmse_sf;
        [~,i_min_rmse_tf] = min(res_box(i).rmse_tf,[],'all','omitnan');
        [~,i_min_rmse_sf] = min(res_box(i).rmse_sf,[],'all','omitnan');
        [~,i_min_rmse]    = min(rmse_both,[],'all','omitnan');
        [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,id_best_tf] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse_tf);
        [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,id_best_sf] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse_sf);
        [i1_best,i2_best,i3_best,i4_best,id_best] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse);
    
        inds_best_tf = [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,id_best_tf];
        inds_best_sf = [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,id_best_sf];
        inds_best2   = [i1_best,i2_best,i3_best,i4_best,id_best];
    else
        rmse_both = w_rmse_t.*res_box(i).rmse_tf + (1-w_rmse_t).*res_box(i).rmse_sf;
        [~,i_min_rmse_tf] = min(squeeze(res_box(i).rmse_tf(:,:,:,:,i_tgt_day)),[],'all','omitnan');
        [~,i_min_rmse_sf] = min(squeeze(res_box(i).rmse_sf(:,:,:,:,i_tgt_day)),[],'all','omitnan');
        [~,i_min_rmse]    = min(squeeze(rmse_both(:,:,:,:,i_tgt_day)),[],'all','omitnan');
        [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf] = ind2sub([dims_ensemble(2:end)],i_min_rmse_tf);
        [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf] = ind2sub([dims_ensemble(2:end)],i_min_rmse_sf);
        [i1_best,i2_best,i3_best,i4_best] = ind2sub([dims_ensemble(2:end)],i_min_rmse);
    
        inds_best_tf = [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf];
        inds_best_sf = [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf];
        inds_best2   = [i1_best,i2_best,i3_best,i4_best];
    end

    best_fjord_params(i).inds_best_tf = inds_best_tf;
    best_fjord_params(i).inds_best_sf = inds_best_sf;
    best_fjord_params(i).inds_best2 = inds_best2;

end

%% Plotting

hf = figure('Name','Best parameters','Position',[40 40 1200 1000]);
tiledlayout('flow');
fjord_names = cell(size(fjord_model));
for i_param=1:length(param_names)
    nexttile; hold on; box on;
    for i_fjord=1:n_fjord_runs
        h1 = scatter(i_fjord,range_params{i_param}(best_fjord_params(i_fjord).inds_best_tf(i_param)),150,lcolor(n_days+1,:),'filled','^','MarkerFaceAlpha',.5);
        h2 = scatter(i_fjord,range_params{i_param}(best_fjord_params(i_fjord).inds_best_sf(i_param)),150,lcolor(n_days+2,:),'filled','v','MarkerFaceAlpha',.5);
        h3 = scatter(i_fjord,range_params{i_param}(best_fjord_params(i_fjord).inds_best2(i_param)),150,lcolor(n_days+3,:),'filled','o','MarkerFaceAlpha',.5);
        % fjord_names{i_fjord} = res_box(i_fjord).name;
        fjord_names{i_fjord} = res_box(i_fjord).id;
    end
    ylabel(param_names{i_param})
    xlim([0 n_fjord_runs+1])
    ylim([0.9*min(range_params{i_param}) 1.1*max(range_params{i_param})])
    set(gca,'Xtick',0:1:n_fjord_runs+1)
    xlabels = get(gca,'XTickLabels');
    xlabels(2:end-1) = fjord_names;
    xlabels{1} = ' '; xlabels{end} = ' ';
    set(gca,'XtickLabels',xlabels,'fontsize',fsize);
    % set(gca,'fontsize',fsize)
    % xtickangle(90);
    if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
        set(gca,'YScale','log')
    end

    if i_param==1
        legend([h1,h2,h3],{'RMSE_T','RMSE_S','RMSE_{both}'},'fontsize',fsize,'Location','Northwest');
    end
end
if n_days > 1
    nexttile; hold on; box on;
    for i_fjord=1:n_fjord_runs
        h1 = scatter(i_fjord,best_fjord_params(i_fjord).inds_best_tf(end),150,lcolor(n_days+1,:),'filled','^','MarkerFaceAlpha',.5);
        h2 = scatter(i_fjord,best_fjord_params(i_fjord).inds_best_sf(end),150,lcolor(n_days+2,:),'filled','v','MarkerFaceAlpha',.5);
        h3 = scatter(i_fjord,best_fjord_params(i_fjord).inds_best2(end),150,lcolor(n_days+3,:),'filled','o','MarkerFaceAlpha',.5);
        fjord_names{i_fjord} = res_box(i_fjord).id;
    end
    ylabel('Time')
    xlim([0 n_fjord_runs+1])
    set(gca,'Xtick',0:1:n_fjord_runs+1)
    xlabels = get(gca,'XTickLabels');
    xlabels(2:end-1) = fjord_names;
    xlabels{1} = ' '; xlabels{end} = ' ';
    set(gca,'XtickLabels',xlabels);
    % xtickangle(90);
    ylim([0 length(tgt_days)+1])
    set(gca,'Ytick',0:1:length(tgt_days)+1)
    ylabels = get(gca,'YTickLabels');
    ylabels(2:end-1) = name_days;
    ylabels{1} = ' '; ylabels{end} = ' ';
    set(gca,'ytickLabels',ylabels);
    set(gca,'fontsize',fsize)
end

% fprintf("Best parametres for (%s) %s: \n",res_box(i).id,res_box(i).name)
%     fprintf("Param\t Temperature\t Salinity\t Both\n")
%     for i_param=1:length(param_names)
%         fprintf("%s \t %.1e\t",param_names{i_param},range_params{i_param}(inds_best_tf(i_param)))
%         fprintf("%.1e\t\t",range_params{i_param}(inds_best_sf(i_param)))
%         fprintf("%.1e\n",range_params{i_param}(inds_best2(i_param)))
%     end
%     fprintf("day \t %s\t\t%s\t\t%s\n",name_days{inds_best_tf(end)},name_days{inds_best_sf(end)},name_days{inds_best2(end)})
%     disp("============================")