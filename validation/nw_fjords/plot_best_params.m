function hf = plot_best_params(fjord_model,ensemble,res_box,param_names,range_params,i_tgt_day)

if nargin < 6, i_tgt_day=1; end
n_fjord_runs = length(fjord_model);
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize=14;
lcolor = lines(3);


%% Finding the best combination of model parameters
name_fields = {'best_t','best_s','best_2'};
name_fields{2,1} = cell(size(fjord_model));
best_fjord_params = struct(name_fields{:});

for i_fjord=1:n_fjord_runs

    % find run with the smallest RMSE
    rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf + (1-w_rmse_t).*res_box(i_fjord).rmse_sf;
    [~,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
    [~,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
    [~,inds_best2]    = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
    
    for i_param=1:length(param_names)
        best_fjord_params(i_fjord).best_t.(param_names{i_param}) = ensemble(i_fjord,inds_best_tf).p.(param_names{i_param});
        best_fjord_params(i_fjord).best_s.(param_names{i_param}) = ensemble(i_fjord,inds_best_sf).p.(param_names{i_param});
        best_fjord_params(i_fjord).best_2.(param_names{i_param}) = ensemble(i_fjord,inds_best2).p.(param_names{i_param});
    end

end

%% Plotting

hf = figure('Name','Best parameters','Position',[40 40 1200 1000]);
tiledlayout('flow');
fjord_names = cell(size(fjord_model));
for i_param=1:length(param_names)
    nexttile; hold on; box on;
    for i_fjord=1:n_fjord_runs
        x_var = i_fjord;
        % x_var = ensemble(i_fjord,1).s.Qsg_max;
        % x_var = ensemble(i_fjord,1).p.L.*1e-3;
        h1 = scatter(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),250,lcolor(1,:),'filled','^','MarkerFaceAlpha',.5);
        h2 = scatter(x_var,best_fjord_params(i_fjord).best_s.(param_names{i_param}),250,lcolor(2,:),'filled','v','MarkerFaceAlpha',.5);
        
        % fjord_names{i_fjord} = res_box(i_fjord).name;
        fjord_names{i_fjord} = res_box(i_fjord).id;
    end
    ylabel(param_names{i_param})
    ylim([0.5*min(range_params{i_param}) 1.1*max(range_params{i_param})])

    xlim([0 n_fjord_runs+1])
    set(gca,'Xtick',0:1:n_fjord_runs+1)
    xlabels = get(gca,'XTickLabels');
    xlabels(2:end-1) = fjord_names;
    xlabels{1} = ' '; xlabels{end} = ' ';
    set(gca,'XtickLabels',xlabels);
    % xlabel('Peak subglacial discharge (m^3s^{-1})','fontsize',fsize);
    % xlabel('Fjord length (km)');
    
    
    set(gca,'fontsize',fsize)
    % xtickangle(90);
    if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
        set(gca,'YScale','log')
        ylim([0.1*min(range_params{i_param}) 10*max(range_params{i_param})])
    end
    hline(range_params{i_param},'--','color',[0.75 0.75 0.75])
    if i_param==1
        % legend([h1,h2,h3],{'RMSE_T','RMSE_S','RMSE_{both}'},'fontsize',fsize,'Location','Northwest');
        legend([h1,h2],{'RMSE_T','RMSE_S'},'fontsize',fsize,'Location','Northwest');
    end
end
