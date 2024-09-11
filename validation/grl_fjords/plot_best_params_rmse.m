function hf = plot_best_params_rmse(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,i_tgt_day)

if nargin < 6, i_tgt_day=1; end
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize    = 14;
n_fjords = length(fjord_IDs);
n_years  = length(fjord_model_yr);
n_params = length(param_names);
lcolor   = cmocean('thermal',n_years);
letters  = lower(char(65:1:65+n_params));
rmse_bnds = [999 -999];
all_rmse = {};

fjord_names = cell([1,length(fjord_IDs)]);
for i=1:length(fjord_names), fjord_names{i} = fjord_IDs(i); end

hf = figure('Name','Best parameters','Position',[40 40 1200 1000]);
ht = tiledlayout('flow');


h_yr = [];
lbl_years = cell([n_years,1]);
for i_yr=1:n_years

    fjord_model = fjord_model_yr{i_yr};
    ensemble    = ensemble_yr{i_yr};
    res_box     = res_box_yr{i_yr};
    n_fjord_runs = length(fjord_model);

    %% Finding the best combination of model parameters
    name_fields = {'best_t','best_s','best_2'};
    name_fields{2,1} = cell(size(fjord_model));
    best_fjord_params = struct(name_fields{:});
    
    for i_fjord=1:n_fjord_runs
    
        % find run with the smallest RMSE
        rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf + (1-w_rmse_t).*res_box(i_fjord).rmse_sf;
        [best_rmse_t,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
        [best_rmse_s,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
        [best_rmse_2,inds_best2]    = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
        
        for i_param=1:n_params
            best_fjord_params(i_fjord).best_t.(param_names{i_param}) = ensemble(i_fjord,inds_best_tf).p.(param_names{i_param});
            best_fjord_params(i_fjord).best_s.(param_names{i_param}) = ensemble(i_fjord,inds_best_sf).p.(param_names{i_param});
            best_fjord_params(i_fjord).best_2.(param_names{i_param}) = ensemble(i_fjord,inds_best2).p.(param_names{i_param});
            best_fjord_params(i_fjord).rmse_t = best_rmse_t;
            best_fjord_params(i_fjord).rmse_s = best_rmse_s;
            best_fjord_params(i_fjord).rmse_2 = best_rmse_2;
            if best_rmse_t < rmse_bnds(1), rmse_bnds(1) = best_rmse_t; end
            if best_rmse_t > rmse_bnds(2), rmse_bnds(2) = best_rmse_t; end
            all_rmse{end+1} = best_rmse_t;
        end
    
    end
    
    %% Plotting
    
    for i_param=1:n_params
        nexttile(i_param); hold on; box on;
        for i_fjord=1:n_fjord_runs
            % fjord_names{i_fjord} = res_box(i_fjord).name;
            % fjord_names{i_fjord} = res_box(i_fjord).id;

            x_var = int32(str2double(res_box(i_fjord).id));% - 64; % converting the fjord ID ('A' to 'N') to an integer
            % x_var = i_fjord;
            % x_var = ensemble(i_fjord,1).s.Qsg_max
            % x_var = ensemble(i_fjord,1).p.L.*1e-3;
            size_marker = 250 - abs(best_fjord_params(i_fjord).rmse_t).*20;
            if size_marker < 0, size_marker = 10; end
            if i_param==1 && i_fjord==1
                h1 = scatter(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,abs(best_fjord_params(i_fjord).rmse_t),'filled','^','MarkerFaceAlpha',.75);
                % h2 = scatter(x_var,best_fjord_params(i_fjord).best_s.(param_names{i_param}),150,lcolor(i_yr,:),'filled','v','MarkerFaceAlpha',.5);
            else
                scatter(x_var,best_fjord_params(i_fjord).best_t.(param_names{i_param}),size_marker,abs(best_fjord_params(i_fjord).rmse_t),'filled','^','MarkerFaceAlpha',.75);
                % scatter(x_var,best_fjord_params(i_fjord).best_s.(param_names{i_param}),150,lcolor(i_yr,:),'filled','v','MarkerFaceAlpha',.5);
            end
        end
        ylabel(param_names{i_param})
        ylim([0.5*min(range_params{i_param}) 1.1*max(range_params{i_param})])
    
        % if i_yr==1
        %     xlim([0 n_fjords+1])
        %     set(gca,'Xtick',0:1:n_fjords+1)
        %     xlabels = get(gca,'XTickLabels');
        %     xlabels(2:end-1) = fjord_names;
        %     xlabels{1} = ' '; xlabels{end} = ' ';
        %     set(gca,'XtickLabels',xlabels,'fontsize',fsize);text(0.02,0.95,sprintf("(%s)",letters(i_param)),'Units','normalized','fontsize',fsize);
        % end
        % xlabel('Fjord length (km)');
        % xlabel('Peak subglacial discharge (m^3s^{-1})');
        
        set(gca,'fontsize',fsize)
        % xtickangle(90);
        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            set(gca,'YScale','log')
            ylim([0.1*min(range_params{i_param}) 10*max(range_params{i_param})])
        end
        hline(range_params{i_param},'--','color',[0.75 0.75 0.75])
        % clim(rmse_bnds);
        clim(prctile(cell2mat(all_rmse),[25,75]));
        xlabel('Fjord','fontsize',16)
    end
    h_yr = [h_yr h1];
    lbl_years{i_yr} = num2str(2015+i_yr);
end
% colormap(cmocean('thermal'))
colormap('parula')
hc = colorbar('eastoutside');
ylabel(hc,'RMSE (^oC)')
% legend([h1,h2,h3],{'RMSE_T','RMSE_S','RMSE_{both}'},'fontsize',fsize,'Location','Northwest');
% legend([h1,h2],{'RMSE_T','RMSE_S'},'fontsize',fsize,'Location','Northwest');
% legend(h_yr,lbl_years,'fontsize',fsize,'Location','best');

