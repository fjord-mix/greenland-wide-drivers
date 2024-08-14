function hf = plot_best_params_dist(fjord_IDs,fjord_model_yr,ensemble_yr,res_box_yr,param_names,range_params,i_tgt_day,which_dist)

if nargin < 6, i_tgt_day=1; end
if nargin < 7, which_dist='Normal'; end
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize    = 14;
n_years  = length(fjord_model_yr);
n_params = length(param_names);
lcolor   = cmocean('thermal',n_years);
param_entries = cell(1,n_years);

fjord_names = cell([1,length(fjord_IDs)]);
for i=1:length(fjord_names), fjord_names{i} = fjord_IDs(i); end

hf = figure('Name','Best parameters','Position',[40 40 1200 1000]);
tiledlayout('flow');


h_yr = [];
lbl_years = cell([n_years,1]);
for i_yr=1:n_years

    fjord_model = fjord_model_yr{i_yr};
    ensemble    = ensemble_yr{i_yr};
    res_box     = res_box_yr{i_yr};
    n_fjords = length(fjord_model);

    %% Finding the best combination of model parameters
    name_fields = {'best_t','best_s','best_2'};
    name_fields{2,1} = cell(size(fjord_model));
    best_fjord_params = struct(name_fields{:});
    param_entry = NaN([n_fjords,n_params]);
    
    for i_fjord=1:n_fjords
    
        % find run with the smallest RMSE
        rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf + (1-w_rmse_t).*res_box(i_fjord).rmse_sf;
        [~,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
        [~,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
        [~,inds_best2]    = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
        
        for i_param=1:n_params
            best_fjord_params(i_fjord).best_t.(param_names{i_param}) = ensemble(i_fjord,inds_best_tf).p.(param_names{i_param});
            best_fjord_params(i_fjord).best_s.(param_names{i_param}) = ensemble(i_fjord,inds_best_sf).p.(param_names{i_param});
            best_fjord_params(i_fjord).best_2.(param_names{i_param}) = ensemble(i_fjord,inds_best2).p.(param_names{i_param});
            param_entry(i_fjord,i_param) = best_fjord_params(i_fjord).best_t.(param_names{i_param});
        end
    
    end
    
    
    %% Plotting
    
    for i_param=1:n_params
        nexttile(i_param); hold on; box on;
        x_bnds = range_params{i_param};
        x_val = linspace(x_bnds(1),x_bnds(2),100);
        if isempty(which_dist)
            h1 = histogram(param_entry(:,i_param),10,'FaceColor',lcolor(i_yr,:),'FaceAlpha',0.5);
        else
            % kern_marginal = uq_KernelMarginals(param_entry(:,i_param),range_params{i_param});
            kern_marginal= fitdist(param_entry(:,i_param),which_dist);
            h1 = plot(x_val,pdf(kern_marginal,x_val),'linewidth',2,'color',lcolor(i_yr,:));
        end
        
        % xtickangle(90);
        if (max(range_params{i_param}) - min(range_params{i_param}) > 1e3) || max(range_params{i_param}) - min(range_params{i_param}) < 1e-3
            set(gca,'XScale','log')
            xlim([0.1*min(range_params{i_param}) 10*max(range_params{i_param})])
        end
        if i_yr==n_years
            xlabel(param_names{i_param})
            set(gca,'fontsize',fsize)
            % xlim(range_params{i_param})
            % vline(range_params{i_param},'--','color',[0.25 0.25 0.25])
            xlim(range_params{i_param});
        end
    end
    h_yr = [h_yr h1];
    lbl_years{i_yr} = sprintf('%d (n=%d)',2015+i_yr,n_fjords);
end

% TODO: add legend pertaining to different years
% legend([h1,h2,h3],{'RMSE_T','RMSE_S','RMSE_{both}'},'fontsize',fsize,'Location','Northwest');
% legend([h1,h2],{'RMSE_T','RMSE_S'},'fontsize',fsize,'Location','Northwest');
legend(h_yr,lbl_years,'fontsize',fsize,'Location','northwest');

