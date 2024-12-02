function [hf_profiles,hf_series] = plot_multiday_ensemble(fjord_model,res_box,res_obs,dims_ensemble,param_names,range_params,tgt_days,name_days,i_tgt_day,plt_series,verbose)

n_fjord_runs = length(fjord_model);
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize=12;

if nargin < 10 || isempty(plt_series),plt_series= 0; end
if nargin < 11 || isempty(verbose),   verbose   = 0; end

lcolor = lines(3+length(tgt_days));
hf_profiles = figure('Name','Temperature profiles','Position',[40 40 1200 400*length(fjord_model)/2]);
tiledlayout(length(fjord_model)/2,2);
if plt_series
    hf_series   = figure('Name','Temperature evolution','Position',[40 40 1200 400*length(fjord_model)/2]);
    tiledlayout(length(fjord_model)/2,2);
end
for i=1:n_fjord_runs

    % find run with the smallest RMSE
    if isempty(i_tgt_day)
        rmse_both = w_rmse_t.*res_box(i).rmse_tf + (1-w_rmse_t).*res_box(i).rmse_sf;
        [~,i_min_rmse_tf] = min(res_box(i).rmse_tf,[],'all','omitnan');
        [~,i_min_rmse_sf] = min(res_box(i).rmse_sf,[],'all','omitnan');
        [~,i_min_rmse]    = min(rmse_both,[],'all','omitnan');
        [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,id_best_tf] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse_tf);
        [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,id_best_sf] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse_sf);
        [i1_best,i2_best,i3_best,i4_best,id_best] = ind2sub([dims_ensemble(2:end),length(tgt_days)],i_min_rmse);
    
        tf_best = res_box(i).ensemble_tf(:,i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,id_best_tf);
        sf_best = res_box(i).ensemble_sf(:,i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,id_best_sf);
    
        tf_best2 = res_box(i).ensemble_tf(:,i1_best,i2_best,i3_best,i4_best,id_best);
        sf_best2 = res_box(i).ensemble_sf(:,i1_best,i2_best,i3_best,i4_best,id_best);
    
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
    
        tf_best = res_box(i).ensemble_tf(:,i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf);
        sf_best = res_box(i).ensemble_sf(:,i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf);
    
        tf_best2 = res_box(i).ensemble_tf(:,i1_best,i2_best,i3_best,i4_best);
        sf_best2 = res_box(i).ensemble_sf(:,i1_best,i2_best,i3_best,i4_best);
    
        inds_best_tf = [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf];
        inds_best_sf = [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf];
        inds_best2   = [i1_best,i2_best,i3_best,i4_best];
    end

    if verbose
        fprintf("Best parametres for (%s) %s: \n",res_box(i).id,res_box(i).name)
        fprintf("Param\t Temperature\t Salinity\t Both\n")
        for i_param=1:length(param_names)
            fprintf("%s \t %.1e\t",param_names{i_param},range_params{i_param}(inds_best_tf(i_param)))
            fprintf("%.1e\t\t",range_params{i_param}(inds_best_sf(i_param)))
            fprintf("%.1e\n",range_params{i_param}(inds_best2(i_param)))
        end
        fprintf("day \t %s\t\t%s\t\t%s\n",name_days{inds_best_tf(end)},name_days{inds_best_sf(end)},name_days{inds_best2(end)})
        disp("============================")
    end
    %% Plotting
    figure(hf_profiles)
    nexttile; hold on; box on; grid on
    text(0.02,1.02,sprintf("(%s) %s (%.0f km long)",res_box(i).id,res_box(i).name, fjord_model(i).p.L/1e3),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
    text(0.02,0.02,sprintf("n=%d",res_box(i).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)

    % figure; hold on
    % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
    % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
    % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);

    hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
    hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
    hbest = plot(tf_best,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    hbest2 = plot(tf_best2,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    hb=[];

    for i_day=1:length(tgt_days)
        hb_d = plot(res_box(i).tf(:,i_day),-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3+i_day,:));
        hb = [hb hb_d];
    end
    % plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    % plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    scatter(0,fjord_model(i).p.zgl,40,'v','filled','MarkerFaceColor','black')
    plot([0 0],[-fjord_model(i).p.H fjord_model(i).p.silldepth],'-k','linewidth',2)

    if mod(i,2) > 0, ylabel('Depth (m)'); end
    if i>n_fjord_runs-2
        xlabel('Temperature (^oC)');  
    else
        set(gca,'xticklabels',[])
    end
    set(gca,'fontsize',fsize)
    xlim([-2 6])
    ylim([-fjord_model(i).p.H 0])
    if i==1 
        string_legend = {"Shelf","Fjord","Best_T","Best_{TS}"};
        for i_day=1:length(tgt_days)
            string_legend{end+1} = sprintf("mean_{%s}",name_days{i_day});
        end
        hl1 = legend([hs, hf, hbest, hbest2, hb],string_legend,'fontsize',fsize,'Location','Southeast'); 
        % title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        hl1.NumColumns=2;
    end

    if plt_series
        figure(hf_series)
        nexttile; hold on; box on; grid on
        text(0.02,0.99,sprintf("(%s) %s (%.0f km long; n=%d)",res_box(i).id,res_box(i).name, fjord_model(i).p.L/1e3,res_box(i).n),'units','normalized','VerticalAlignment','top','fontsize',fsize)
        % text(0.02,0.02,sprintf("n=%d",res_box(i).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)
        if length(res_box(i).t) == length(res_box(i).Tupper)
            
            hu = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
            hi = plot(res_box(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
            hl = plot(res_box(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
        
            % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
           
            for i_day=1:length(tgt_days)
                xline(tgt_days(i_day),'linestyle','--','color',lcolor(3+i_day,:),'linewidth',2)
            end
        end
        ylim([-2 5.5])
        if mod(i,2) > 0, ylabel('Temperature (^oC)'); end
        if i>n_fjord_runs-2 
            xlabel('Model time (days)'); 
        else
            set(gca,'xticklabels',[]);
       end
        set(gca,'fontsize',fsize)
        if i==1%n_fjords
            hl2 = legend([hu,hi,hl],{"0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','Northeast'); 
            % title(hl2,'Time series')
            hl2.NumColumns=1;
        end
    end
end



end