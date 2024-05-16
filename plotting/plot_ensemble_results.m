function [hf_t,hf_s] = plot_ensemble_results(tgt_day,res_box,res_obs,fjord_model,mitgcm,dims_ensemble,param_names,range_params,n_fjord_runs)

    lcolor = lines(3);
    fsize=12;
    w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
    hf_t = figure('Name','Temperature comparison','Position',[40 40 1200 400*length(fjord_model)]);
    tiledlayout(length(fjord_model),2);

    hf_s = figure('Name','Salinity comparison','Position',[40 40 1200 400*length(fjord_model)]);
    tiledlayout(length(fjord_model),2);
    for i=1:n_fjord_runs
    
        % find run with the smallest RMSE
        rmse_both = w_rmse_t.*res_box(i).rmse_tf + (1-w_rmse_t).*res_box(i).rmse_sf;
        [~,i_min_rmse_tf] = min(res_box(i).rmse_tf,[],'all','omitnan');
        [~,i_min_rmse_sf] = min(res_box(i).rmse_sf,[],'all','omitnan');
        [~,i_min_rmse]    = min(rmse_both,[],'all','omitnan');
        [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,i5_best_tf] = ind2sub(dims_ensemble(2:end),i_min_rmse_tf);
        [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,i5_best_sf] = ind2sub(dims_ensemble(2:end),i_min_rmse_sf);
        [i1_best,i2_best,i3_best,i4_best,i5_best] = ind2sub(dims_ensemble(2:end),i_min_rmse);
    
        tf_best = res_box(i).ensemble_tf(:,i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,i5_best_tf);
        sf_best = res_box(i).ensemble_sf(:,i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,i5_best_sf);

        tf_best2 = res_box(i).ensemble_tf(:,i1_best,i2_best,i3_best,i4_best,i5_best);
        sf_best2 = res_box(i).ensemble_sf(:,i1_best,i2_best,i3_best,i4_best,i5_best);

        inds_best_tf = [i1_best_tf,i2_best_tf,i3_best_tf,i4_best_tf,i5_best_tf];
        inds_best_sf = [i1_best_sf,i2_best_sf,i3_best_sf,i4_best_sf,i5_best_sf];
        inds_best2   = [i1_best,i2_best,i3_best,i4_best,i5_best];
        fprintf("Best parametres for (%s) %s: \n",res_box(i).id,res_box(i).name)
        fprintf("Param\t Temperature\t Salinity\t Both\n")
        for i_param=1:length(param_names)
            fprintf("%s \t %.1e\t",param_names{i_param},range_params{i_param}(inds_best_tf(i_param)))
            fprintf("%.1e\t\t",range_params{i_param}(inds_best_sf(i_param)))
            fprintf("%.1e\n",range_params{i_param}(inds_best2(i_param)))
        end
        disp("============================")
        
        %% Plotting temperature
        figure(hf_t);
        nexttile; hold on; box on; grid on
        text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i).id,res_box(i).name, fjord_model(i).p.L/1e3),'units','normalized','fontsize',fsize)
        text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',fsize)
    
        % figure; hold on
        % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
        % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
        % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);
    
        hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
        hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
        hm = plot(mitgcm(i).Tprofile,-mitgcm(i).z,'linewidth',1.5,'color','k');
        plot(tf_best,-res_obs(i).zf,'linewidth',2.5,'color',lcolor(3,:),'linestyle',':');
        hb = plot(tf_best2,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    
        he = plot(res_box(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5]);
        plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        scatter(0,fjord_model(i).p.zgl,40,'v','filled','MarkerFaceColor','black')
        plot([0 0],[-fjord_model(i).p.H fjord_model(i).p.silldepth],'-k','linewidth',2)
        xlabel('Temperature (^oC)'); ylabel('Depth (m)');
        set(gca,'fontsize',fsize)
        xlim([-2 6])
        ylim([-fjord_model(i).p.H 0])
    
        nexttile; hold on; box on; grid on
        if length(res_box(i).t) == length(res_box(i).Tupper)
            hbl = plot(mitgcm(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k'); % the first two lines are dummy for the legend entries
            hml = plot(mitgcm(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k','LineStyle','--');
        
            hu = plot(mitgcm(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
            hi = plot(mitgcm(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
            hl = plot(mitgcm(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
        
            plot(mitgcm(i).t,mitgcm(i).Tupper,'linewidth',1.5,'color',lcolor(1,:),'LineStyle','--');
            plot(mitgcm(i).t,mitgcm(i).Tinter,'linewidth',1.5,'color',lcolor(2,:),'LineStyle','--');
            plot(mitgcm(i).t,mitgcm(i).Tlower,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
            % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
           
        end
        ylim([-2 5.5])
        ylabel('Temperature (^oC)'); xlabel('Model time (days)');
        set(gca,'fontsize',fsize)
        if i==1 
            hl1 = legend([hs, hf, he, hb, hm], {"Shelf","Fjord","Box model (ens)","Box model (best)","MITgcm"},'fontsize',fsize,'Location','Southeast'); 
            title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        end
        if i==1%n_fjords
            hl2 = legend([hbl,hml,hu,hi,hl],{"Box model","MITgcm","0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','Northeast'); 
            title(hl2,'Time series')
            hl2.NumColumns=1;
        end

        %% Plotting salinity
        figure(hf_s);
        nexttile; hold on; box on; grid on
        text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i).id,res_box(i).name, fjord_model(i).p.L/1e3),'units','normalized','fontsize',fsize)
        text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',fsize)
    
        hs = plot(res_obs(i).ss,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
        hf = plot(res_obs(i).sf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
        hm = plot(mitgcm(i).Sprofile,-mitgcm(i).z,'linewidth',1.5,'color','k');
        plot(sf_best,-res_obs(i).zf,'linewidth',2.5,'color',lcolor(3,:),'linestyle',':');
        hb = plot(sf_best2,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    
        he = plot(res_box(i).sf,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5]);
        plot(res_box(i).sfmin,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        plot(res_box(i).sfmax,-res_obs(i).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        scatter(35,fjord_model(i).p.zgl,40,'v','filled','MarkerFaceColor','black')
        plot([35 35],[-fjord_model(i).p.H fjord_model(i).p.silldepth],'-k','linewidth',2)
        xlabel('Salinity'); ylabel('Depth (m)');
        set(gca,'fontsize',fsize)
        xlim([30.5 35.5])
        ylim([-fjord_model(i).p.H 0])
    
        nexttile; hold on; box on; grid on
        if length(res_box(i).t) == length(res_box(i).Tupper)
            hbl = plot(mitgcm(i).t,res_box(i).Supper,'linewidth',1.5,'color','k'); % the first two lines are dummy for the legend entries
            hml = plot(mitgcm(i).t,res_box(i).Supper,'linewidth',1.5,'color','k','LineStyle','--');
        
            hu = plot(mitgcm(i).t,res_box(i).Supper,'linewidth',1.5,'color',lcolor(1,:));
            hi = plot(mitgcm(i).t,res_box(i).Sinter,'linewidth',1.5,'color',lcolor(2,:));
            hl = plot(mitgcm(i).t,res_box(i).Slower,'linewidth',1.5,'color',lcolor(3,:));
        
            plot(mitgcm(i).t,mitgcm(i).Supper,'linewidth',1.5,'color',lcolor(1,:),'LineStyle','--');
            plot(mitgcm(i).t,mitgcm(i).Sinter,'linewidth',1.5,'color',lcolor(2,:),'LineStyle','--');
            plot(mitgcm(i).t,mitgcm(i).Slower,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');           
        end
        % ylim([-2 5.5])
        ylabel('Salinity'); xlabel('Model time (days)');
        set(gca,'fontsize',fsize)
        if i==1 
            hl1 = legend([hs, hf, he, hb, hm], {"Shelf","Fjord","Box model (ens)","Box model (best)","MITgcm"},'fontsize',fsize,'Location','Southeast'); 
            title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        end
        if i==1%n_fjords
            hl2 = legend([hbl,hml,hu,hi,hl],{"Box model","MITgcm","0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','Northeast'); 
            title(hl2,'Time series')
            hl2.NumColumns=1;
        end
    end
end
