function [hf_t,hf_s] = plot_sensitivity_to_param(res_box,res_obs,fjord_model,param_names,param_ranges)


hf_t = figure('Name','Temperature sensitivity','Position',[40 40 900 250*length(param_names)]);
ht_t = tiledlayout('flow');
hf_s = figure('Name','Salinity sensitivity','Position',[40 40 900 250*length(param_names)]);
ht_s = tiledlayout('flow');
handle_fjords = [];
handle_fjords_s = [];
lcolor = lines(length(res_box));
lbl_fjords = cell([1,length(res_box)]);
for i_param=1:length(param_names)

    %% Plotting temperature
    figure(hf_t)
    nexttile; hold on; box on
    for i_fjord=1:length(res_box)
        ensemble_tf = res_box(i_fjord).ensemble_tf;
        if isempty(ensemble_tf), continue; end % we skip any empty entries

        % set up an array with only the OTHER parameter-related dimensions, which will be averaged
        dims_to_avg = [];
        for i_dim=2:length(size(ensemble_tf))
            if i_dim~=i_param+1, dims_to_avg = [dims_to_avg, i_dim]; end
        end
    
        % average the ensemble over the other dimensions
        tf = squeeze(mean(ensemble_tf,dims_to_avg,'omitnan'));
        tfmean = mean(tf,2,'omitnan');
        tfmin = min(tf,[],2,'omitnan');
        tfmax = max(tf,[],2,'omitnan');
    
        % if i_param==1
        %     text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i_fjord).id,res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','fontsize',14)
        % end
        if i_fjord==1
            text(0.02,0.95,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_ranges{i_param}(1),param_ranges{i_param}(end)),'units','normalized','fontsize',14)
        end
        havg = plot(tfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5]);
        hbnd = plot(tfmin,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        plot(tfmax,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');

        plot(tfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:));
        plot(tfmin,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:),'LineStyle','--');
        plot(tfmax,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:),'LineStyle','--');
        hfjd = plot(tfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:));

        scatter(0-0.1*i_fjord,fjord_model(i_fjord).p.zgl,40,'v','filled','MarkerFaceColor',lcolor(i_fjord,:))
        plot([0-0.1*i_fjord 0-0.1*i_fjord],[-fjord_model(i_fjord).p.H fjord_model(i_fjord).p.silldepth],'-','linewidth',2,'color',lcolor(i_fjord,:))
    
        set(gca,'fontsize',14)
        xlim([-1 2.5])
        % ylim([-fjord_model(i_fjord).p.H 0])
        if i_param==1
            handle_fjords = [handle_fjords hfjd];
            lbl_fjords{i_fjord} = res_box(i_fjord).name;
        end
    end

    %% Plotting salinity
    figure(hf_s)
    nexttile; hold on; box on
    for i_fjord=1:length(res_box)
        ensemble_sf = res_box(i_fjord).ensemble_sf;
        if isempty(ensemble_sf), continue; end % we skip any empty entries

        % set up an array with only the OTHER parameter-related dimensions, which will be averaged
        dims_to_avg = [];
        for i_dim=2:length(size(ensemble_sf))
            if i_dim~=i_param+1, dims_to_avg = [dims_to_avg, i_dim]; end
        end
    
        % average the ensemble over the other dimensions
        sf = squeeze(mean(ensemble_sf,dims_to_avg,'omitnan'));
        sfmean = mean(sf,2,'omitnan');
        sfmin = min(sf,[],2,'omitnan');
        sfmax = max(sf,[],2,'omitnan');
    
        if i_fjord==1
            text(0.98,0.95,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_ranges{i_param}(1),param_ranges{i_param}(end)),'units','normalized','fontsize',14,'HorizontalAlignment','right')
        end
        havg_s = plot(sfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5]);
        hbnd_s = plot(sfmin,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');
        plot(sfmax,-res_obs(i_fjord).zf,'linewidth',1.5,'color',[0.5 0.5 0.5],'LineStyle','--');

        plot(sfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:));
        plot(sfmin,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:),'LineStyle','--');
        plot(sfmax,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:),'LineStyle','--');
        hfjd_s = plot(sfmean,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(i_fjord,:));

        scatter(35-0.1*i_fjord,fjord_model(i_fjord).p.zgl,40,'v','filled','MarkerFaceColor',lcolor(i_fjord,:))
        plot([35-0.1*i_fjord 35-0.1*i_fjord],[-fjord_model(i_fjord).p.H fjord_model(i_fjord).p.silldepth],'-','linewidth',2,'color',lcolor(i_fjord,:))
    
        set(gca,'fontsize',14)
        xlim([33 35])
        % ylim([-fjord_model(i_fjord).p.H 0])
        if i_param==1
            handle_fjords_s = [handle_fjords_s hfjd_s];
            lbl_fjords{i_fjord} = res_box(i_fjord).name;
        end
    end

end
handle_fjords = [handle_fjords havg hbnd];
handle_fjords_s = [handle_fjords_s havg_s hbnd_s];
lbl_fjords = lbl_fjords(~cellfun(@isempty,lbl_fjords));
lbl_fjords{end+1} = 'Mean';
lbl_fjords{end+1} = 'Spread';

figure(hf_t)
nexttile; axis off
hl = legend(gca,handle_fjords,lbl_fjords,'fontsize',14,'Location','North');
% hl.Layout.Tile='southeast';
xlabel(ht_t,'Temperature (^oC)','fontsize',14);
ylabel(ht_t,'Depth (m)','fontsize',14);
ht_t.TileSpacing='compact';
ht_t.Padding='compact';

figure(hf_s)
nexttile; axis off
hl = legend(gca,handle_fjords_s,lbl_fjords,'fontsize',14,'Location','North');
% hl.Layout.Tile='southeast';
xlabel(ht_s,'Salinity','fontsize',14);
ylabel(ht_s,'Depth (m)','fontsize',14);
ht_s.TileSpacing='compact';
ht_s.Padding='compact';

end