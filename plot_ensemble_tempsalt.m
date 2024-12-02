function hf_ts = plot_ensemble_tempsalt(fjord_model,ensemble,res_box,res_obs,n_runs,tgt_days,i_tgt_day,which_fjords)

if isempty(tgt_days)
    n_days = 1;
else
    n_days = length(tgt_days);
end

n_fjord_runs = length(fjord_model);
w_rmse_t     = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)

fsize   = 14;
lcolor  = lines(3+length(tgt_days));
letters = lower(char(65:65+length(which_fjords)*2));

fig_width = 900; 
fig_height = 300*length(which_fjords);
% rmse_table(size(fjord_model)) = struct("tf_rpm",[],"sf_rpm",[],"ts_rpm",[],"tf_gcm",[],"sf_gcm",[],"ts_gcm",[]);
% rmse_table = cell(size(fjord_model));

hf_ts    = figure('Name','Temperature and salinity profiles','Position',[40 40 fig_width fig_height]);
ht       = tiledlayout(length(which_fjords),2,'TileSpacing','loose','padding','compact');
i_iter   = 0;
i_letter =1;
for i_fjord=1:n_fjord_runs

    % find run with the smallest RMSE
    if isempty(res_box(i_fjord).rmse_tf)
        continue
    end
    [~,tf_best,sf_best] = get_best_profiles_rmse(res_box,i_fjord,n_runs,w_rmse_t,tgt_days,i_tgt_day);
    
    % if isempty(res_box(i_fjord).rmse_tf)
    %     continue
    % elseif isempty(i_tgt_day)
    %     rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf./mean(res_box(i_fjord).tf,'omitnan') + (1-w_rmse_t).*res_box(i_fjord).rmse_sf./mean(res_box(i_fjord).sf,'omitnan');
    % 
    %     [rmse_table(i_fjord).tf_rpm,i_min_rmse_tf] = min(res_box(i_fjord).rmse_tf,[],'all','omitnan');
    %     [rmse_table(i_fjord).sf_rpm,i_min_rmse_sf] = min(res_box(i_fjord).rmse_sf,[],'all','omitnan');
    %     [rmse_table(i_fjord).ts_rpm,i_min_rmse]    = min(rmse_both,[],'all','omitnan');
    % 
    %     [irun_best_tf,id_best_tf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_tf);
    %     [irun_best_sf,id_best_sf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_sf);
    %     [irun_best,id_best] = ind2sub([n_runs,length(tgt_days)],i_min_rmse);
    % 
    %     tf_best = res_box(i_fjord).ensemble_tf(:,irun_best_tf,id_best_tf);
    %     sf_best = res_box(i_fjord).ensemble_sf(:,irun_best_sf,id_best_sf);
    % 
    %     tf_best2 = res_box(i_fjord).ensemble_tf(:,irun_best,id_best);
    %     sf_best2 = res_box(i_fjord).ensemble_sf(:,irun_best,id_best);
    % 
    %     inds_best_tf = [irun_best_tf,id_best_tf];
    %     inds_best_sf = [irun_best_sf,id_best_sf];
    %     inds_best2   = [irun_best,id_best];
    % else
    %     rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf./mean(res_box(i_fjord).tf,1,'omitnan') + (1-w_rmse_t).*res_box(i_fjord).rmse_sf./mean(res_box(i_fjord).sf,1,'omitnan');
    % 
    %     [rmse_table(i_fjord).tf_rpm,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
    %     [rmse_table(i_fjord).sf_rpm,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
    %     [rmse_table(i_fjord).ts_rpm,inds_best2]   = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
    % 
    %     tf_best = res_box(i_fjord).ensemble_tf(:,inds_best_tf,i_tgt_day);
    %     sf_best = res_box(i_fjord).ensemble_sf(:,inds_best_sf,i_tgt_day);
    % 
    %     tf_best2 = res_box(i_fjord).ensemble_tf(:,inds_best2,i_tgt_day);
    %     sf_best2 = res_box(i_fjord).ensemble_sf(:,inds_best2,i_tgt_day);
    % end
    % 
    % % adds Shelf-fjord RMSE to the table
    % rmse_table(i_fjord).rmse_ts = res_box(i_fjord).rmse_ts;
    % rmse_table(i_fjord).rmse_ss = res_box(i_fjord).rmse_ss;


    %% Plotting temperature
    if nargin > 7
        plot_fjord=1;
        for i_tgt_fjords=1:length(which_fjords)
            if strcmp(which_fjords{i_tgt_fjords},res_box(i_fjord).id) == 1
                plot_fjord=1;
                break
            else
                plot_fjord=0;
            end
        end
        if ~plot_fjord, continue; end
    end
    i_iter = i_iter+1;
    nexttile; hold on; box on; grid on
    text(0.98,0.98,sprintf("(%s)",letters(i_letter)),'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right','FontSize',fsize)
    text(0.02,1.02,sprintf("%s) %s (%.0f km)",res_box(i_fjord).id,res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
    % text(0.98,0.02,sprintf("n=%.1f %%",res_box(i_fjord).n),'Units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',fsize)
    i_letter=i_letter+1;

    % Observed fjord profile
    hf = plot(res_obs(i_fjord).tf,-res_obs(i_fjord).zf,'linewidth',2.5,'color',lcolor(2,:));
    hb=[];
 
    % RPM profile(s)
    % forcing
    if isfield(res_box(i_fjord),'Tforc') % we use this in case we are not using the observed profile as the shelf forcing
        hs = plot(res_box(i_fjord).Tforc,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(1,:),'LineStyle','-');
    else
        hs = plot(res_obs(i_fjord).ts,-res_obs(i_fjord).zs,'linewidth',2.5,'color',lcolor(1,:)); % plots observed if Tforc doesnt exist
    end
    % results
    for i_day=1:n_days
        % res_box(i_fjord).zf = res_box(i_fjord).zf';

        y2 = [-res_box(i_fjord).zf; flip(-res_box(i_fjord).zf)];
        inBetween = [res_box(i_fjord).tfmin(:,i_day); flip(res_box(i_fjord).tfmax(:,i_day))];
        try
            hfjd = patch(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
            hfjd = [];
        end
        try
            inBetween = [res_box(i_fjord).tf1sl(:,i_day); flip(res_box(i_fjord).tf1su(:,i_day))];
            fill(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
            hfjd = [];
        end
        hb_d = plot(res_box(i_fjord).tf(:,i_day),-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3+i_day,:));
        hb = [hb hb_d];
    end

    % Best run
    hbest = plot(tf_best,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3,:));
    % hbest2 = plot(tf_best2,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3,:),'LineStyle','--');
    
    % Depiction of GL and sill
    scatter(0,-fjord_model(i_fjord).p.Hgl,60,'v','filled','MarkerFaceColor','black')
    if fjord_model(i_fjord).p.sill
        plot([0 0],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)
    end

    set(gca,'fontsize',fsize)
    % xlim([-2.5 9])
    xlim([-4 7.5])
    ylim([-fjord_model(i_fjord).p.H 0])
    if i_iter==length(which_fjords)
        xlabel('Temperature (^oC)','fontsize',fsize+2);  
    end
    
    %% Plotting salinity    
    nexttile; hold on; box on; grid on
    text(0.98,0.98,sprintf("(%s)",letters(i_letter)),'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right','FontSize',fsize)
    i_letter=i_letter+1;
    % Observed shelf and fjord profiles
    
    plot(res_obs(i_fjord).sf,-res_obs(i_fjord).zf,'linewidth',2.5,'color',lcolor(2,:));        

    % forcing
    if isfield(res_box(i_fjord),'Sforc')
        plot(res_box(i_fjord).Sforc,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(1,:),'LineStyle','-');
    else
        plot(res_obs(i_fjord).ss,-res_obs(i_fjord).zs,'linewidth',2.5,'color',lcolor(1,:));
    end
    for i_day=1:n_days
        y2 = [-res_box(i_fjord).zf; flip(-res_box(i_fjord).zf)];
        inBetween = [res_box(i_fjord).sfmin(:,i_day); flip(res_box(i_fjord).sfmax(:,i_day))];
        try
            patch(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
        end
        try
            inBetween = [res_box(i_fjord).sf1sl(:,i_day); flip(res_box(i_fjord).sf1su(:,i_day))];
            fill(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
        end
        plot(res_box(i_fjord).sf(:,i_day),-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3+i_day,:));
    end

    % Best run
    plot(sf_best,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3,:),'LineStyle','-');
    % plot(sf_best2,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(3,:));

    % Depiction of GL and sill
    scatter(33,-fjord_model(i_fjord).p.Hgl,60,'v','filled','MarkerFaceColor','black')
    if fjord_model(i_fjord).p.sill
        plot([33 33],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)
    end

    set(gca,'fontsize',fsize)
    xlim([30.5 35])
    ylim([-fjord_model(i_fjord).p.H 0])

    %% Legend
    if i_iter==1 
        % string_legend = {"Shelf","Fjord","Best_T","Best_{TS}"};
        if exist('hforc','var')
            string_legend = {"RPM Forcing","Shelf","Fjord","FjordRPM_{best}"};
        else
            string_legend = {"Shelf","Fjord","FjordRPM_{best}"};
        end

        if length(tgt_days)==1
            % string_legend{end+1} = sprintf("mean_{%s}",name_days{i_tgt_day});
            string_legend{end+1} = sprintf("FjordRPM_{mean}");
        else
            for i_day=1:length(tgt_days)
                % string_legend{end+1} = sprintf("mean_{%s}",name_days{i_day});
                string_legend{end+1} = sprintf("FjordRPM_{mean}");
            end
        end
        % hl1 = legend([hs, hf, hbest, hbest2, hb, hm],string_legend,'fontsize',fsize-4,'Location','best'); 
        if exist('hforc','var')
            leg_handles = [hforc, hs, hf, hbest, hb];
        else
            leg_handles = [hs, hf, hbest, hb];
        end
        hl1 = legend(leg_handles,string_legend,'fontsize',fsize,'Location','best'); 
        % title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        hl1.NumColumns=1;
        % old_pos = hl1.Position;
        % hl1.Position = old_pos+[0.2 0 0 0];
    end
    if i_iter==length(which_fjords)
        xlabel('Salinity','fontsize',fsize+2);
    end
end
ylabel(ht,'Depth (m)','fontsize',fsize+2);

end