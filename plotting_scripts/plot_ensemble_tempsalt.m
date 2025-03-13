function hf_ts = plot_ensemble_tempsalt(fjord_model_yr,ensemble_yr,res_box_yr,res_obs_yr,n_runs,tgt_days,i_tgt_day,which_fjords,which_yr)

if isempty(tgt_days)
    n_days = 1;
else
    n_days = length(tgt_days);
end

fjord_model = fjord_model_yr{which_yr};
res_box = res_box_yr{which_yr};
res_obs = res_obs_yr{which_yr};

n_fjord_runs = length(fjord_model);
w_rmse_t     = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)

fsize   = 14;
lcolor  = lines(3+length(tgt_days));
letters = lower(char(65:65+length(which_fjords)*2+10));

fig_width = 1200; %900; 
fig_height = 300*length(which_fjords);
% rmse_table(size(fjord_model)) = struct("tf_rpm",[],"sf_rpm",[],"ts_rpm",[],"tf_gcm",[],"sf_gcm",[],"ts_gcm",[]);
% rmse_table = cell(size(fjord_model));

hf_ts    = figure('Name','Temperature and salinity profiles','Position',[40 40 fig_width fig_height]);
ht       = tiledlayout(length(which_fjords),3,'TileSpacing','loose','padding','compact');
i_iter   = 0;
i_letter =1;
i_row=1;
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
    nexttile(i_row); hold on; box on; grid on
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
    % scatter(0,-fjord_model(i_fjord).p.Hgl,60,'v','filled','MarkerFaceColor','black')
    plot([-5 5],[-fjord_model(i_fjord).p.Hgl -fjord_model(i_fjord).p.Hgl],'linestyle','--','color',[0.3 0.3 0.3])
    if fjord_model(i_fjord).p.sill
        % plot([0 0],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)
        plot([-5 5],[-fjord_model(i_fjord).p.Hsill -fjord_model(i_fjord).p.Hsill],'linestyle','-','color',[0.3 0.3 0.3])
    end

    set(gca,'fontsize',fsize)
    % xlim([-2.5 9])
    % xlim([-4 7.5])
    % ylim([-fjord_model(i_fjord).p.H 0])
    xlim ([-2 4])  % as per Tom's suggestion
    ylim([-600 0]) % as per Tom's suggestion
    if i_iter==length(which_fjords)
        xlabel('Temperature (^oC)','fontsize',fsize+2);  
    end
    
    %% Plotting salinity    
    nexttile(i_row+1); hold on; box on; grid on
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
    % scatter(33,-fjord_model(i_fjord).p.Hgl,60,'v','filled','MarkerFaceColor','black')
    plot([30 36],[-fjord_model(i_fjord).p.Hgl -fjord_model(i_fjord).p.Hgl],'linestyle','--','color',[0.3 0.3 0.3])
    if fjord_model(i_fjord).p.sill
        % plot([33 33],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)
        plot([30 36],[-fjord_model(i_fjord).p.Hsill -fjord_model(i_fjord).p.Hsill],'linestyle','-','color',[0.3 0.3 0.3])
    end

    set(gca,'fontsize',fsize)
    xlim([30.5 35])
    % ylim([-fjord_model(i_fjord).p.H 0])
    ylim([-600 0]) % as per Tom's suggestion

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
    i_row  = i_row+3;
end
ylabel(ht,'Depth (m)','fontsize',fsize+2);

%% Adding the scatter plots

n_years    = length(fjord_model_yr);
first_time = 1;
all_r_temp = {};
all_r_salt = {};

for i_year=n_years:-1:1
    n_fjord_runs = length(fjord_model_yr{i_year});
    res_box      = res_box_yr{i_year};
    res_obs      = res_obs_yr{i_year};
    for i_fjord=1:n_fjord_runs
        
        % Getting the profiles for comparing
        [~,tf_best,sf_best] = get_best_profiles_rmse(res_box,i_fjord,n_runs,w_rmse_t,tgt_days,i_tgt_day);

        % interpolaring observations to box model depths
        tf_obs = interp1(res_obs(i_fjord).zf,res_obs(i_fjord).tf,res_box(i_fjord).zf,'linear');
        sf_obs = interp1(res_obs(i_fjord).zf,res_obs(i_fjord).sf,res_box(i_fjord).zf,'linear');
        ts_obs = res_box(i_fjord).Tforc;
        ss_obs = res_box(i_fjord).Sforc;

        dt_obs = tf_obs  - ts_obs;
        dt_rpm = tf_best - ts_obs;
        ds_obs = sf_obs  - ss_obs;
        ds_rpm = sf_best - ss_obs;
        

        if first_time
            nexttile(3); hold on; box on; grid on
            text(0.02,0.98,sprintf("(%s)",letters(i_letter)),'Units','normalized','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
            i_letter=i_letter+1;
            xlabel('T_{obs}') ; ylabel('T_{rpm}')
            % xlabel('\Delta T_{obs}') ; ylabel('\Delta T_{rpm}')
            plot([-10, 10],[-10, 10],'-r')

            nexttile(6); hold on; box on; grid on
            text(0.02,0.98,sprintf("(%s)",letters(i_letter)),'Units','normalized','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
            i_letter=i_letter+1;
            xlabel('S_{obs}') ; ylabel('S_{rpm}')
            % xlabel('\Delta S_{obs}') ; ylabel('\Delta S_{rpm}')
            % plot([10, 36],[10, 36],'-r')
            plot([-10, 36],[-10, 36],'-r')
            first_time = 0;

        end

        nexttile(3); hold on; box on; grid on
        scatter(tf_obs,tf_best,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);
        % scatter(dt_obs,dt_rpm,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);
        
        nexttile(6); hold on; box on; grid on
        scatter(sf_obs,sf_best,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);
        % scatter(ds_obs,ds_rpm,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);

        all_r_temp{end+1} = corr(tf_obs,tf_best,'type','pearson','rows','complete');
        all_r_salt{end+1} = corr(sf_obs,sf_best,'type','pearson','rows','complete');
        % all_r_temp{end+1} = corr(dt_obs,dt_rpm,'type','pearson','rows','complete');
        % all_r_salt{end+1} = corr(ds_obs,ds_rpm,'type','pearson','rows','complete');
    end % fjords
end % years

nexttile(3);
xlim([-2, 6]); ylim([-2, 6])
% xlim([-10, 5]); ylim([-10, 5])
hcb = colorbar;
ylabel(hcb,'Depth (m)','FontSize',fsize)
set(gca,'fontsize',fsize)

nexttile(6);
xlim([25, 36]); ylim([25, 36])
% xlim([-2, 2]); ylim([-2, 2])
set(gca,'fontsize',fsize)
colormap(flip(cmocean('deep')))

nexttile(9)
text(0.02,0.98,sprintf("(%s)",letters(i_letter)),'Units','normalized','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
bin_edges = -0.1:0.05:1.0;
hold on; box on
histogram(cell2mat(all_r_temp),bin_edges,'Normalization','count','FaceAlpha',0.5);
histogram(cell2mat(all_r_salt),bin_edges,'Normalization','count','FaceAlpha',0.5);
legend('Temperature','Salinity','Location','north')
xlabel('FjordRPM_{best} fit to obs.'); ylabel('Count')
set(gca,'fontsize',fsize)


end