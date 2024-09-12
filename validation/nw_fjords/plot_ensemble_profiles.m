function [hf_profiles,hf_series] = plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days,name_days,i_tgt_day,mitgcm,plt_salt,plt_series,verbose,which_fjords,plt_rmse)

n_fjord_runs = length(fjord_model);
w_rmse_t = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)
fsize=16;

if nargin < 10 || isempty(mitgcm),    plt_mitgcm = 0; else, plt_mitgcm=1; end
if nargin < 11 || isempty(plt_salt),  plt_salt   = 0; end
if nargin < 12 || isempty(plt_series),plt_series = 0; end
if nargin < 13 || isempty(verbose),   verbose    = 0; end


if exist('rmse_table',"var"),       clear res_obs; end
rmse_table(size(fjord_model)) = struct("tf_rpm",[],"sf_rpm",[],"ts_rpm",[],"tf_gcm",[],"sf_gcm",[],"ts_gcm",[]);

lcolor = lines(3+length(tgt_days));
if length(fjord_model) > 1
    fig_width = 1500; 
    fig_height = 400*length(fjord_model)/2;
else
    fig_width = 500; 
    fig_height = 300;
end

hf_profiles = figure('Name','Temperature profiles','Position',[40 40 fig_width fig_height]);
ht_temp = tiledlayout("flow");
% tiledlayout(length(fjord_model)/2,2);

if plt_salt
    hfs_profiles = figure('Name','Salinity profiles','Position',[40 40 fig_width fig_height]);
    ht_salt = tiledlayout("flow");
end
if plt_series
    hf_series   = figure('Name','Temperature evolution','Position',[40 40 fig_width fig_height]);
    % tiledlayout(length(fjord_model)/2,2);
    tiledlayout("flow");
end


for i_fjord=1:n_fjord_runs

    % find run with the smallest RMSE
    if isempty(res_box(i_fjord).rmse_tf)
        continue
    elseif isempty(i_tgt_day)
        rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf./mean(res_box(i_fjord).tf,'omitnan') + (1-w_rmse_t).*res_box(i_fjord).rmse_sf./mean(res_box(i_fjord).sf,'omitnan');

        [rmse_table(i_fjord).tf_rpm,i_min_rmse_tf] = min(res_box(i_fjord).rmse_tf,[],'all','omitnan');
        [rmse_table(i_fjord).sf_rpm,i_min_rmse_sf] = min(res_box(i_fjord).rmse_sf,[],'all','omitnan');
        [rmse_table(i_fjord).ts_rpm,i_min_rmse]    = min(rmse_both,[],'all','omitnan');

        [irun_best_tf,id_best_tf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_tf);
        [irun_best_sf,id_best_sf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_sf);
        [irun_best,id_best] = ind2sub([n_runs,length(tgt_days)],i_min_rmse);
    
        tf_best = res_box(i_fjord).ensemble_tf(:,irun_best_tf,id_best_tf);
        sf_best = res_box(i_fjord).ensemble_sf(:,irun_best_sf,id_best_sf);
    
        tf_best2 = res_box(i_fjord).ensemble_tf(:,irun_best,id_best);
        sf_best2 = res_box(i_fjord).ensemble_sf(:,irun_best,id_best);
    
        inds_best_tf = [irun_best_tf,id_best_tf];
        inds_best_sf = [irun_best_sf,id_best_sf];
        inds_best2   = [irun_best,id_best];
    else
        rmse_both = w_rmse_t.*res_box(i_fjord).rmse_tf./mean(res_box(i_fjord).tf,1,'omitnan') + (1-w_rmse_t).*res_box(i_fjord).rmse_sf./mean(res_box(i_fjord).sf,1,'omitnan');

        [rmse_table(i_fjord).tf_rpm,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
        [rmse_table(i_fjord).sf_rpm,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
        [rmse_table(i_fjord).ts_rpm,inds_best2]   = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
    
        tf_best = res_box(i_fjord).ensemble_tf(:,inds_best_tf,i_tgt_day);
        sf_best = res_box(i_fjord).ensemble_sf(:,inds_best_sf,i_tgt_day);
    
        tf_best2 = res_box(i_fjord).ensemble_tf(:,inds_best2,i_tgt_day);
        sf_best2 = res_box(i_fjord).ensemble_sf(:,inds_best2,i_tgt_day);
    end

    % compute RMSE for equivalent MITgcm runs
    if plt_mitgcm
        for i_gcm=1:length(mitgcm)
            if strcmp(mitgcm(i_gcm).id,res_box(i_fjord).id{1})
                tprofile_gcm = interp1(mitgcm(i_gcm).z,mitgcm(i_gcm).Tprofile,res_obs(i_fjord).zf,'linear','extrap');
                sprofile_gcm = interp1(mitgcm(i_gcm).z,mitgcm(i_gcm).Sprofile,res_obs(i_fjord).zf,'linear','extrap');
                rmse_table(i_fjord).tf_gcm = rmse(tprofile_gcm,res_obs(i_fjord).tf,'omitnan');%./mean(res_obs(i_fjord).tf,'omitnan');
                rmse_table(i_fjord).sf_gcm = rmse(sprofile_gcm,res_obs(i_fjord).sf,'omitnan');%./mean(res_obs(i_fjord).sf,'omitnan');
                rmse_table(i_fjord).ts_gcm = w_rmse_t.*rmse_table(i_fjord).tf_gcm./mean(res_obs(i_fjord).tf,'omitnan') + (1-w_rmse_t).*rmse_table(i_fjord).sf_gcm/mean(res_obs(i_fjord).sf,'omitnan');
            end
        end
    end

    % adds Shelf-fjord RMSE to the table
    rmse_table(i_fjord).rmse_ts = res_box(i_fjord).rmse_ts;
    rmse_table(i_fjord).rmse_ss = res_box(i_fjord).rmse_ss;

    if verbose
        fprintf("Best parametres for (%s) %s: \n",res_box(i_fjord).id,res_box(i_fjord).name)
        fprintf("Param\t Temperature\t Salinity\t Both\n")
        for i_param=1:length(param_names)
            fprintf("%s \t %.1e\t",param_names{i_param},ensemble(i_fjord,inds_best_tf).p.(param_names{i_param}))
            fprintf("%.1e\t\t",ensemble(i_fjord,inds_best_sf).p.(param_names{i_param}))
            fprintf("%.1e\n",ensemble(i_fjord,inds_best2).p.(param_names{i_param}))
        end
        % fprintf("day \t %s\t\t%s\t\t%s\n",name_days{inds_best_tf(end)},name_days{inds_best_sf(end)},name_days{inds_best2(end)})
        disp("============================")
    end
    %% Plotting temperature
    if nargin > 13
        plot_fjord=1;
        for i_tgt_fjords=1:length(which_fjords)
            if strcmp(which_fjords{i_tgt_fjords},res_box(i_fjord).id{1}) == 1
                plot_fjord=1;
                break
            else
                plot_fjord=0;
            end
        end
        if ~plot_fjord, continue; end
    end
    figure(hf_profiles)
    nexttile; hold on; box on; grid on
    text(0.02,1.02,sprintf("%s) %s (%.0f km)",res_box(i_fjord).id,res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
    text(0.98,0.02,sprintf("n=%.1f %%",res_box(i_fjord).n),'Units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',fsize-2)

    % Observed shelf and fjord profiles
    hs = plot(res_obs(i_fjord).ts,-res_obs(i_fjord).zs,'linewidth',2.5,'color',lcolor(1,:));
    hf = plot(res_obs(i_fjord).tf,-res_obs(i_fjord).zf,'linewidth',2.5,'color',lcolor(2,:));
    hbest = plot(tf_best,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3,:));
    hbest2 = plot(tf_best2,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3,:),'LineStyle','--');
    hb=[];
    % 
    % % RPM profile(s)
    % forcing
    if isfield(res_box(i_fjord),'Tforc')
        hforc = plot(res_box(i_fjord).Tforc,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(1,:),'LineStyle','--');
    end
    % results
    for i_day=1:length(i_tgt_day)
        % res_box(i_fjord).zf = res_box(i_fjord).zf';

        y2 = [-res_box(i_fjord).zf; flip(-res_box(i_fjord).zf)];
        inBetween = [res_box(i_fjord).tfmin(:,i_tgt_day(i_day)); flip(res_box(i_fjord).tfmax(:,i_tgt_day(i_day)))];
        try
            hfjd = patch(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
            hfjd = [];
        end
        try
            inBetween = [res_box(i_fjord).tf1sl(:,i_tgt_day(i_day)); flip(res_box(i_fjord).tf1su(:,i_tgt_day(i_day)))];
            fill(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
        catch
            warning('Ensemble dimensions do not support shading')
            hfjd = [];
        end
        hb_d = plot(res_box(i_fjord).tf(:,i_tgt_day(i_day)),-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(3+i_day,:));
        hb = [hb hb_d];
    end

    % GCM profile (if exists for that fjord)
    if plt_mitgcm
        for i_gcm=1:length(mitgcm)
            if strcmp(mitgcm(i_gcm).id,res_box(i_fjord).id{1})
                hm = plot(mitgcm(i_gcm).Tprofile,-mitgcm(i_gcm).z,'-k','linewidth',2.5);
            end
        end
    else
        hm = [];
    end

    scatter(0,-fjord_model(i_fjord).p.Hgl,40,'v','filled','MarkerFaceColor','black')
    plot([0 0],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)

    set(gca,'fontsize',fsize)
    % xlim([-2.5 9])
    xlim([-2 7.5])
    ylim([-fjord_model(i_fjord).p.H 0])
    if i_fjord==1 
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
        string_legend{end+1} = 'MITgcm';
        % hl1 = legend([hs, hf, hbest, hbest2, hb, hm],string_legend,'fontsize',fsize,'Location','east'); 
        if exist('hforc','var')
            leg_handles = [hforc, hs, hf, hbest, hb, hm];
        else
            leg_handles = [hs, hf, hbest, hb, hm];
        end
        hl1 = legend(leg_handles,string_legend,'fontsize',fsize,'Location','best'); 
        % title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        hl1.NumColumns=1;
        % old_pos = hl1.Position;
        % hl1.Position = old_pos+[0.2 0 0 0];
    end

    %% Plotting salinity
    if plt_salt
        figure(hfs_profiles)
        nexttile; hold on; box on; grid on
        text(0.02,1.02,sprintf("(%s) %s (%.0f km long)",res_box(i_fjord).id{1},res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
        % text(0.02,0.02,sprintf("n=%.1f %%",res_box(i_fjord).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)
    
        % Observed shelf and fjord profiles
        plot(res_obs(i_fjord).ss,-res_obs(i_fjord).zs,'linewidth',1.5,'color',lcolor(1,:));
        plot(res_obs(i_fjord).sf,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(2,:));
        plot(sf_best,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
        plot(sf_best2,-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(3,:));
        
    
        % RPM profile(s)
        % forcing
        if isfield(res_box(i_fjord),'Sforc')
            plot(res_box(i_fjord).Sforc,-res_box(i_fjord).zf,'linewidth',2.5,'color',lcolor(1,:),'LineStyle','--');
        end
        for i_day=1:length(i_tgt_day)
            y2 = [-res_box(i_fjord).zf'; flip(-res_box(i_fjord).zf')];
            try
                inBetween = [res_box(i_fjord).sfmin(:,i_tgt_day(i_day)); flip(res_box(i_fjord).sfmax(:,i_tgt_day(i_day)))];
                fill(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
                inBetween = [res_box(i_fjord).sf1sl(:,i_tgt_day(i_day)); flip(res_box(i_fjord).sf1su(:,i_tgt_day(i_day)))];
                fill(inBetween, y2, lcolor(3+i_day,:),'edgecolor','none','facealpha',0.1);
            catch
                warning('Ensemble dimensions do not support shading')
                hfjd = [];
            end
            plot(res_box(i_fjord).sf(:,i_tgt_day(i_day)),-res_obs(i_fjord).zf,'linewidth',1.5,'color',lcolor(3+i_day,:));
        end
    
        % GCM profile (if exists for that fjord)
        if plt_mitgcm
            for i_gcm=1:length(mitgcm)
                if strcmp(mitgcm(i_gcm).id,res_box(i_fjord).id{1})
                    plot(mitgcm(i_gcm).Sprofile,-mitgcm(i_gcm).z,'-k','linewidth',1.5);
                end
            end
        end
        scatter(33,-fjord_model(i_fjord).p.Hgl,40,'v','filled','MarkerFaceColor','black')
        plot([33 33],[-fjord_model(i_fjord).p.H -fjord_model(i_fjord).p.Hsill],'-k','linewidth',2)
    
        set(gca,'fontsize',fsize)
        xlim([30.5 35])
        ylim([-fjord_model(i_fjord).p.H 0])
    end

    %% Plotting series
    if plt_series
        figure(hf_series)
        nexttile; hold on; box on; grid on
        text(0.02,0.99,sprintf("(%s) %s (%.0f km; n=%d)",res_box(i_fjord).id{1},res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3,res_box(i_fjord).n),'units','normalized','VerticalAlignment','top','fontsize',fsize)
        % text(0.02,0.02,sprintf("n=%d",res_box(i).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)
        if length(res_box(i_fjord).t) == length(res_box(i_fjord).Tupper)
            
            hu = plot(res_box(i_fjord).t,res_box(i_fjord).Tupper,'linewidth',1.5,'color',lcolor(1,:));
            hi = plot(res_box(i_fjord).t,res_box(i_fjord).Tinter,'linewidth',1.5,'color',lcolor(2,:));
            hl = plot(res_box(i_fjord).t,res_box(i_fjord).Tlower,'linewidth',1.5,'color',lcolor(3,:));
           
            for i_day=1:length(tgt_days)
                xline(tgt_days(i_day),'linestyle','--','color',lcolor(3+i_day,:),'linewidth',2)
            end
        end
        ylim([-2 5.5])
        if mod(i_fjord,2) > 0, ylabel('Temperature (^oC)'); end
        if i_fjord>n_fjord_runs-2 
            xlabel('Model time (days)'); 
        else
            set(gca,'xticklabels',[]);
       end
        set(gca,'fontsize',fsize)
        if i_fjord==1%n_fjords
            hl2 = legend([hu,hi,hl],{"0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','best'); 
            % title(hl2,'Time series')
            hl2.NumColumns=1;
        end
    end
end

if plt_salt
    xlabel(ht_salt,'Salinity','fontsize',fsize+2);  
    ylabel(ht_salt,'Depth (m)','fontsize',fsize+2);
end
xlabel(ht_temp,'Temperature (^oC)','fontsize',fsize+2);  
ylabel(ht_temp,'Depth (m)','fontsize',fsize+2);

%% RMSE table
fprintf('\n')
fprintf("RMSE table\n")
fprintf('Fjord | RPM_t  | RPM_s  | RPM_ts | GCM_t  | GCM_s  | GCM_ts |\n')
fprintf('-------------------------------------------------------------------\n')
for i_fjord=1:n_fjord_runs
    fprintf("%s     | %.4f | %.4f | %.4f ",res_box(i_fjord).id,rmse_table(i_fjord).tf_rpm,rmse_table(i_fjord).sf_rpm,rmse_table(i_fjord).ts_rpm)
    if ~isempty(rmse_table(i_fjord).tf_gcm)
        fprintf('| %.4f | %.4f | %.4f ',rmse_table(i_fjord).tf_gcm,rmse_table(i_fjord).sf_gcm,rmse_table(i_fjord).ts_gcm)
    end
    fprintf('|\n')
end
fprintf('-------------------------------------------------------------------\n')

%% RMSE scatter plot

if nargin > 14 && plt_rmse
    mk_sz=100;
    fjord_names = cell(size(fjord_model));
    figure; hold on; box on;
    % h_tf = scatter(-1,0,mk_sz,'^','MarkerEdgeColor',[0 0 0]);
    % h_sf = scatter(-1,0,mk_sz,'v','MarkerEdgeColor',[0 0 0]);
    % h_ts = scatter(-1,0,mk_sz,'square','MarkerEdgeColor',[0 0 0]);
    h_rpm = scatter(-1,0,mk_sz,'filled','o','MarkerFaceColor',lcolor(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.95); 
    h_shf = scatter(-1,0,mk_sz,'filled','o','MarkerFaceColor',lcolor(1,:),'MarkerEdgeColor','none');
    if plt_mitgcm
        h_gcm = scatter(-1,0,mk_sz,'filled','o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
    end
    
    for i_fjord=1:n_fjord_runs
        scatter(i_fjord,rmse_table(i_fjord).tf_rpm,mk_sz,'o','MarkerFaceColor',lcolor(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
        % scatter(i_fjord,rmse_table(i_fjord).sf_rpm,mk_sz,'v','MarkerFaceColor',lcolor(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.95);
        % scatter(i_fjord,rmse_table(i_fjord).ts_rpm,mk_sz,'square','MarkerFaceColor',lcolor(3,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.95);
        scatter(i_fjord,rmse_table(i_fjord).rmse_ts,mk_sz,'o','MarkerFaceColor',lcolor(1,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
        if ~isempty(rmse_table(i_fjord).tf_gcm) && plt_mitgcm
            scatter(i_fjord,rmse_table(i_fjord).tf_gcm,mk_sz,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
            % scatter(i_fjord,rmse_table(i_fjord).sf_gcm,mk_sz,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
            % scatter(i_fjord,rmse_table(i_fjord).ts_gcm,mk_sz,'square','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
        end
        fjord_names{i_fjord} = res_box(i_fjord).id;
    end
    
    % ylim([-6 6])
    xlim([0 n_fjord_runs+1])
    hline(0,'color',[0.5 0.5 0.5],'linestyle','--')
    set(gca,'Xtick',0:1:n_fjord_runs+1)
    xlabels = get(gca,'XTickLabels');
    xlabels(2:end-1) = fjord_names;
    xlabels{1} = ' '; xlabels{end} = ' ';
    set(gca,'XtickLabels',xlabels,'FontSize',14);
    ylabel('RMSE rel. fjord cast (^oC)');
    xlabel('Fjord')
    if plt_mitgcm
        leg_handles = [h_rpm,h_gcm,h_shf];
        leg_labels  = {'FjordRPM','MITgcm','shelf'};
    else
        leg_handles = [h_rpm,h_shf];
        leg_labels  = {'FjordRPM','shelf'};
    end
    legend(leg_handles,leg_labels,'Location','best');
end
end