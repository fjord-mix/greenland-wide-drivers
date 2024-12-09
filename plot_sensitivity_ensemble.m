function [hf_t,hf_e] = plot_sensitivity_ensemble(X,ensemble,res_box,res_obs,param_names,which_fjords,plt_fw)

if nargin < 7, plt_fw = 0; hf_e=[]; end

if nargin > 4
    n_fjords = length(which_fjords);
else
    n_fjords = size(ensemble,1);
end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end
n_params   = size(X,2);
letters=lower(char(65:65+n_params*n_fjords));

lcolor = lines(4);
fsize = 14;
fig_width = 1200;
fig_height = 250*length(which_fjords);
%% Applying the RMSE filter - speculative, for now trying out to see if it works

% rmse_threshold = 0.5;
% for i_fjord=1:size(ensemble,1)
%     for i_run=1:size(ensemble,2)
%         if ~isempty(ensemble(i_fjord,i_run).p)
%             if res_box(i_fjord).rmse_tf(i_run,2) > rmse_threshold
%                 ensemble(i_fjord,i_run).s = []; % we clear the outputs structure, which will make L112 skip this ensemble member
%             end
%         end
%     end
% end

%% Finding the low/mid/high ranges for the different parameters
key_param_bnd = {'low ','mid ','high '};
ls_bnds = {':','-.','-'};
n_cols = 2;

% profiles: low, mid, and high
param_bnds_profile = NaN([4,n_params]);
for i_param=1:n_params
    param_bnds_profile(:,i_param) = prctile(X(:,i_param),[0,33,66,100]);
    % param_bnds_profile(:,i_param) = quantile(X(:,i_param),4); % no qualitative changes compared with using percentiles, but the latter makes the differences clearer
end

% scatter: 10 bins of percentiles
bins       = 0:10:100;
param_bnds_prctile = NaN([length(bins),n_params]);
for i_param=1:n_params
    param_bnds_prctile(:,i_param) = prctile(X(:,i_param),bins);
end

% loops over all ensemble entries
struct_fields = param_names;
struct_fields{2,1} = cell(size(ensemble));
mask_bnds_profile = struct(struct_fields{:});
mask_bnds_prctile = struct(struct_fields{:});
for i_fjord=1:size(ensemble,1)
    for i_run=1:size(ensemble,2)
        if ~isempty(ensemble(i_fjord,i_run).p)
            for i_param=1:n_params
                % for each parameter, assigns a mask of {1|2|3} if the value is {low|mid|high} according to the defined param_bnds
                if ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds_profile(1,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds_profile(2,i_param)
                    mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) = 1;
                elseif ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds_profile(2,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds_profile(3,i_param)
                    mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) = 2;
                elseif ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds_profile(3,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds_profile(4,i_param)
                    mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) = 3;
                end
                % for each parameter, assigns a mask of 1:10 if the value is is in the 1:100th percentile
                for i_bin=1:length(bins)-1
                    if ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds_prctile(i_bin,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds_prctile(i_bin+1,i_param)
                        mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) = i_bin;
                    end
                end
            end
        end
    end
end

%% Preparing the figures
handle_fjords = [];
lbl_fjords = cell([1,length(param_names)]);


%% Plotting temperature
hf_t = figure('Name','Temperature sensitivity','Position',[40 40 fig_width fig_height]);
ht_t = tiledlayout(2*n_fjords,n_cols*length(param_names));
no_legend = 1;
i_plt_fjord=0;
i_panel=1;
for i_fjord=1:size(ensemble,1)
    if nargin > 4
        for i_tgt_fjords=1:n_fjords
            if strcmp(which_fjords{i_tgt_fjords},res_box(i_fjord).id) == 1
                plot_fjord=1;
                break
            else
                plot_fjord=0;
            end
        end
        if ~plot_fjord, continue; end
    end
    % find best run with the smallest RMSE
    if isempty(res_box(i_fjord).rmse_tf)
        tf_best = NaN(size(res_box(i_fjord).zf));
    else 
        [~,tf_best,~] = get_best_profiles_rmse(res_box,i_fjord,size(ensemble,2));
    end
    


    i_plt_fjord=i_plt_fjord+1;
    i_plt = 1+(i_plt_fjord-1)*2*(n_cols*n_params);
    i_plt_sub = i_plt+(n_cols*n_params);
    for i_param=1:length(param_names)
        ha_main = nexttile(i_plt,[2 n_cols]); hold on; box on
        text(0.99,1.01,['(',letters(i_panel),') ',param_names{i_param}],'HorizontalAlignment','right','VerticalAlignment','bottom','Units','normalized','fontsize',fsize)
        i_panel=i_panel+1;

        base_gl_and_sill_t = 1;
        base_gl_and_sill_p = -2;

        for i_bnd=1:length(key_param_bnd)
            tf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            znb_ensemble = NaN([1,size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,2);

                    % same time as T profile
                    % znb_ensemble(i_run) = ensemble(i_fjord,i_run).s.znb(2);
                    
                    % mean of 30 days around peak discharge of last year
                    % (mean of 30 days after peak discharge of last year is not that different from using the September conditions)
                    % [~,i_peak] = max(ensemble(i_fjord,i_run).s.Qsg(end-365:end));
                    % znb_ensemble(i_run) = mean(ensemble(i_fjord,i_run).s.znb_t(end-365+i_peak-15:end-365+i_peak+15)); 

                    % mean over melt season of the last year
                    i_discharge = ensemble(i_fjord,i_run).s.Qsg > 0 & ensemble(i_fjord,i_run).s.t > ensemble(i_fjord,i_run).s.t(end)-365;
                    znb_ensemble(i_run) = mean(ensemble(i_fjord,i_run).s.znb_t(i_discharge),'omitnan'); 

                end
                Hsill = ensemble(i_fjord,i_run).p.Hsill;
                Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                H     = ensemble(i_fjord,i_run).p.H;
                has_sill = ensemble(i_fjord,i_run).p.sill;
            end


            depths = -res_box(i_fjord).zf;
            tfmean = mean(tf_ensemble,2,'omitnan');
            if i_bnd==1
                % Plotting shelf forcing observation
                plot(res_box(i_fjord).Tforc,depths,'linewidth',1.0,'color',lcolor(1,:),'linestyle','-');

                % Plotting fjord observation
                % plot(res_obs(i_fjord).tf,-res_obs(i_fjord).zf,'linewidth',1.0,'color',lcolor(2,:),'linestyle','-');
    
                % Plotting fjord ensemble best run
                % plot(tf_best,depths,'linewidth',1.0,'color',lcolor(3,:),'linestyle','-');
            end

            % Plot sensitivity profile
            plot(tfmean,depths,'linewidth',2.0,'color',lcolor(4,:),'linestyle',ls_bnds{i_bnd});
            if i_bnd==length(key_param_bnd)
                hp = plot(tfmean,depths,'linewidth',2.0,'color',lcolor(4,:),'linestyle',ls_bnds{i_bnd});
            end

            % add depth-range of plume neutral buoyancy
            % scatter(base_gl_and_sill_t+0.1*i_bnd,mean(znb_ensemble,'omitnan'),60,lcolor(i_param,:),'x')
            % plot([base_gl_and_sill_t+0.1*i_bnd,base_gl_and_sill_t+0.1*i_bnd],[mean(znb_ensemble,'omitnan')-2*std(znb_ensemble,'omitnan'),...
            %                                     mean(znb_ensemble,'omitnan')+2*std(znb_ensemble,'omitnan')],...
            %                                     'linewidth',1.7,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd})
            scatter(base_gl_and_sill_p+0.2*i_bnd,mean(znb_ensemble,'omitnan'),60,lcolor(4,:),'x')
            plot([base_gl_and_sill_p+0.2*i_bnd,base_gl_and_sill_p+0.2*i_bnd],[mean(znb_ensemble,'omitnan')-2*std(znb_ensemble,'omitnan'),...
                                                mean(znb_ensemble,'omitnan')+2*std(znb_ensemble,'omitnan')],...
                                                'linewidth',1.7,'color',lcolor(4,:),'linestyle',ls_bnds{i_bnd})
        end
        % add depictions of GL and sill depths
        scatter(base_gl_and_sill_t,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if has_sill
            plot([base_gl_and_sill_t base_gl_and_sill_t],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
        end

        set(gca,'fontsize',fsize)
        xlim([-3 3])
        ylim([-H 0])
        if i_param==1
            text(0.02,1.01,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
        end
        if no_legend==1
            handle_fjords = [handle_fjords hp];
            lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
        end

        % Adding inset scatter plot
        ax2 = axes(ht_t);
        ax2.Layout.Tile=i_plt_sub;
        ax2.Layout.TileSpan=[1 n_cols-1];
        ax2.Box = 'on';
        hold on;
        % nexttile(i_plt_sub,[1 1]); box on;
        for i_bin=1:length(bins)-1
            tf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            tf_full_ens = NaN(size(tf_ensemble));
            % find profiles that fit into that bin for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) == i_bin
                    tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,2);
                end
                tf_full_ens(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,2);
            end
            tfplot = mean(tf_ensemble,[1,2],'omitnan') - mean(tf_full_ens,[1,2],'omitnan');
            scatter(i_bin*10,tfplot,50,'filled','MarkerFaceColor',lcolor(4,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        end
        set(gca,'YAxisLocation','right','XAxisLocation','top','fontsize',8)
        xlabel([param_names{i_param},' percentile'],'fontsize',10)
        % xlabel('percentile','fontsize',10)
        ylabel('T diff. (^oC)','fontsize',10)
        % xticks = get(gca,'XTickLabel');
        % yticks = get(gca,'YTickLabel');
        % set(gca,'XTickLabel',{'',xticks{2:end-1},''})
        set(gca,'XTick',[25, 50, 75])
        ytickangle(45)
        % set(gca,'YTickLabel',{'',yticks{2:end}})
        i_plt = i_plt+n_cols;
        i_plt_sub=i_plt_sub+n_cols;
    end
end
xlabel(ht_t,'Temperature (^oC)','fontsize',fsize);
ylabel(ht_t,'Depth (m)','fontsize',fsize);
ht_t.TileSpacing='compact';
ht_t.Padding='compact';

% if no_legend==1
    % legend(ha_main,handle_fjords,param_names,'fontsize',10,'Location','southeast');
    % no_legend = 0;
% end
%% Plotting shelf export
if plt_fw
hf_e = figure('Name','Shelf export sensitivity','Position',[40 40 fig_width fig_height]);
ht_e = tiledlayout(2*n_fjords,n_cols*length(param_names));
no_legend=1;
i_plt_fjord=0;
i_panel=1;
for i_fjord=1:size(ensemble,1)
    if nargin > 4
        for i_tgt_fjords=1:n_fjords
            if strcmp(which_fjords{i_tgt_fjords},res_box(i_fjord).id) == 1
                plot_fjord=1;
                break
            else
                plot_fjord=0;
            end
        end
        if ~plot_fjord, continue; end
    end
    i_plt_fjord=i_plt_fjord+1;
    i_plt = 1+(i_plt_fjord-1)*2*(n_cols*n_params);
    i_plt_sub = i_plt+(n_cols*n_params)+1;
    for i_param=1:length(param_names)    
        ha_main = nexttile(i_plt,[2 n_cols]); hold on; box on
        text(0.98,1.01,['(',letters(i_panel),')'],'HorizontalAlignment','right','VerticalAlignment','bottom','Units','normalized','fontsize',fsize)
        i_panel=i_panel+1;
        base_gl_and_sill_ef = 0;
        for i_bnd=1:length(key_param_bnd)
            
            ef_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            znb_ensemble = NaN([1,size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    Sref = 35.0;
                    % ef_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.QVsfinal(:,2).*((Sref-ensemble(i_fjord,i_run).s.Sfinal(:,2))/Sref);
                    ef_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.fw_profile_export;%/trapz(ensemble(i_fjord,i_run).s.fw_profile_export,ensemble(i_fjord,i_run).s.z);
                    znb_ensemble(i_run) = ensemble(i_fjord,i_run).s.znb(1);
                end
                Hsill = ensemble(i_fjord,i_run).p.Hsill;
                Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                H     = ensemble(i_fjord,i_run).p.H;
                has_sill = ensemble(i_fjord,i_run).p.sill;
            end
            vline(0,'linewidth',0.5,'linestyle',':','color',[0.7 0.7 0.7])

            % take mean and min/max for that subset
            depths = -res_box(i_fjord).zf;
            efmean = -mean(ef_ensemble,2,'omitnan');
            
            plot(efmean,depths,'linewidth',1.5,'color',lcolor(4,:),'linestyle',ls_bnds{i_bnd});
            if i_bnd==length(key_param_bnd)
                hp = plot(efmean,depths,'linewidth',1.5,'color',lcolor(4,:),'linestyle',ls_bnds{i_bnd});
            end

            % add plume
            % scatter(base_gl_and_sill_ef+1*i_bnd,mean(znb_ensemble,'omitnan'),60,[0 0 0],'x')
            % plot([base_gl_and_sill_ef+1*i_bnd,base_gl_and_sill_ef+1*i_bnd],[mean(znb_ensemble,'omitnan')-2*std(znb_ensemble,'omitnan'),...
            %                                     mean(znb_ensemble,'omitnan')+2*std(znb_ensemble,'omitnan')],...
            %                                     'linewidth',1.7,'color',[0 0 0],'linestyle',ls_bnds{i_bnd})
        end
        
        % add depictions of GL and sill depths
        scatter(base_gl_and_sill_ef,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if has_sill
            plot([base_gl_and_sill_ef base_gl_and_sill_ef],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
        end
        
        set(gca,'fontsize',fsize)
        ylim([-H 0])
        if no_legend==1
            handle_fjords = [handle_fjords hp];
            lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
        end
        if i_param==1
            text(gca,0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',fsize)
        end

        % Adding inset scatter plot
        ax2 = axes(ht_e);
        ax2.Layout.Tile=i_plt_sub;
        ax2.Layout.TileSpan=[1 n_cols-1];
        ax2.Box = 'on';
        hold on;
        for i_bin=1:length(bins)-1
            fz_ensemble = NaN([1,size(ensemble,2)]);
            fz_full_ens = NaN(size(fz_ensemble));
            ef_ensemble_peaks = NaN(size(ef_ensemble));
            % find profiles that fit into that bin for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) == i_bin
                    % fz_ensemble(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
                    % ef_ensemble_peaks(:,i_run) = ensemble(i_fjord,i_run).s.QVsfinal(:,1).*((Sref-ensemble(i_fjord,i_run).s.Sfinal(:,1))/Sref);
                    ef_ensemble_peaks(:,i_run) = ensemble(i_fjord,i_run).s.fw_profile_export;
                    ef_ensemble_peaks(:,i_run) = ef_ensemble_peaks(:,i_run)/trapz(ef_ensemble_peaks(:,i_run),ensemble(i_fjord,i_run).s.z); % normalising the flux magnitude
                end
                fz_full_ens(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
            end

            % finds peaks in the FW export profile. We want those associated with export (i.e., positive & above the grounding line)
            ef_peaks = -mean(ef_ensemble_peaks,2,'omitnan');
            ef_peaks(depths < -Hgl) = 0;
            ef_peaks(ef_peaks < 0) = 0;
            
            % [ef_peaks_rsmp,depths_rsmp] = resample(ef_peaks,-round(depths),1);
            [val_peaks,z_peaks,w_peaks,p_peaks] = findpeaks(ef_peaks,-round(depths));%,'NPeaks',5);
            area_peak = w_peaks.*val_peaks; % Although this is not properly the area, it worked much better than the numerical integral approach below

            % functions to estimate area undr the curve of each peak
            % sfcn = @(b,x) b(4) + b(1).*x.*(exp(b(2).*x) + exp(b(3).*x));
            % AUC = @(b,x) - b(4).*(x(1) - x(end)) - b(1).*((exp(b(2).*x(1)).*(b(2)*x(1) - 1))./b(2).^2 - (exp(b(2).*x(end)).*(b(2)*x(end) - 1))./b(2).^2) - b(1)*((exp(b(3)*x(1)).*(b(3).*x(1) - 1))./b(3)^2 - (exp(b(3)*x(end))*(b(3)*x(end) - 1))./b(3)^2);
            % for k = 1:numel(z_peaks)
            %     % ixrng = -(z_peaks(k)+20:-1:z_peaks(k)-200);
            %     ixrng = max(1,z_peaks(k)-100) : min(z_peaks(k)+100,numel(depths_rsmp));
            %     drv = depths_rsmp(ixrng) - depths_rsmp(ixrng(1));
            %     B(k,:) = fminsearch(@(b)norm(ef_peaks_rsmp(ixrng) - sfcn(b,drv)), [z_peaks(k)*10 -rand(1,2)*10 0.01]);
            %     sfit{k} = [drv+depths_rsmp(ixrng(1)), sfcn(B(k,:),drv), ixrng(:)];
            %     % area_peak(k) = AUC(B(k,:),drv);
            %     area_peak(k,:) = trapz(drv, sfit{k}(:,2));
            % end

            % we are interested in the widest peak, i.e., the one associated with the most freshwater export
            if isempty(val_peaks)
                    fzplot=0;
            else
                fzplot = -z_peaks(area_peak == max(area_peak)); 
            end

            % fzplot = mean(fz_ensemble,[1,2],'omitnan'); % take mean for that subset
            scatter(i_bin*10,round(fzplot,2),50,'filled','MarkerFaceColor',lcolor(4,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        end
        set(gca,'YAxisLocation','left','XAxisLocation','top','fontsize',8)
        xlabel([param_names{i_param},' percentile'],'fontsize',10)
        % xlabel('percentile','fontsize',10)
        ylabel('z export (m)','fontsize',10)
        % xticks = get(gca,'XTickLabel');
        % yticks = get(gca,'YTickLabel');
        % set(gca,'XTickLabel',{'',xticks{2:end-1},''})
        set(gca,'XTick',[25, 50, 75])
        ytickangle(-45)
        % set(gca,'YTickLabel',{'',yticks{2:end}})
        i_plt = i_plt+n_cols;
        i_plt_sub=i_plt_sub+n_cols;
    end
end
% xlabel(ht_e,'Shelf-fjord freshwater flux (m^3s^{-1})','fontsize',fsize);
xlabel(ht_e,'Mean exported freshwater flux (m^3s^{-1})','fontsize',fsize);
ylabel(ht_e,'Depth (m)','fontsize',fsize);
ht_e.TileSpacing='compact';
ht_e.Padding='compact';
% legend(ha_main,handle_fjords,param_names,'fontsize',10,'Location','southwest');
% no_legend = 0;
end
end