function [hf_t,hf_e] = plot_sensitivity_ensemble(X,ensemble,res_box,param_names,which_fjords)


if nargin > 4
    n_fjords = length(which_fjords);
else
    n_fjords = size(ensemble,1);
end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end
n_params   = size(X,2);
letters=lower(char(65:65+n_params*n_fjords));

%% Finding the low/mid/high ranges for the different parameters
key_param_bnd = {'low ','mid ','high '};
ls_bnds = {':','-.','-'};
n_cols = 2;

% profiles: low, mid, and high
n_std      = 1; % how many standard deviations away from the mean we want our central interval to span
param_bnds_profile = NaN([4,n_params]);
for i_param=1:n_params
    std_param = std(X(:,i_param));
    avg_param = mean(X(:,i_param));
    param_bnds_profile(:,i_param) = [min(X(:,i_param)),avg_param-n_std.*std_param,avg_param+n_std.*std_param,max(X(:,i_param))]';
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
lcolor = lines(length(param_names));
lbl_fjords = cell([1,length(param_names)]);


%% Plotting temperature
hf_t = figure('Name','Temperature sensitivity','Position',[40 40 1200 900]);
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
    i_plt_fjord=i_plt_fjord+1;
    i_plt = 1+(i_plt_fjord-1)*2*(n_cols*4);
    i_plt_sub = i_plt+(n_cols*4);
    for i_param=1:length(param_names)
        ha_main = nexttile(i_plt,[2 n_cols]); hold on; box on
        text(0.98,1.01,['(',letters(i_panel),')'],'HorizontalAlignment','right','VerticalAlignment','bottom','Units','normalized','fontsize',12)
        i_panel=i_panel+1;
        base_gl_and_sill_t = 0;
        for i_bnd=1:length(key_param_bnd)
            tf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            znb_ensemble = NaN([1,size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,2);
                    znb_ensemble(i_run) = ensemble(i_fjord,i_run).s.znb(2);
                end
                Hsill = ensemble(i_fjord,i_run).p.Hsill;
                Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                H     = ensemble(i_fjord,i_run).p.H;
                has_sill = ensemble(i_fjord,i_run).p.sill;
            end

            depths = -res_box(i_fjord).zf;
            tfmean = mean(tf_ensemble,2,'omitnan');
    
            plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            if i_bnd==length(key_param_bnd)
                hp = plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            end
            % add depth-range of plume neutral buoyancy
            plot([base_gl_and_sill_t+0.1*i_bnd,base_gl_and_sill_t+0.1*i_bnd],[mean(znb_ensemble,'omitnan')-2*std(znb_ensemble,'omitnan'),...
                                                mean(znb_ensemble,'omitnan')+2*std(znb_ensemble,'omitnan')],...
                                                'linewidth',1.7,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd})
        end
        % add depictions of GL and sill depths
        scatter(base_gl_and_sill_t,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if has_sill
            plot([base_gl_and_sill_t base_gl_and_sill_t],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
        end

        set(gca,'fontsize',14)
        % xlim([-2 5])
        ylim([-H 0])
        if i_param==1
            text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
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
            scatter(i_bin*10,tfplot,50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
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
% if no_legend==1
    legend(ha_main,handle_fjords,param_names,'fontsize',10,'Location','southeast');
    % no_legend = 0;
% end
%% Plotting shelf export
hf_e = figure('Name','Shelf export sensitivity','Position',[40 40 1200 900]);
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
    i_plt = 1+(i_plt_fjord-1)*2*(n_cols*4);
    i_plt_sub = i_plt+(n_cols*4);
    for i_param=1:length(param_names)    
        ha_main = nexttile(i_plt,[2 n_cols]); hold on; box on
        text(0.98,1.01,['(',letters(i_panel),')'],'HorizontalAlignment','right','VerticalAlignment','bottom','Units','normalized','fontsize',12)
        i_panel=i_panel+1;
        for i_bnd=1:length(key_param_bnd)
            
            ef_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_profile(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    Sref = 35.0;
                    ef_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.QVsfinal(:,1).*((Sref-ensemble(i_fjord,i_run).s.Sfinal(:,1))/Sref);
                end
                Hsill = ensemble(i_fjord,i_run).p.Hsill;
                Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                H     = ensemble(i_fjord,i_run).p.H;
                has_sill = ensemble(i_fjord,i_run).p.sill;
            end
            vline(0,'linewidth',0.5,'linestyle',':','color',[0.7 0.7 0.7])

            % take mean and min/max for that subset
            depths = -res_box(i_fjord).zf;
            efmean = mean(ef_ensemble,2,'omitnan');
            
            plot(efmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            if i_bnd==length(key_param_bnd)
                hp = plot(efmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            end
        end
        
        % add depictions of GL and sill depths
        base_gl_and_sill_t = 0;
        scatter(base_gl_and_sill_t,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if has_sill
            plot([base_gl_and_sill_t base_gl_and_sill_t],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
        end
        
        set(gca,'fontsize',14)
        ylim([-H 0])
        if no_legend==1
            handle_fjords = [handle_fjords hp];
            lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
        end
        if i_param==1
            text(gca,0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
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
            % find profiles that fit into that bin for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds_prctile(i_fjord,i_run).(param_names{i_param}) == i_bin
                    fz_ensemble(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
                end
                fz_full_ens(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
            end
            fzplot = mean(fz_ensemble,[1,2],'omitnan'); % take mean for that subset
            scatter(i_bin*10,round(fzplot,2),50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        end
        set(gca,'YAxisLocation','right','XAxisLocation','top','fontsize',8)
        xlabel([param_names{i_param},' percentile'],'fontsize',10)
        % xlabel('percentile','fontsize',10)
        ylabel('z export (m)','fontsize',10)
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
legend(ha_main,handle_fjords,param_names,'fontsize',10,'Location','southeast');
% no_legend = 0;


figure(hf_t)
xlabel(ht_t,'Temperature (^oC)','fontsize',14);
ylabel(ht_t,'Depth (m)','fontsize',14);
ht_t.TileSpacing='compact';
ht_t.Padding='compact';

figure(hf_e)
xlabel(ht_e,'Shelf-fjord freshwater flux (m^3s^{-1})','fontsize',14);
ylabel(ht_e,'Depth (m)','fontsize',14);
ht_e.TileSpacing='compact';
ht_e.Padding='compact';

end