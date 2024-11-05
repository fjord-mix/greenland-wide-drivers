function [hf_t,hf_s,hf_e] = plot_sensitivity_profiles_v3(X,ensemble,res_box,param_names,i_day,plt_salt,plt_exp,which_fjords)

if nargin < 5, i_day=1; end
if nargin < 6, plt_salt=0; end
if nargin < 7, plt_exp=0; end
if nargin > 7
    n_fjords = length(which_fjords);
else
    n_fjords = size(ensemble,1);
end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end
no_legend = 1;

%% Finding the low/mid/high ranges for the different parameters
key_param_bnd = {'low ','mid ','high '};
ls_bnds = {':','-.','-'};

n_std      = 1; % how many standard deviations away from the mean we want our central interval to span
n_params   = size(X,2);
param_bnds = NaN([4,n_params]);
for i_param=1:n_params
    std_param = std(X(:,i_param));
    avg_param = mean(X(:,i_param));
    param_bnds(:,i_param) = [min(X(:,i_param)),avg_param-n_std.*std_param,avg_param+n_std.*std_param,max(X(:,i_param))]';
end

% loops over all ensemble entries
struct_fields = param_names;
struct_fields{2,1} = cell(size(ensemble));
mask_bnds = struct(struct_fields{:});
for i_fjord=1:size(ensemble,1)
    for i_run=1:size(ensemble,2)
        % for each parameter, assigns a mask of {1|2|3} if the value is {low|mid|high} according to the defined param_bnds
        if ~isempty(ensemble(i_fjord,i_run).p)
            for i_param=1:n_params
                if ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds(1,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds(2,i_param)
                    mask_bnds(i_fjord,i_run).(param_names{i_param}) = 1;
                elseif ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds(2,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds(3,i_param)
                    mask_bnds(i_fjord,i_run).(param_names{i_param}) = 2;
                elseif ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds(3,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds(4,i_param)
                    mask_bnds(i_fjord,i_run).(param_names{i_param}) = 3;
                end
            end
        end
    end
end

%% Preparing the figures
handle_fjords = [];
lcolor = lines(length(param_names));
lbl_fjords = cell([1,length(param_names)]);


hf_t = figure('Name','Temperature sensitivity','Position',[40 40 800 600]);
ht_t = tiledlayout(n_fjords,length(param_names));

if plt_salt
    hf_s = figure('Name','Salinity sensitivity','Position',[40 40 800 600]);
    ht_s = tiledlayout(n_fjords,length(param_names));
end
if plt_exp
    hf_e = figure('Name','Shelf export sensitivity','Position',[40 40 800 600]);
    ht_e = tiledlayout(n_fjords,length(param_names));
end
for i_fjord=1:size(ensemble,1)
    if nargin > 7
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
    %% Plotting temperature

    figure(hf_t)
    % plot(res_obs(i_fjord).ts,-res_obs(i_fjord).zs,'color',[0 0 0])
    
    for i_param=1:length(param_names)    
        nexttile; hold on; box on
        base_gl_and_sill_t = 0;
        for i_bnd=1:length(key_param_bnd)
            % tf_ensemble = NaN([length(ensemble(i_fjord,1).s.z),size(ensemble,2)]);
            tf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            znb_ensemble = NaN([1,size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,i_day);
                    znb_ensemble(i_run) = ensemble(i_fjord,i_run).s.znb(i_day);
                end
                Hsill = ensemble(i_fjord,i_run).p.Hsill;
                Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                H     = ensemble(i_fjord,i_run).p.H;
                has_sill = ensemble(i_fjord,i_run).p.sill;
            end
            
            % take mean and min/max for that subset
            % depths = -ensemble(i_fjord,1).s.z;
            depths = -res_box(i_fjord).zf;
            tfmean = mean(tf_ensemble,2,'omitnan');
            % tfmin  = min(tf_ensemble,[],2,'omitnan');
            % tfmax  = max(tf_ensemble,[],2,'omitnan');
    
            % plot
            % y2 = [depths; flip(depths)];
            % inBetween = [tfmin; flip(tfmax)];
            % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor',lcolor(i_param,:),'facealpha',0.1);
            % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor','none','facealpha',0.1);
            plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            if i_bnd==length(key_param_bnd)
                hp = plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
            end
            % add depth of plume neutral buoyancy
            % scatter(base_gl_and_sill_t,mean(znb_ensemble,'omitnan'),40,'s','filled','MarkerFaceColor',[0 0 0])
            plot([base_gl_and_sill_t+0.1*i_bnd,base_gl_and_sill_t+0.1*i_bnd],[mean(znb_ensemble,'omitnan')-2*std(znb_ensemble,'omitnan'),...
                                                mean(znb_ensemble,'omitnan')+2*std(znb_ensemble,'omitnan')],...
                                                'linewidth',1.7,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd})
            % errorbar(base_gl_and_sill_t+0.1*i_bnd,mean(znb_ensemble,'omitnan'),2*std(znb_ensemble,'omitnan'),...
            %                                      '.','linewidth',1.,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd})
        end
        % add depictions of GL and sill depths
        scatter(base_gl_and_sill_t,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if has_sill
            plot([base_gl_and_sill_t base_gl_and_sill_t],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
        end

        set(gca,'fontsize',14)
        % xlim([-2 5])
        ylim([-H 0])
        if no_legend==1
            handle_fjords = [handle_fjords hp];
            lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
        end
        if i_param==1
            % text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
            text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
        end
    end
    if no_legend==1
        legend(gca,handle_fjords,param_names,'fontsize',10,'Location','best');
        % hl{i_bnd}.NumColumns=2;
        no_legend = 0;
    end

    %% Plotting salinity
    if plt_salt
        figure(hf_s)
        
        for i_param=1:length(param_names)    
            nexttile; hold on; box on
            for i_bnd=1:length(key_param_bnd)
                % sf_ensemble = NaN([length(ensemble(i_fjord,1).s.z),size(ensemble,2)]);
                sf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
                % find profiles that fit into that interval for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bnd
                        sf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                    end
                    Hsill = ensemble(i_fjord,i_run).p.Hsill;
                    Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                    H     = ensemble(i_fjord,i_run).p.H;
                    has_sill = ensemble(i_fjord,i_run).p.sill;
                end
                
                % take mean and min/max for that subset
                % depths = -ensemble(i_fjord,1).s.z;
                depths = -res_box(i_fjord).zf;
                sfmean = mean(sf_ensemble,2,'omitnan');
                % tfmin  = min(sf_ensemble,[],2,'omitnan');
                % tfmax  = max(sf_ensemble,[],2,'omitnan');
        
                % plot
                % y2 = [depths; flip(depths)];
                % inBetween = [sfmin; flip(tfmax)];
                % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor',lcolor(i_param,:),'facealpha',0.1);
                % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor','none','facealpha',0.1);
                plot(sfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
                if i_bnd==length(key_param_bnd)
                    hp = plot(sfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:),'linestyle',ls_bnds{i_bnd});
                end
            end
            % add depictions of GL and sill depths
            base_gl_and_sill_t = 0;
            scatter(base_gl_and_sill_t,-Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
            if has_sill
                plot([base_gl_and_sill_t base_gl_and_sill_t],[-H -Hsill],'-','linewidth',2,'color',[0 0 0])
            end
            set(gca,'fontsize',14)
            % xlim([-2 5])
            ylim([-H 0])
            if no_legend==1
                handle_fjords = [handle_fjords hp];
                lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
            end
            if i_param==1
                % text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
                text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
            end
        end
        if no_legend==1
            legend(gca,handle_fjords,param_names,'fontsize',10,'Location','best');
            % hl{i_bnd}.NumColumns=2;
            no_legend = 0;
        end
    end

    %% Plotting shelf export
    if plt_exp
        figure(hf_e)
        no_legend=0;
        for i_param=1:length(param_names)    
            nexttile; hold on; box on
            for i_bnd=1:length(key_param_bnd)
                
                ef_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
                % find profiles that fit into that interval for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bnd
                        % ef_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.QVsfinal(:,i_day);
                        Sref = 35.0;
                        ef_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.QVsfinal(:,i_day).*((Sref-ensemble(i_fjord,i_run).s.Sfinal(:,i_day))/Sref);
                    end
                    
                    Hsill = ensemble(i_fjord,i_run).p.Hsill;
                    Hgl   = ensemble(i_fjord,i_run).p.Hgl;
                    H     = ensemble(i_fjord,i_run).p.H;
                    has_sill = ensemble(i_fjord,i_run).p.sill;
                end
                vline(0,'linewidth',0.5,'linestyle',':','color',[0.7 0.7 0.7])
                
                % take mean and min/max for that subset
                % depths = -ensemble(i_fjord,1).s.z;
                depths = -res_box(i_fjord).zf;
                efmean = mean(ef_ensemble,2,'omitnan');
                % efmin  = min(tf_ensemble,[],2,'omitnan');
                % efmax  = max(tf_ensemble,[],2,'omitnan');
        
                % plot
                % y2 = [depths; flip(depths)];
                % inBetween = [efmin; flip(tfmax)];
                % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor',lcolor(i_param,:),'facealpha',0.1);
                % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor','none','facealpha',0.1);
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
            % xlim([-2 5])
            ylim([-H 0])
            if no_legend==1
                handle_fjords = [handle_fjords hp];
                lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
            end
            if i_param==1
                % text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
                text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
            end
        end
        if no_legend==1
            legend(gca,handle_fjords,param_names,'fontsize',10,'Location','best');
            % hl{i_bnd}.NumColumns=2;
            no_legend = 0;
        end
    end

    
end
figure(hf_t)
% nexttile; axis off
% hl.Layout.Tile='southeast';
xlabel(ht_t,'Temperature (^oC)','fontsize',14);
ylabel(ht_t,'Depth (m)','fontsize',14);
ht_t.TileSpacing='compact';
ht_t.Padding='compact';

if plt_salt
    figure(hf_s)
    xlabel(ht_s,'Salinity','fontsize',14);
    ylabel(ht_s,'Depth (m)','fontsize',14);
    ht_s.TileSpacing='compact';
    ht_s.Padding='compact';
else
    hf_s = [];
end

if plt_exp
    figure(hf_e)
    % xlabel(ht_e,'Shelf-fjord volume flux (m^3s^{-1})','fontsize',14);
    xlabel(ht_e,'Shelf-fjord freshwater flux (m^3s^{-1})','fontsize',14);
    ylabel(ht_e,'Depth (m)','fontsize',14);
    ht_s.TileSpacing='compact';
    ht_s.Padding='compact';
else
    hf_e = [];
end
end