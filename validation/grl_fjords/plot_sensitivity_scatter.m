function [hf_t,hf_s,hf_f,hf_z] = plot_sensitivity_scatter(X,ensemble,res_box,param_names,i_day,plt_salt,plt_fw,plt_fwz,which_fjords)

if nargin < 5, i_day=1; end
if nargin < 6, plt_salt=0; end
if nargin < 7, plt_fw=0; end
if nargin < 8, plt_fwz=0; end
if nargin > 8
    n_fjords = length(which_fjords);
else
    n_fjords = size(ensemble,1);
end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end

%% Finding the value bins for the different parameters
n_runs     = size(X,1);
n_params   = size(X,2);
bins       = 0:10:100;
param_bnds = NaN([length(bins),n_params]);
for i_param=1:n_params
    param_bnds(:,i_param) = prctile(X(:,i_param),bins);
end

% loops over all ensemble entries
struct_fields = param_names;
struct_fields{2,1} = cell(size(ensemble));
mask_bnds = struct(struct_fields{:});
for i_fjord=1:size(ensemble,1)
    for i_run=1:size(ensemble,2)
        % for each parameter, assigns a mask of 1:10 if the value is is in the 1:100th percentile
        if ~isempty(ensemble(i_fjord,i_run).p)
            for i_param=1:n_params
                for i_bin=1:length(bins)-1
                    if ensemble(i_fjord,i_run).p.(param_names{i_param}) > param_bnds(i_bin,i_param) && ensemble(i_fjord,i_run).p.(param_names{i_param}) <= param_bnds(i_bin+1,i_param)
                        mask_bnds(i_fjord,i_run).(param_names{i_param}) = i_bin;
                    end
                end
            end
        end
    end
end

%% Preparing the figures
lcolor = lines(n_params);

hf_t = figure('Name','Temperature sensitivity','Position',[40 40 900 700]);
ht_t = tiledlayout(n_fjords,n_params);

if plt_salt
    hf_s = figure('Name','Salinity sensitivity','Position',[40 40 900 700]);
    ht_s = tiledlayout(n_fjords,n_params);
end
if plt_fw
    hf_f = figure('Name','Freshwater export sensitivity','Position',[40 40 900 700]);
    ht_f = tiledlayout(n_fjords,n_params);
end
if plt_fwz
    hf_z = figure('Name','Freshwater export depth sensitivity','Position',[40 40 900 700]);
    ht_z = tiledlayout(n_fjords,n_params);
end
i_iter=0;
for i_fjord=1:size(ensemble,1)
    if nargin > 8
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
    i_iter=i_iter+1;

    %% Plotting temperature
    figure(hf_t)
    for i_param=1:n_params    
        nexttile; hold on; box on
        for i_bin=1:length(bins)-1
            tf_ensemble = NaN([length(res_box(i_fjord).zf),size(ensemble,2)]);
            tf_full_ens = NaN(size(tf_ensemble));
            % find profiles that fit into that bin for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bin
                    tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,i_day);
                end
                tf_full_ens(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,i_day);
            end
            
            % tfplot = mean(tf_ensemble,[1,2],'omitnan'); % take mean for that subset
            % tfplot = std(mean(tf_ensemble,1,'omitnan'),'omitnan');
            % tfplot = max(mean(tf_ensemble,2,'omitnan'))- min(mean(tf_ensemble,2,'omitnan'));
            tfplot = mean(tf_ensemble,[1,2],'omitnan') - mean(tf_full_ens,[1,2],'omitnan');
            scatter(i_bin*10,tfplot,50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        end
        set(gca,'fontsize',14)
        if i_param==1
            text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
        end
        if i_iter==length(which_fjords)
            xlabel([param_names{i_param},' percentile'],'fontsize',16)
        end
    end

    %% Plotting salinity
    if plt_salt
        figure(hf_s)
        for i_param=1:n_params
            nexttile; hold on; box on
            for i_bin=1:length(bins)-1
                sf_ensemble = NaN([length(ensemble(i_fjord,1).s.z),size(ensemble,2)]);
                sf_full_ens = NaN(size(sf_ensemble));
                % find profiles that fit into that interval for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bin
                        sf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                    end
                    sf_full_ens(:,i_run) = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                end
                
                % sfplot = mean(sf_ensemble,[1,2],'omitnan'); % take mean for that subset
                % sfplot = std(mean(sf_ensemble,1,'omitnan'),'omitnan');
                % sfplot = max(mean(sf_ensemble,2,'omitnan'))- min(mean(sf_ensemble,2,'omitnan'));
                sfplot = 1e3.*(mean(sf_ensemble,[1,2],'omitnan') - mean(sf_full_ens,[1,2],'omitnan'));
                scatter(i_bin*10,sfplot,50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5); % plot
            
            end
            set(gca,'fontsize',14)
            if i_param==1
                text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
            end
            if i_iter==length(which_fjords)
                xlabel([param_names{i_param},' percentile'],'fontsize',16)
            end
        end
    end

    %% Plotting freshwater export
    if plt_fw
        figure(hf_f)
        for i_param=1:n_params    
            nexttile; hold on; box on
            for i_bin=1:length(bins)-1
                fw_ensemble = NaN([1,size(ensemble,2)]);
                fw_full_ens = NaN(size(fw_ensemble));
                % find profiles that fit into that bin for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bin
                        % fw_ensemble(i_run) = ensemble(i_fjord,i_run).s.fw_export(i_day);
                        fw_ensemble(i_run) = mean(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                        % fw_ensemble(i_run) = min(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                    end
                    % fw_full_ens(i_run) = ensemble(i_fjord,i_run).s.fw_export(i_day);
                    fw_full_ens(i_run) = mean(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                    % fw_full_ens(i_run) = min(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                end
                
                fwplot = mean(fw_ensemble,[1,2],'omitnan'); % take mean for that subset
                % fwplot = std(mean(fw_ensemble,1,'omitnan'),'omitnan');
                % fwplot = max(mean(fw_ensemble,2,'omitnan'))- min(mean(fw_ensemble,2,'omitnan'));
                % fwplot = mean(fw_ensemble,[1,2],'omitnan') - mean(fw_full_ens,[1,2],'omitnan');
                scatter(i_bin*10,-1e-3.*fwplot,50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            end
            set(gca,'fontsize',14)
            if i_param==1
                text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
            end
            if i_iter==length(which_fjords)
                xlabel([param_names{i_param},' percentile'],'fontsize',16)
            end
        end
    end
    %% Plotting max export depth
    if plt_fwz
        figure(hf_z)
        for i_param=1:n_params    
            nexttile; hold on; box on
            for i_bin=1:length(bins)-1
                fz_ensemble = NaN([1,size(ensemble,2)]);
                fz_full_ens = NaN(size(fz_ensemble));
                % find profiles that fit into that bin for that cast
                for i_run=1:size(ensemble,2)
                    if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                    if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bin
                        % fz_ensemble(i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_day);
                        fz_ensemble(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
                        % fz_ensemble(i_run) = mean(ensemble(i_fjord,i_run).s.z_max_export_t(end-365:end));
                        % i_fz = find(ensemble(i_fjord,i_run).s.fw_export_t == fw_full_ens(i_run),1,'last');
                        % fz_ensemble(i_run) = ensemble(i_fjord,i_run).s.z_max_export_t(i_fz);
                    end
                    % fz_full_ens(i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_day);
                    fz_full_ens(i_run) = ensemble(i_fjord,i_run).s.z_max_export;
                    % fz_full_ens(i_run) = mean(ensemble(i_fjord,i_run).s.z_max_export_t(end-365:end));
                    % i_fz = find(ensemble(i_fjord,i_run).s.fw_export_t == fw_full_ens(i_run),1,'last');
                    % fz_full_ens(i_run) = ensemble(i_fjord,i_run).s.z_max_export_t(i_fz);
                end
                
                fzplot = mean(fz_ensemble,[1,2],'omitnan'); % take mean for that subset
                % fzplot = std(mean(fz_ensemble,1,'omitnan'),'omitnan');
                % fzplot = max(mean(fz_ensemble,2,'omitnan'))- min(mean(fz_ensemble,2,'omitnan'));
                % fzplot = mean(fz_ensemble,[1,2],'omitnan') - mean(fz_full_ens,[1,2],'omitnan');
                scatter(i_bin*10,round(fzplot,2),50,'filled','MarkerFaceColor',lcolor(i_param,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            end
            set(gca,'fontsize',14)
            % ylim([-150 -10])
            if i_param==1
                text(0.02,1.075,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',12)
            end
            if i_iter==length(which_fjords)
                xlabel([param_names{i_param},' percentile'],'fontsize',16)
            end
        end
    end
end

% ylabel(ht_t,'Mean temperature (^oC)','fontsize',16)
% ylabel(ht_t,'\sigma temperature (^oC)','fontsize',16)
% ylabel(ht_t,'Temperature range (^oC)','fontsize',16)
ylabel(ht_t,'Departure from mean temperature (^oC)','fontsize',16)

if plt_salt
    % ylabel(ht_s,'Mean salinity','fontsize',16)
    % ylabel(ht_s,'\sigma salinity','fontsize',16)
    % ylabel(ht_s,'Salinity range','fontsize',16)
    ylabel(ht_s,'Departure from mean salinity (10^{-3} PSU)','fontsize',16)
end
if plt_fw
    ylabel(ht_f,'Mean freshwater export (mSv)','fontsize',16)
    % ylabel(ht_f,'Peak freshwater export (mSv)','fontsize',16)
end
if plt_fwz
    ylabel(ht_z,'Average maximum export depth (m)','fontsize',16)
    % ylabel(ht_z,'Departure from mean maximum export depth (m)','fontsize',16)
end

end