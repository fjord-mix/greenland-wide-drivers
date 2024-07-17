function [hf_t,hf_s] = plot_sensitivity_profiles_v2(X,ensemble,res_box,param_names,i_day,plt_salt,figs_path,i_yr)

if nargin < 5, i_day=1; end
if nargin < 6, plt_salt=0; end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end
n_runs = size(X,1);

%% Finding the low/mid/high ranges for the different parameters
key_param_bnd = {'low ','mid ','high '};

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
hf_t = cell([3,1]);
ht_t = cell(size(hf_t));
hl   = cell(size(hf_t));
if plt_salt
    hf_s = cell(size(hf_t));
    ht_s = cell(size(hf_t));
end

handle_fjords = [];
lcolor = lines(length(param_names));
lbl_fjords = cell([1,length(param_names)]);

for i_bnd = 1:length(key_param_bnd)
hf_t{i_bnd} = figure('Name',sprintf('Temperature sensitivity %s',key_param_bnd{i_bnd}),'Position',[40 40 900 250*length(param_names)]);
ht_t{i_bnd} = tiledlayout('flow');
% ht_t{i_bnd} = tiledlayout(size(ensemble,1)/2,2);

if plt_salt
    hf_s{i_bnd} = figure('Name',sprintf('Salinity sensitivity %s',key_param_bnd{i_bnd}),'Position',[40 40 900 250*length(param_names)]);
    ht_s{i_bnd} = tiledlayout('flow');
end
for i_fjord=1:size(ensemble,1)
    %% Plotting temperature

    figure(hf_t{i_bnd})
    nexttile(i_fjord); hold on; box on
    
    for i_param=1:length(param_names)    
        tf_ensemble = NaN([length(ensemble(i_fjord,1).s.z),size(ensemble,2)]);
        % find profiles that fit into that interval for that cast
        for i_run=1:size(ensemble,2)
            if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
            if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bnd
                tf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Tfinal(:,i_day);
            end
        end
        
        % take mean and min/max for that subset
        depths = -ensemble(i_fjord,1).s.z;
        tfmean = mean(tf_ensemble,2,'omitnan');
        tfmin  = min(tf_ensemble,[],2,'omitnan');
        tfmax  = max(tf_ensemble,[],2,'omitnan');

        % plot
        y2 = [depths; flip(depths)];
        inBetween = [tfmin; flip(tfmax)];
        % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor',lcolor(i_param,:),'facealpha',0.1);
        hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor','none','facealpha',0.1);
        plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:));
        
        % add depictions of GL and sill depths
        base_gl_and_sill_t = 4;
        scatter(base_gl_and_sill_t,-ensemble(i_fjord,1).p.Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
        if ensemble(i_fjord,1).p.sill
            plot([base_gl_and_sill_t base_gl_and_sill_t],[-ensemble(i_fjord,1).p.H -ensemble(i_fjord,1).p.Hsill],'-','linewidth',2,'color',[0 0 0])
        end
        set(gca,'fontsize',14)
        xlim([-2 5])
        ylim([-ensemble(i_fjord,1).p.H 0])
        if i_fjord==1
            handle_fjords = [handle_fjords hfjd];
            lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
        end
        if i_param==1
            % text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
            text(0.02,0.09,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',14)
        end
    end
    if i_fjord==1
        legend(gca,handle_fjords,strcat(key_param_bnd{i_bnd},param_names),'fontsize',12,'Location','West');
        % hl{i_bnd}.NumColumns=2;
    end

    %% Plotting salinity
    if plt_salt
        figure(hf_s{i_bnd})
        nexttile(i_fjord); hold on; box on
        
        for i_param=1:length(param_names)    
            sf_ensemble = NaN([length(ensemble(i_fjord,1).s.z),size(ensemble,2)]);
            % find profiles that fit into that interval for that cast
            for i_run=1:size(ensemble,2)
                if isempty(ensemble(i_fjord,i_run).s), continue; end % we skip any empty entries
                if mask_bnds(i_fjord,i_run).(param_names{i_param}) == i_bnd
                    sf_ensemble(:,i_run) = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                end
            end
            
            % take mean and min/max for that subset
            depths = -ensemble(i_fjord,1).s.z;
            sfmean = mean(sf_ensemble,2,'omitnan');
            sfmin  = min(sf_ensemble,[],2,'omitnan');
            sfmax  = max(sf_ensemble,[],2,'omitnan');
    
            % plot
            y2 = [depths; flip(depths)];
            inBetween = [sfmin; flip(sfmax)];
            % hfjd = patch(inBetween, y2, lcolor(i_param,:),'edgecolor',lcolor(i_param,:),'facealpha',0.1);
            patch(inBetween, y2, lcolor(i_param,:),'edgecolor','none','facealpha',0.1);
            plot(sfmean,depths,'linewidth',1.5,'color',lcolor(i_param,:));
            
            % add depictions of GL and sill depths
            base_gl_and_sill_s = 33;
            scatter(base_gl_and_sill_s,-ensemble(i_fjord,1).p.Hgl,40,'v','filled','MarkerFaceColor',[0 0 0])
            if ensemble(i_fjord,1).p.sill
                plot([base_gl_and_sill_s base_gl_and_sill_s],[-ensemble(i_fjord,1).p.H -ensemble(i_fjord,1).p.Hsill],'-','linewidth',2,'color',[0 0 0])
            end
            set(gca,'fontsize',14)
            xlim([30.5 35])
            ylim([-ensemble(i_fjord,1).p.H 0])
            % if i_fjord==1
            %     handle_fjords = [handle_fjords hfjd];
            %     lbl_fjords{i_fjord} = param_names{i_param};%res_box(i_fjord).id;
            % end
            if i_param==1
                % text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
                text(0.02,0.09,sprintf("(%s) %s",res_box(i_fjord).id,res_box(i_fjord).name),'units','normalized','fontsize',14)
            end
        end
        if i_fjord==1
            legend(gca,handle_fjords,strcat(key_param_bnd{i_bnd},param_names),'fontsize',12,'Location','West');
            % hl{i_bnd}.NumColumns=2;
        end
    end

    
end
figure(hf_t{i_bnd})
% nexttile; axis off
% hl.Layout.Tile='southeast';
xlabel(ht_t{i_bnd},'Temperature (^oC)','fontsize',14);
ylabel(ht_t{i_bnd},'Depth (m)','fontsize',14);
ht_t{i_bnd}.TileSpacing='compact';
ht_t{i_bnd}.Padding='compact';

if plt_salt
    figure(hf_s{i_bnd})
    xlabel(ht_s{i_bnd},'Salinity','fontsize',14);
    ylabel(ht_s{i_bnd},'Depth (m)','fontsize',14);
    ht_s{i_bnd}.TileSpacing='compact';
    ht_s{i_bnd}.Padding='compact';
end

if nargin > 7
    figure(hf_t{i_bnd})
    exportgraphics(gcf,[figs_path,'sensitivity_profiles_temp_',key_param_bnd{i_bnd},'_',num2str(2015+i_yr),'_n',num2str(n_runs),'.png'],'Resolution',300)

    if plt_salt
        figure(hf_s{i_bnd})
        exportgraphics(gcf,[figs_path,'sensitivity_profiles_salt_',key_param_bnd{i_bnd},'_',num2str(2015+i_yr),'_n',num2str(n_runs),'.png'],'Resolution',300)
    end
end
end


end