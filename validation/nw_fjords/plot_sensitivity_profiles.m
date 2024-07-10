function [hf_t,hf_s] = plot_sensitivity_profiles(X,ensemble,res_box,param_names,param_ranges,i_day)

if nargin < 5, i_day=1; end
if size(X,2) ~=length(param_names), error('input parameter matrix must me [n_runs,n_params], and param_names must have entries for each param'); end

%% Getting how many indices the ensemble has
% the section below should yield something like 'i1,i2,i3,i4,i5'
% n_indices=length(param_names);
% indices=[];
% for i_dim=1:n_indices-1
%     indices=[indices,'i',num2str(i_dim),','];
% end
% indices=[indices,'i',num2str(i_dim+1)];
% indices_split=strsplit(indices,',');
% for i_dim=1:length(indices_split)
%     eval(sprintf("%s=ceil(length(param_ranges(%d))/2.);",indices_split{i_dim},i_dim));
% end

%% Finding the low/mid/high ranges for the different parameters
key_param_bnd = {'low','mid','high'};

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
% hf_s = cell(size(hf_t));
% ht_s = cell(size(hf_t));

handle_fjords = [];
lcolor = lines(size(ensemble,1));
lbl_fjords = cell([1,size(ensemble,1)]);

for i_bnd = 1:length(key_param_bnd)
hf_t{i_bnd} = figure('Name',sprintf('Temperature sensitivity %s',key_param_bnd{i_bnd}),'Position',[40 40 900 250*length(param_names)]);
ht_t{i_bnd} = tiledlayout('flow');
% hf_s{i_bnd} = figure('Name',sprintf('Salinity sensitivity %s',key_param_bnd{i_bnd}),'Position',[40 40 900 250*length(param_names)]);
% ht_s{i_bnd} = tiledlayout('flow');
for i_param=1:length(param_names)
    nexttile; hold on; box on

    %% Plotting temperature
    
    for i_fjord=1:size(ensemble,1)
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
        hfjd = patch(inBetween, y2, lcolor(i_fjord,:),'edgecolor','none','facealpha',0.1);
        plot(tfmean,depths,'linewidth',1.5,'color',lcolor(i_fjord,:));
        
        % add depictions of GL and sill depths
        base_gl_and_sill_t = 2.5;
        scatter(base_gl_and_sill_t-0.1*i_fjord,-ensemble(i_fjord,1).p.Hgl,40,'v','filled','MarkerFaceColor',lcolor(i_fjord,:))
        if ensemble(i_fjord,1).p.sill
            plot([base_gl_and_sill_t-0.1*i_fjord base_gl_and_sill_t-0.1*i_fjord],[-ensemble(i_fjord,1).p.H -ensemble(i_fjord,1).p.Hsill],'-','linewidth',2,'color',lcolor(i_fjord,:))
        end
        set(gca,'fontsize',14)
        xlim([-2 5])
        % ylim([-ensemble(i_fjord,1).p.H 0])
        if i_param==1
            handle_fjords = [handle_fjords hfjd];
            lbl_fjords{i_fjord} = res_box(i_fjord).id;
        end
        if i_fjord==1
            text(0.02,0.05,sprintf("%s = [%.1e,%.1e]",param_names{i_param},param_bnds(i_bnd,i_param),param_bnds(i_bnd+1,i_param)),'units','normalized','fontsize',14)
        end
    end


end
figure(hf_t{i_bnd})
% nexttile; axis off
hl = legend(gca,handle_fjords,lbl_fjords,'fontsize',12,'Location','West');
hl.NumColumns=2;
% hl.Layout.Tile='southeast';
xlabel(ht_t{i_bnd},'Temperature (^oC)','fontsize',14);
ylabel(ht_t{i_bnd},'Depth (m)','fontsize',14);
ht_t{i_bnd}.TileSpacing='compact';
ht_t{i_bnd}.Padding='compact';

end

% handle_fjords = [handle_fjords havg hbnd];
% handle_fjords_s = [handle_fjords_s havg_s hbnd_s];
% lbl_fjords = lbl_fjords(~cellfun(@isempty,lbl_fjords));
% lbl_fjords{end+1} = 'Mean';
% lbl_fjords{end+1} = 'Spread';


end