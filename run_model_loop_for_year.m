function [file_out,tgt_days] = run_model_loop_for_year(which_year,fjords_digitised,fjords_centreline,fjord_matrix,folder_ctd_casts,X,param_names,n_years,tgt_days,dt_in_h,dt_plume_h,plot_ensemble)

warning('off','all')
run setup_paths % Configuring paths
% letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
clear_empty = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array

%% Initialise all needed variables

if nargin < 1
    which_year = 2016;
end
time_start = datetime(which_year,01,01);
time_end   = datetime(which_year,12,31);

n_runs = size(X,1);
input_dt   = 30;          % time step of model input (in days)
if isempty(dt_in_h), dt_in_h    = 3.; end          % time step in hours
model_dt   = dt_in_h/24.; % time step in days for the model (e.g., 2h/24 ~ 0.083 days)
n_years    = 10;           % how many years we want to run
tgt_days   = [n_years*365-60,n_years*365-90];  % which days of the run we want vertical profiles for
name_days  = {'post','peak'};

t = 0:model_dt:365*n_years; % we want to repeat the yearly data <n_years> times

iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile



if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(height(fjord_matrix)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);

for i_fjord=1:height(fjord_matrix)
    eval(sprintf('id_cast_shelf = num2str(fjord_matrix.shelf_%d(i_fjord));',which_year));
    eval(sprintf('id_cast_fjord = num2str(fjord_matrix.fjord_%d(i_fjord));',which_year));
    if ~strcmp(id_cast_shelf,'NaN') && ~strcmp(id_cast_fjord,'NaN')
        % digitised_id = find([fjords_digitised.id] == fjord_matrix.ID(i_fjord));
        % centreline_id = find([fjords_centreline.id] == fjord_matrix.ID(i_fjord));

        % x_fjord = fjords_digitised(digitised_id).x;
        % y_fjord = fjords_digitised(digitised_id).y;

        % default parameters
        p=default_parameters();
        p.fixedthickness=1;
        p.N=60;
        p.dt=model_dt;
        p.t_save = t(1):1:t(end); % save at daily resolution
        p.run_plume_every = floor(24*dt_plume_h/dt_in_h); % run the plume model every dt_plume_h hours

        % find the ocean forcing and fjord profile
        omg_data_shelf = dir([folder_ctd_casts,'/*',id_cast_shelf,'*.nc']);
        omg_data_fjord = dir([folder_ctd_casts,'/*',id_cast_fjord,'*.nc']);

        % get the metadata
        m.ID = num2str(fjord_matrix.ID(i_fjord));
        m.name = fjord_matrix.name(i_fjord);
        m.lon=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lon');
        m.lat=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lat');

        % fjord obervations for reference and model fitting
        tf = movmean(ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'temperature'),5);
        sf = movmean(ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'salinity'),5);
        zf = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'depth');
        is_downcast = zeros(size(zf));
        max_depth=0;
        for k=1:length(is_downcast)
            if zf(k) > max_depth 
                is_downcast(k) = 1;
                max_depth = zf(k);
            end
        end

        tf_nonans = tf(~isnan(tf) & ~isnan(sf) & is_downcast);
        sf_nonans = sf(~isnan(tf) & ~isnan(sf) & is_downcast);
        zf_nonans = zf(~isnan(tf) & ~isnan(sf) & is_downcast);
        if isempty(tf_nonans) || isempty(sf_nonans) % we disregard this fjord because we actually do not have the data
            fprintf('Fjord %s has only NaNs in temperature or salinity!\n',m.ID)
            continue
        end
        
        % get geometry
        p.L = fjords_centreline([fjords_centreline.id] == fjord_matrix.ID(i_fjord)).length.*1e3;
        p.W = fjords_digitised([fjords_digitised.id] == fjord_matrix.ID(i_fjord)).area/p.L;
        p.Hsill = abs(fjord_matrix.sill_depth(i_fjord));
        p.Hgl   = abs(fjord_matrix.gl_depth(i_fjord));

        if p.Hsill == -999 % if the fjord has no sill
            p.sill=0;        % we ignore the sill
            p.Hsill = p.Hgl; % just so the RPM initial checks won't fail
        else
            p.sill=1;
        end
        
        p.H = max(p.Hgl,p.Hsill); % we treat the fjord to be as deep as the deepest point we know exists in the fjord
        if max(zf) > p.H % if the cast inside the fjord goes deeper than what we set as current depth, then the fjord is as deep as the cast
            p.H = ceil(max(zf));
        end

        if p.H == p.Hsill,p.sill=0; end % if we made the fjord as deep as the sill, it essentially has no sill

        if (p.H-p.Hsill < 10)                     % if p.H and p.Hsill are too close to each other, 
            p.Hsill = p.Hsill-(10-(p.H-p.Hsill)); % we bump the sill upwards a wee bit (workaround so the layer division below the sill works)
        end
        if p.H < p.Hsill % the sill cannot be deeper than the fjord itself
            p.Hsill = p.H;
            p.sill = 0;
        end
        
        % get subglacial discharge (separate file for better readability
        run set_subglacial_discharge.m
        if max(f.Qsg) > 4e3
            disp('Subglacial runoff data is abnormaly high! Skipping...')
            continue
        end

        % get the shelf-profile forcing
        z_shelf = ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'depth');
        t_shelf = movmean(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'temperature'),5);
        s_shelf = movmean(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'salinity'),5);
        % removes upcast and fluctuations at subsurface (e.g., rosette bobbing up and down with the swell)
        is_downcast = zeros(size(z_shelf));
        max_depth=0;
        for k=1:length(is_downcast)
            if z_shelf(k) > max_depth 
                is_downcast(k) = 1;
                max_depth = z_shelf(k);
            end
        end
        
        % remove NaNs from cast
        z_shelf_nonans = z_shelf(~isnan(t_shelf) & ~isnan(s_shelf) & is_downcast);
        t_shelf_nonans = t_shelf(~isnan(t_shelf) & ~isnan(s_shelf) & is_downcast);
        s_shelf_nonans = s_shelf(~isnan(t_shelf) & ~isnan(s_shelf) & is_downcast);
        
        % we only carry on if the shelf data actually exists
        if ~isempty(t_shelf_nonans) && ~isempty(s_shelf_nonans)

            % check if we need to extend the profiles downwards
            % [~,i_bottom] = max(z_shelf_nonans);
            if max(z_shelf_nonans) < p.H 
                dz_missing = p.H - max(z_shelf_nonans);
                zs = [z_shelf_nonans; max(z_shelf_nonans)+dz_missing];
                Ts = interp1(z_shelf_nonans,t_shelf_nonans,zs,'nearest','extrap');
                Ss = interp1(z_shelf_nonans,s_shelf_nonans,zs,'nearest','extrap');
            else
                zs = z_shelf_nonans;
                Ts = t_shelf_nonans;
                Ss = s_shelf_nonans;
            end

            % check if we need to extend the profiles upwards (up to 1m depth)
            if min(zs) > 1
                if min(zs) > 50 % however, we do not want to excessively interpolate data
                    fprintf('Fjord %s has no data until %.2f m. Skipping...\n',m.ID,min(zs))
                    continue
                end
                zs_with_surf = [1; zs];
                Ts = interp1(zs,Ts,zs_with_surf,'nearest','extrap');
                Ss = interp1(zs,Ss,zs_with_surf,'nearest','extrap');
                zs = zs_with_surf;
            end
        
            zs = flip(zs);
            Ts = flip(Ts);
            Ss = flip(Ss);
            f.zs = -zs;
            f.ts = [t(1),t(end)];
            f.Ts = repmat(Ts,1,length(f.ts));
            f.Ss = repmat(Ss,1,length(f.ts));

            % icebergs
            % nu0 controls how much the total iceberg area extends below the surface
            % we do not vary this parameter, as previous sensitivity
            % analyses showed it to have a smaller impact than A0 or M0
            p.nu0=25;
    
            % initial conditions
            a.H0 = ones([p.N],1).*p.H/p.N;
            ints = [0;cumsum(a.H0)];
            z_box = 0.5*(ints(1:end-1)+ints(2:end));
            [T_clim, S_clim, ~] = bin_forcings(f, a.H0, t);
            a.T0 = T_clim(:,1);
            a.S0 = S_clim(:,1);
            if p.Hsill < p.H % we ignore the shelf part in the initial conditions, and extrapolate with values at sill depth instead
                % figure; hold on; plot(a.T0,-z_box); hline(-p.Hsill);
                [~,i_sill] = min(abs(p.Hsill-z_box));
                a.T0(i_sill+1:end) = a.T0(i_sill);
                a.S0(i_sill+1:end) = a.S0(i_sill);
                % plot(a.T0,-z_box)
            end
    
    
            % create the model structure
            fjord_model(i_fjord).t = t;
            fjord_model(i_fjord).m = m;
            fjord_model(i_fjord).p = p;
            fjord_model(i_fjord).a = a;
            fjord_model(i_fjord).f = f;
            fjord_model(i_fjord).c.tf = tf_nonans;
            fjord_model(i_fjord).c.sf = sf_nonans;
            fjord_model(i_fjord).c.zf = zf_nonans;
            fjord_model(i_fjord).c.ts = Ts;
            fjord_model(i_fjord).c.ss = Ss;
            fjord_model(i_fjord).c.zs = zs;
            fprintf('Done pre-processing fjord %d.\n',fjord_matrix.ID(i_fjord))
            clear m p a f
        else
            fprintf('Fjord %d has no shelf forcing data. Skipping...\n',fjord_matrix.ID(i_fjord))
        end

    else
        fprintf('No %d data for fjord %d.\n',which_year,fjord_matrix.ID(i_fjord))
    end
end
idx = arrayfun(clear_empty,fjord_model);
fjord_model(idx)=[]; % remove the empty elements
n_fjords = length(fjord_model);
fprintf('Inputs processing complete. Total of fjords for %d: %d\n',which_year,n_fjords)

%% Running the model
warning('on','all')
dims_ensemble = [n_fjords,n_runs];
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);
ensemble = struct(ensemble_fields{:});
ensemble(dims_ensemble) = struct("p",[],"t",[],"m",[],"s",[]);
fjords_crashed = {};

parpool(3)
% run_counter=0;
for i_fjord=1:length(fjord_model)
tic
% for i_run=1:n_runs
parfor i_run=1:n_runs
    % run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    
    % Sets the model parameters according to the LHS we created
    for i_param=1:size(X,2)
        cur_fjord.p.(param_names{i_param}) = X(i_run,i_param);
        % if strcmp(param_names{i_param},'A0') % switching from "congestion factor" to iceberg submerged area
        %     cur_fjord.p.(param_names{i_param}) = X(i_run,i_param).*cur_fjord.p.L.*cur_fjord.p.W.*cur_fjord.p.Hgl;
        % end
    end
    cur_fjord.a.I0 = cur_fjord.p.A0*iceberg_fun(cur_fjord.p.nu0, abs(cur_fjord.p.Hgl), -cumsum(cur_fjord.a.H0)+cur_fjord.a.H0/2);
    
    % cur_fjord.p.plot_runtime=1; % for debuggging purposes if needed
        try
            cur_fjord.s = run_model(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
            fprintf('Fjord %d run %d complete. ',i_fjord,i_run)
    
            % we will not save the entire solution for disk management purposes
            ensemble(i_fjord,i_run).p   = cur_fjord.p;
            ensemble(i_fjord,i_run).t   = cur_fjord.t;
            ensemble(i_fjord,i_run).m   = cur_fjord.m;
            ensemble(i_fjord,i_run).s.t = cur_fjord.s.t;
            ensemble(i_fjord,i_run).s.Qsg_max = max(cur_fjord.s.Qsg);

            Tfinal  =NaN([cur_fjord.p.N,length(tgt_days)]);
            Sfinal  =NaN(size(Tfinal));
            QVsfinal=NaN(size(Tfinal));
            QVpfinal=NaN(size(Tfinal));
            QMpfinal=NaN(size(Tfinal));
            
                         % we use a constant value for an easier comparison between different fjords. 
            Sref = 35.0; % If we were concerned about the FW quantity itself, this should be a fjord-specific value 
                         % (e.g., shelf salinity at sill depth). Varing
                         % Sref by 2 PSU alters the results by < ~5%

            fw_export=NaN([1,length(tgt_days)]);
            i_max_export=NaN(size(fw_export));
            inb=NaN(size(i_max_export));
            for i_day=1:length(tgt_days)
                % 10-day avg centered at the target day
                tgt_day = tgt_days(i_day);
                Tfinal(:,i_day)   = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)),2); 
                Sfinal(:,i_day)   = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)),2);
                QVsfinal(:,i_day) = mean(cur_fjord.s.QVs(:,(tgt_day-5:tgt_day+5)),2);
                QVpfinal(:,i_day) = mean(squeeze(cur_fjord.s.QVp(:,(tgt_day-5:tgt_day+5))),2);
                QMpfinal(:,i_day) = mean(squeeze(cur_fjord.s.QMp(:,(tgt_day-5:tgt_day+5))),2);
                fw_export(i_day)  = sum(QVsfinal(:,i_day).*((Sref-Sfinal(:,i_day))/Sref));
                inb(i_day)        = ceil(mean(cur_fjord.s.knb(tgt_day-5:tgt_day+5)));
            end
            QVs_mean = mean(cur_fjord.s.QVs(:,(end-364:end)),2,'omitnan');
            QVp_mean = mean(squeeze(cur_fjord.s.QVp(:,(end-364:end))),2,'omitnan');
            QMp_mean = mean(squeeze(cur_fjord.s.QMp(:,(end-364:end))),2,'omitnan');
            S_mean   = mean(cur_fjord.s.S(:,(end-364:end)),2,'omitnan');
            T_mean   = mean(cur_fjord.s.T(:,(end-364:end)),2,'omitnan');
            fw_mean_export_profile = QVs_mean.*((Sref-S_mean)/Sref);
            fw_mean_discharge_profile = QVp_mean.*((Sref-S_mean)/Sref);

            % we need to ignore the cases when there is no plume (i.e., knb_t==0 and inb==0)
            knb_t = cur_fjord.s.knb;
            inb_mask = ones(size(inb));
            knb_mask = ones(size(knb_t));
            inb_mask(inb==0)   = 0;
            knb_mask(knb_t==0) = 0;
            inb(inb==0)        = 1; % so it doesnt crash when getting the depths
            knb_t(knb_t==0)    = 1;

            ensemble(i_fjord,i_run).s.Tfinal   = Tfinal;
            ensemble(i_fjord,i_run).s.Sfinal   = Sfinal;
            ensemble(i_fjord,i_run).s.QVsfinal = QVsfinal;
            ensemble(i_fjord,i_run).s.QVpfinal = QVpfinal;
            ensemble(i_fjord,i_run).s.QMpfinal = QMpfinal;
            ensemble(i_fjord,i_run).s.Tmean    = T_mean;
            ensemble(i_fjord,i_run).s.Smean    = S_mean;
            ensemble(i_fjord,i_run).s.QMpmean  = QMp_mean;

            ensemble(i_fjord,i_run).s.Tforc = mean(cur_fjord.s.Ts(:,(tgt_day-5:tgt_day+5)),2); 
            ensemble(i_fjord,i_run).s.Sforc = mean(cur_fjord.s.Ss(:,(tgt_day-5:tgt_day+5)),2); 

            ymean_fw_export = mean(cur_fjord.s.QVs(:,end-364:end).*((Sref-cur_fjord.s.S(:,end-364:end))/Sref),2,'omitnan');
            [~,i_max_export] = min(ymean_fw_export); % we use "min" because QVs < 0 means water is leaving the layer towards the shelf

            fw_export_t  = sum(cur_fjord.s.QVs.*((Sref-cur_fjord.s.S)/Sref),1,'omitnan');
            [~,i_max_export_t] = min(cur_fjord.s.QVs.*((Sref-cur_fjord.s.S)/Sref),[],1,'omitnan');
            
            ensemble(i_fjord,i_run).s.fw_profile_export    = fw_mean_export_profile;
            ensemble(i_fjord,i_run).s.fw_profile_discharge = fw_mean_discharge_profile;
            ensemble(i_fjord,i_run).s.fw_export            = fw_export;
            ensemble(i_fjord,i_run).s.fw_export_t          = fw_export_t;
            ensemble(i_fjord,i_run).s.z_max_export         = cur_fjord.s.z(i_max_export);
            ensemble(i_fjord,i_run).s.z_max_export_t       = cur_fjord.s.z(i_max_export_t);
            ensemble(i_fjord,i_run).s.Qsg                  = cur_fjord.s.Qsg;

            ensemble(i_fjord,i_run).s.znb                  = cur_fjord.s.z(inb);
            ensemble(i_fjord,i_run).s.znb_t                = cur_fjord.s.z(knb_t);
            ensemble(i_fjord,i_run).s.znb(inb_mask==0)        = NaN;
            ensemble(i_fjord,i_run).s.znb_t(knb_mask==0)      = NaN;

            zf_obs = cur_fjord.c.zf';
            z_box = -cur_fjord.s.z;
            Tregular = interp1(z_box,cur_fjord.s.T,zf_obs,'nearest','extrap');
            Sregular = interp1(z_box,cur_fjord.s.S,zf_obs,'nearest','extrap');
            
            range_upper = zf_obs < 50;
            range_inter = zf_obs > 50 & zf_obs < 250;
            range_lower = zf_obs > 250 & zf_obs < 500;
            
            ensemble(i_fjord,i_run).s.z      = z_box;
            ensemble(i_fjord,i_run).s.Tupper = mean(Tregular(range_upper,:),1,'omitnan'); % avg temp of the upper 50 m
            ensemble(i_fjord,i_run).s.Tinter = mean(Tregular(range_inter,:),1,'omitnan'); % avg temp of the 50 - 250 m layer
            ensemble(i_fjord,i_run).s.Tlower = mean(Tregular(range_lower,:),1,'omitnan'); % avg temp of the 250 - 500 m layer

            ensemble(i_fjord,i_run).s.Supper = mean(Sregular(range_upper,:),1,'omitnan'); % avg temp of the upper 50 m
            ensemble(i_fjord,i_run).s.Sinter = mean(Sregular(range_inter,:),1,'omitnan'); % avg temp of the 50 - 250 m layer
            ensemble(i_fjord,i_run).s.Slower = mean(Sregular(range_lower,:),1,'omitnan'); % avg temp of the 250 - 500 m layer
            fprintf('Output interpolation complete. ')
        catch ME
            fprintf('Fjord %d run %d failed: %s. ',i_fjord,i_run,ME.message)
            % cur_fjord.m.error = ME.message;
            % cur_fjord.m.stack = ME.stack;
            % fjords_crashed{end+1} = cur_fjord;
        end
        fprintf('\n')
end % for i_run
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
file_out = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year),'_dtp',num2str(dt_plume_h),'h_dtm',num2str(dt_in_h),'h.mat'];
save(file_out,'-v7.3','ensemble','fjord_model');
if ~isempty(fjords_crashed)
    save([file_out,'_crashed.mat'],'-v7.3','fjords_crashed');
end
disp('Outputs saved.')

%% Sanity-check post processing and plot
% load(file_out);

if plot_ensemble
    [res_obs,res_box] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
    disp('Postprocessing ensemble done.')
    plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2);
    % plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2,[],0,0,0);
    exportgraphics(gcf,[figs_path,'profiles_GRL_temp_',num2str(which_year),'_n',num2str(n_runs),'.png'],'Resolution',300)
    
    % plot_sensitivity_profiles_v3(X,ensemble,res_box,res_obs,param_names,2,0,[],which_year-2015);%, which_fj_sens{which_year-2015});
    % plot_sensitivity_profiles_v3(X,ensemble,res_box,res_obs,param_names,2,0,[],which_year-2015, which_fj_sens{which_year-2015});
    % exportgraphics(gcf,[figs_path,'sensitivity_profiles_temp_',num2str(which_year),'_n',num2str(n_runs),'.png'],'Resolution',300)
end
end