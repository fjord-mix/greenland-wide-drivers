%% Driver file for simulating the compiled fjords
clearvars
run setup_paths % Configuring paths 

which_year    = 2020;
plot_ensemble = 0;   % whether we want the ensemble to be plotted at the end of `run_model_loop_for_year`
n_pts        = 10; % number of runs per fjord
dt_in_h       = 1;   % model time step in hours
n_years       = 10;  % how many years we want to run
tgt_days      = [n_years*365-180,n_years*365-105];  % which days of the run we want vertical profiles for

time_start = datetime(which_year,01,01);
time_end   = datetime(which_year,12,31);

input_dt   = 30;          % time step of model input (in days)
model_dt   = dt_in_h/24.; % time step in days for the model (e.g., 2h/24 ~ 0.083 days)

t = 0:model_dt:365*n_years; % we want to repeat the yearly data <n_years> times
iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile

%% Compiling data 

file_fjords_compiled = [data_path,'/martim/fjords_digitisation/fjords_gl_sill_depths_reduced_v2.xlsx'];
folder_ctd_casts     = [data_path,'/greenland_common/obs/OMG_all_casts'];
file_lengths         = [data_path,'/martim/fjords_digitisation/fjords_centreline.shp'];
file_fjords          = [data_path,'/martim/fjords_digitisation/fjords_grl.shp'];

fjords_digitised  = shaperead(file_fjords);
fjords_centreline = shaperead(file_lengths);

fjord_matrix = readtable(file_fjords_compiled);
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

%% Define parameter space
param_names  = {'A0','wp','C0'};
param_units  = {'m^2','m','s^{-1}'};
which_fjords = {'0','28','89'};

sermilik_max_bergs = 3e8;  % maximum submerged iceberg area within Sermilik fjord
% sermilik_area      = 1.1850e09; % area of Sermilik fjord (W*L)
% sermilik_vagl      = 7.7028e11; % volume above grounding line of Sermilik fjord (W*L*Hgl)
% iceberg_congestion = sermilik_max_bergs/sermilik_area;

range_params = {[0,1.2*sermilik_max_bergs],...  % A0 (if log scale, starts at 1)
                [10,700],... % wp 
                log10([5e1,5e5])};  % C0

X = NaN([n_pts,length(range_params)]);
for i_param=1:length(range_params)
    X(:,i_param) = linspace(range_params{i_param}(1),range_params{i_param}(2),n_pts);
end
X(:,3) = 10.^X(:,3); % reverting from log quantities to the ones we actually need
range_params{3} = 10.^range_params{3};
disp('Parameter space created.')
% plot_lhs(X,param_names,param_units,1); % quick check of input params distribution

%% Running the model

fjord_model([length(which_fjords),1]) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
dims_ensemble   = [length(which_fjords),n_pts,n_pts,n_pts];
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);

ensemble = struct(ensemble_fields{:});
rmse_tf  = NaN(size(ensemble));

warning('off','all')
for i_fjord=1:size(fjord_model,2)

    % check if we want to run this particular fjord
    for i_find_fjord=1:height(fjord_matrix)
        if strcmp(which_fjords{i_fjord},num2str(fjord_matrix.ID(i_find_fjord))) == 1
            run_fjord=1;
            break
        else
            run_fjord=0;
        end
    end
    if ~run_fjord, continue; end

    % process inputs
    eval(sprintf('id_cast_shelf = num2str(fjord_matrix.shelf_%d(i_fjord));',which_year));
    eval(sprintf('id_cast_fjord = num2str(fjord_matrix.fjord_%d(i_fjord));',which_year));
    
    % default parameters
    p=default_parameters();
    p.fixedthickness=1;
    p.N=60;
    p.dt=model_dt;
    p.t_save = t(1):1:t(end); % save at daily resolution

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
    
    % check if we need to extend the profiles downwards
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
    fprintf('Inputs processing complete.\n')

    i_run=0;
    % loop over desired model parameters
    for i_A0=1:n_pts
    for i_wp=1:n_pts
    for i_C0=1:n_pts
        i_run = i_run+1;
        
        cur_fjord = fjord_model(i_fjord);

        % setting up the specific param combination for this run
        cur_fjord.p.A0=X(i_A0,1);
        cur_fjord.p.wp=X(i_wp,2);
        cur_fjord.p.C0=X(i_C0,3);
        cur_fjord.a.I0 = cur_fjord.p.A0*iceberg_fun(cur_fjord.p.nu0, abs(cur_fjord.p.Hgl), -cumsum(cur_fjord.a.H0)+cur_fjord.a.H0/2);
    
        % running the model
        try
            cur_fjord.s = run_model(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
            fprintf('Fjord %s run %d complete. ',cur_fjord.m.ID,i_run)
        
            % getting only the interesting inputs (to save disk space)
            ensemble(i_fjord,i_A0,i_wp,i_C0).p   = cur_fjord.p;
            ensemble(i_fjord,i_A0,i_wp,i_C0).t   = cur_fjord.t;
            ensemble(i_fjord,i_A0,i_wp,i_C0).m   = cur_fjord.m;
            ensemble(i_fjord,i_A0,i_wp,i_C0).s.t = cur_fjord.s.t;
            ensemble(i_fjord,i_A0,i_wp,i_C0).s.Qsg_max = max(cur_fjord.s.Qsg);
    
            % 10-day avg centered at the target day
            tgt_day = tgt_days(2);
            ensemble(i_fjord,i_A0,i_wp,i_C0).s.Tfinal = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)),2);
            ensemble(i_fjord,i_A0,i_wp,i_C0).s.Tforc  = mean(cur_fjord.s.Ts(:,(tgt_day-5:tgt_day+5)),2); 
            
            zf_obs = cur_fjord.c.zf';
            tf_obs = fjord_model(i_fjord).c.tf;
            
            z_box = -cur_fjord.s.z;
            Tregular = interp1(z_box,ensemble(i_fjord,i_A0,i_wp,i_C0).s.Tfinal,zf_obs,'nearest','extrap');
            Tobs = interp1(zf_obs,tf_obs,z_box,'nearest');
            
            min_depth_rmse = 0;

            % using observation profile as reference
            depths_rmse = (zf_obs > min_depth_rmse) & (zf_obs < ensemble(i_fjord,i_A0,i_wp,i_C0).p.Hgl);
            tf_obs_ref = tf_obs(depths_rmse);
            tf_box_ref = Tregular(depths_rmse);

            % using FjordRPM profile as reference
            % depths_rmse = (z_box > min_depth_rmse) & (z_box < ensemble(i_fjord,i_A0,i_wp,i_C0).p.Hgl);
            % tf_obs_ref = Tobs(depths_rmse);
            % tf_box_ref = ensemble(i_fjord,i_A0,i_wp,i_C0).s.Tfinal(depths_rmse);
            
            rmse_tf(i_fjord,i_A0,i_wp,i_C0) = rmse(tf_box_ref',tf_obs_ref,'omitnan');
            
            fprintf('Output interpolation complete. \n')
        catch ME
            fprintf('fjord %s run %d failed: %s. ',cur_fjord.m.ID,i_run,ME.message)
            cur_fjord.m.error = ME.message;
            cur_fjord.m.stack = ME.stack;
            fjords_crashed{end+1} = cur_fjord;
        end
    
    end % i_A0
    end % i_wp
    end % i_C0
end % fjords

%% saving the inputs and outputs

file_out = [outs_path,'rpm_GRL_misfit_fjords_',num2str(n_pts),'pts_',num2str(which_year),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h.mat'];
save(file_out,'-v7.3','X','param_names','param_units','range_params','fjord_matrix','ensemble','fjord_model','rmse_tf');


%% Plotting misfits

hf = figure('Name','Mistfit per parameter','Position',[40 40 600 500]);
ht = tiledlayout(3,1);
% for i_fjord=1:1

    % Plot A0
    nexttile(1+i_fjord-1); hold on; box on;
    for i_wp=1:n_pts
        for i_C0=1:n_pts
            plot(X(:,1),squeeze(rmse_tf(i_fjord,:,i_wp,i_C0)),'linestyle','-')
        end
    end
    text(0.01,1.01,sprintf("%s",fjord_model(i_fjord).m.name{1}),'horizontalAlignment','left','verticalAlignment','bottom','units','normalized','fontsize',16)
    xlabel('A0')
    set(gca,'XScale','log','fontsize',16)
    
    % Plot wp
    nexttile(2+i_fjord-1); hold on; box on;
    for i_A0=1:n_pts
        for i_C0=1:n_pts
            plot(X(:,2),squeeze(rmse_tf(i_fjord,i_A0,:,i_C0)),'linestyle','-')
        end
    end
    xlabel('wp')
    set(gca,'fontsize',16)

    % Plot C0
    nexttile(3+i_fjord-1); hold on; box on;
    for i_A0=1:n_pts
        for i_wp=1:n_pts
            plot(X(:,3),squeeze(rmse_tf(i_fjord,i_A0,i_wp,:)),'linestyle','-')
        end
    end
    xlabel('C0')
    set(gca,'XScale','log','fontsize',16)
    
% end
ylabel(ht,'RMSE (^oC)','fontsize',16);