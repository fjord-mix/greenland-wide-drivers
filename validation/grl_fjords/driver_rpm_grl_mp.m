%% Driver file for FjordRPM forced by MP 0.25 degree resolution
clearvars
warning('off','all')
run setup_paths % Configuring paths

which_year = 2018; % {2016,2017,2018,2019,2020}
time_start = datetime(2013,01,01);
time_end   = datetime(2018,12,31);

n_runs     = 10;         % number of runs per fjord
input_dt   = 30;          % time step of model input (in days)
dt_in_h    = 3.;          % time step in hours
model_dt   = dt_in_h/24.; % time step in days for the model (e.g., 2h/24 ~ 0.083 days)

% tgt_days   = [935,1010];  % which days of the run we want vertical profiles for
tgt_days   = [1280,1400];  % which days of the run we want vertical profiles for
name_days  = {'peak','post'};


runtime_axis = time_start:calendarDuration(0,0,0,dt_in_h,0,0):time_end;
t = [juliandate(time_start):model_dt:juliandate(time_end)] - juliandate(time_start);


fun = @(s) all(structfun(@isempty,s));              % tiny function to get rid of empty entries in array
iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile

%% Define parameter space
param_names = {'A0','wp','C0','K0'};

range_params = {[0,3e8],...    % A0
                [10,700],...   % wp 
                [1e2,5e4],...  % C0
                [1e-4,1e-3]};  % K0

rng('default') % set the seed for reproducibility
uqlab
for i_param=1:length(param_names)
    iOpts.Marginals(i_param).Type       = 'uniform';
    iOpts.Marginals(i_param).Parameters = range_params{i_param};
    iOpts.Marginals(i_param).Bounds     = range_params{i_param};
    iOpts.Marginals(i_param).Name       = param_names{i_param};
end

input = uq_createInput(iOpts); % create probability functions
X = uq_getSample(input,n_runs,'LHS');  % training dataset
disp('Parameter space created.')

%% Compile data

% [datasets,~,~,~,~] = compile_datasets(data_path); % we only want the ocean and Qsg from the compilation
% we set up the dataset structure to use with the data-compilation functions for extracting the ocean forcing
datasets.ocean           = [data_path,'/greenland/ocean/MP_0p25deg_ocn_grl_julian.nc'];
datasets.full_ocn.lat    = ncread(datasets.ocean,'latitude');
datasets.full_ocn.lon    = ncread(datasets.ocean,'longitude'); 
datasets.full_ocn.depth  = ncread(datasets.ocean,'depth');
datasets.full_ocn.time   = ncread(datasets.ocean,'time');
datasets.full_ocn.thetao = ncread(datasets.ocean,'thetao');
datasets.full_ocn.so     = ncread(datasets.ocean,'so');    
datasets.utility.earth_radius = 6378137.0;
datasets.utility.eccentricity = 0.08181919;
datasets.utility.lat_true     = 70;
datasets.utility.lon_posy     = -45;


% We use the manually digitised compilation instead
file_fjords_compiled = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_gl_sill_depths_reduced.xlsx'];
folder_ctd_casts     = [data_path,'/greenland/obs/OMG_all_casts'];
file_lengths = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_centreline.shp'];
file_fjords = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_grl.shp'];

fjords_digitised  = shaperead(file_fjords); % contains centroids and area for each fjord
fjords_centreline = shaperead(file_lengths); % contains centreline length for each fjord

fjord_matrix = readtable(file_fjords_compiled); % contains gate IDs for subglacial discharge, GL depth, and shelf/fjord profiles for each fjord
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(height(fjord_matrix)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);

for i_fjord=1:height(fjord_matrix)
    eval(sprintf('id_cast_shelf = num2str(fjord_matrix.shelf_%d(i_fjord));',which_year));
    eval(sprintf('id_cast_fjord = num2str(fjord_matrix.fjord_%d(i_fjord));',which_year));
    if ~strcmp(id_cast_shelf,'NaN') && ~strcmp(id_cast_fjord,'NaN')
        digitised_id = find([fjords_digitised.id] == fjord_matrix.ID(i_fjord));

        x_fjord = fjords_digitised(digitised_id).x;
        y_fjord = fjords_digitised(digitised_id).y;

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
        m.lat=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lon');

        % fjord obervations for model fitting
        tf = movmean(ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'temperature'),5);
        sf = movmean(ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'salinity'),5);
        zf = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'depth');
        time_fjord_raw = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'time'); % seconds since 1970-01-01T00:00:00Z
        if length(time_fjord_raw) > 1
            time_fjord_raw = time_fjord_raw(1);
        end
        time_fjord = datetime(1970,1,1,0,0,0) + duration(0,0,time_fjord_raw);            % convert to julian date
        time_fjord_in_model = juliandate(time_fjord) - juliandate(time_start);           % convert to the model axis
        [~,i_tsave_fjord_cast] = min(abs(time_fjord_in_model-p.t_save));                 % finds closest point in "output time" to compare later
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
        p.L = fjords_centreline(digitised_id).length.*1e3;
        p.W = fjords_digitised(digitised_id).area/p.L;
        p.Hsill = abs(fjord_matrix.sill_depth(i_fjord));
        p.Hgl   = abs(fjord_matrix.gl_depth(i_fjord));

        if p.Hsill == -999 % if the fjord has no sill
            p.sill=0;        % we ignore the sill
            p.Hsill = p.H; % just so the RPM checks won't fail, and we can get a reasonable ocean profile
        else
            p.sill=1;
        end
        p.H = max(p.Hgl,p.Hsill);
        if p.H == p.Hsill,p.sill=0; end % if we made the fjord as deep as the sill, it essentially has no sill

        if p.H-p.Hsill < 10
            p.Hsill = p.Hsill-(10-p.H-p.Hsill); % workaround so the layer division below the sill works
        end
        if p.H < p.Hsill % the sill cannot be deeper than the fjord itself
            p.Hsill = p.H;
            p.sill = 0;
        end
        if max(zf) > p.H % if the cast inside the fjord goes deeper than what we set as current depth, then the fjord is as deep as the cast
            p.H = ceil(max(zf));
            % p.Hgl = p.H; % since we take the fjord profiles to be as close as possible to the glacier front, this assumption is reasonable for most cases
        end

        % shelf observations for comparison with ocean forcing (possibly bias correction as well?)
        z_shelf = ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'depth');
        t_shelf = movmean(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'temperature'),5);
        s_shelf = movmean(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'salinity'),5);
        time_shelf_raw = ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'time'); % seconds since 1970-01-01T00:00:00Z
        time_shelf = datetime(1970,1,1,0,0,0) + duration(0,0,time_shelf_raw);            % convert to julian date
        time_shelf_in_model = juliandate(time_shelf) - juliandate(time_start);           % convert to the model axis
        [~,i_tsave_shelf_cast] = min(abs(time_shelf_in_model-p.t_save));                 % finds closest point in "output time" to compare later

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
            [~,i_bottom] = max(z_shelf_nonans);
            if max(z_shelf_nonans) < p.H 
                dz_missing = p.H - max(z_shelf_nonans);
                z_shelf = [z_shelf_nonans; max(z_shelf_nonans)+dz_missing];
                t_shelf = interp1(z_shelf_nonans,t_shelf_nonans,z_shelf,'nearest','extrap');
                s_shelf = interp1(z_shelf_nonans,s_shelf_nonans,z_shelf,'nearest','extrap');
            else
                z_shelf = z_shelf_nonans;
                t_shelf = t_shelf_nonans;
                s_shelf = s_shelf_nonans;
            end
    
            % check if we need to extend the profiles upwards (up to 1m depth)
            if min(z_shelf) > 1
                if min(z_shelf) > 50 % however, we do not want to excessively interpolate data
                    fprintf('Fjord %s has no shelf observation until %.2f m. Skipping...\n',m.ID,min(zs))
                    continue
                end
                zs_with_surf = [1; z_shelf];
                c.ts = interp1(z_shelf,t_shelf,zs_with_surf,'nearest','extrap');
                c.ss = interp1(z_shelf,s_shelf,zs_with_surf,'nearest','extrap');
                c.zs = zs_with_surf;
            else
                c.ts = t_shelf;
                c.ss = s_shelf;
                c.zs = z_shelf;
            end
        else
            fprintf('Fjord %d has no shelf data to compare. Skipping...\n',fjord_matrix.ID(i_fjord))
            continue
        end

        % get subglacial discharge (separate file for better readability)
        run set_subglacial_discharge.m
        
        % get the shelf-profile forcing
        ocean = get_ocean_profiles(datasets, x_fjord,y_fjord,-p.Hsill); % gets closest ocean reanalysis cell at least as deep as the sill
        [Ts,Ss] = get_ocean_for_period(ocean,runtime_axis);
        f.Ts = flip(Ts,1);
        f.Ss = flip(Ss,1);
        f.zs = -flip(ocean.depth);
        f.ts = juliandate(runtime_axis) - juliandate(runtime_axis(1));

        % TODO: bias correction?

        % icebergs
        % nu0 controls how much the total iceberg area extends below the surface
        % we do not vary this parameter, as previous sensitivity
        % analyses showed it to have a smaller impact than A0 or M0
        p.nu0=25; 

        % initial conditions
        a.H0 = ones([p.N],1).*p.H/p.N;
        [T_clim, S_clim, ~] = bin_forcings(f, a.H0, t);
        a.T0 = T_clim(:,1);
        a.S0 = S_clim(:,1);


        % create the model structure
        fjord_model(i_fjord).t = t;
        fjord_model(i_fjord).m = m;
        fjord_model(i_fjord).p = p;
        fjord_model(i_fjord).a = a;
        fjord_model(i_fjord).f = f;
        fjord_model(i_fjord).c.tf = tf_nonans;
        fjord_model(i_fjord).c.sf = sf_nonans;
        fjord_model(i_fjord).c.zf = zf_nonans;
        fjord_model(i_fjord).c.i_fjd_cast = i_tsave_fjord_cast;
        fjord_model(i_fjord).c.ts = c.ts;
        fjord_model(i_fjord).c.ss = c.ss;
        fjord_model(i_fjord).c.zs = c.zs;
        fjord_model(i_fjord).c.i_shf_cast = i_tsave_shelf_cast;
        fprintf('Done pre-processing fjord %d.\n',fjord_matrix.ID(i_fjord))
        clear m p a f c
        

    else
        fprintf('No %d data for fjord %d.\n',which_year,fjord_matrix.ID(i_fjord))
    end
end
warning('on','all')
idx = arrayfun(fun,fjord_model);
fjord_model(idx)=[]; % remove the empty elements
n_fjords = length(fjord_model);
fprintf('Inputs processing complete. Total of fjords for %d: %d\n',which_year,n_fjords)

%% Run the model
dims_ensemble = [n_fjords,size(X,1)];
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);
ensemble = struct(ensemble_fields{:});
ensemble(dims_ensemble) = struct("p",[],"t",[],"m",[],"s",[]);
fjords_crashed = {};

run_counter=0;
for i_fjord=1:length(fjord_model)
tic
for i_run=1:n_runs
    run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    
    % Sets the model parameters according to the LHS we created
    for i_param=1:size(X,2)
        cur_fjord.p.(param_names{i_param}) = X(i_run,i_param);
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
            Tfinal=NaN([cur_fjord.p.N,1]);
            Sfinal=NaN(size(Tfinal));
            % Hfinal=NaN(size(Tfinal));
            
            % 10-day avg centered at the target day
            tgt_day = cur_fjord.c.i_fjd_cast;
            Tfinal  = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)),2); 
            Sfinal  = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)),2);
            ensemble(i_fjord,i_run).s.Tfinal = Tfinal; 
            ensemble(i_fjord,i_run).s.Sfinal = Sfinal; 
            
            Tforc = mean(cur_fjord.s.Ts(:,(tgt_day-5:tgt_day+5)),2);
            Sforc = mean(cur_fjord.s.Ss(:,(tgt_day-5:tgt_day+5)),2);
            ensemble(i_fjord,i_run).s.Tforc = Tforc; 
            ensemble(i_fjord,i_run).s.Sforc = Sforc; 

            zf_obs = cur_fjord.c.zf';
            % Tregular=NaN([length(zf_obs),length(cur_fjord.s.t)]);
            % Sregular=NaN(size(Tregular));
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
            fprintf('run %d failed: %s. ',run_counter,ME.message)
            cur_fjord.m.error = ME.message;
            cur_fjord.m.stack = ME.stack;
            fjords_crashed{end+1} = cur_fjord;
        end
        fprintf('\n')
end % for i_run
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
file_out = [outs_path,'rpmMP_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h'];
save(file_out,'-v7.3','ensemble','fjord_model');
if ~isempty(fjords_crashed)
    save([file_out,'_crashed'],'-v7.3','fjords_crashed');
end
disp('Outputs saved.')

%% Postprocessing ensemble
load(file_out);

[res_obs,res_box] = postprocess_ensemble(fjord_model,ensemble);
disp('Postprocessing ensemble done.')
plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(1),[],1);
% exportgraphics(gcf,[figs_path,'profiles_GRLMP_temp_',num2str(which_year),'_n',num2str(n_runs),'.png'],'Resolution',300)