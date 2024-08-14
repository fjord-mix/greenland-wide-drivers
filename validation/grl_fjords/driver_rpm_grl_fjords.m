%% Driver file for simulating the compiled fjords
clearvars
warning('off','all')
run setup_paths % Configuring paths
letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};

%% Initialise all needed variables

which_year = 2016; % {2016,2017,2018,2019,2020}
time_start = datetime(which_year,01,01);
time_end   = datetime(which_year,12,31);

n_runs     = 400;          % number of runs per fjord
input_dt   = 30;          % time step of model input (in days)
dt_in_h    = 2.;          % time step in hours
model_dt   = dt_in_h/24.; % time step in days for the model (e.g., 2h/24 ~ 0.083 days)
n_years    = 4;           % how many years we want to run
% tgt_days   = [935,1010];  % which days of the run we want vertical profiles for
tgt_days   = [1280,1400];  % which days of the run we want vertical profiles for
name_days  = {'peak','post'};

t = 0:model_dt:365*n_years; % we want to repeat the yearly data <n_years> times

fun = @(s) all(structfun(@isempty,s));              % tiny function to get rid of empty entries in array
iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile

%% Define parameter space
param_names = {'C0','wp','K0','A0'};

range_params = {[1e2,1e4],...    % C0 
                [10,700],...     % P0 no crashes with [10,400]
                [1e-4,1e-3],...  % K0 or do we stick to [1e4,1e-3]?
                [0,3e8]};        % A0

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

%% Compiling data 

% Read compilation of all fjords around Greenland (requires data-compilation repository)
% only needed if not using fjord geometries from Cowton et al.
% run compile_process_fjords 
% fjord_ids = [4,9,17,20,22,23,24,25,28,29,30,31,24,37]; % These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
file_fjords_compiled = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_gl_sill_depths.xlsx'];
folder_ctd_casts     = [data_path,'/greenland/obs/OMG_all_casts'];
file_lengths = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_centreline.shp'];
file_fjords = [data_path,'/greenland/FjordMIX/fjords_digitisation/fjords_grl.shp'];


fjords_digitised  = shaperead(file_fjords);
fjords_centreline = shaperead(file_lengths);

fjord_matrix = readtable(file_fjords_compiled);
fjord_matrix(fjord_matrix.gl_depth < 50,:) = [];
fjord_matrix(isnan(fjord_matrix.qsg_id1),:) = [];

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(height(fjord_matrix)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);

for i_fjord=1:height(fjord_matrix)
    eval(sprintf('id_cast_shelf = num2str(fjord_matrix.shelf_%d(i_fjord));',which_year));
    eval(sprintf('id_cast_fjord = num2str(fjord_matrix.fjord_%d(i_fjord));',which_year));
    if ~strcmp(id_cast_shelf,'NaN') && ~strcmp(id_cast_fjord,'NaN')
        digitised_id = find([fjords_digitised.id] == fjord_matrix.ID(i_fjord));

        % default parameters
        p=default_parameters();
        p.fixedthickness=1;
        p.N=60;
        p.dt=model_dt;
        p.t_save = t(1):1:t(end); % save at daily resolution

        % get the metadata
        m.ID = num2str(fjord_matrix.ID(i_fjord));
        m.name = fjord_matrix.name(i_fjord);
        
        % get geometry
        p.L = fjords_centreline(digitised_id).length.*1e3;
        p.W = fjords_digitised(digitised_id).area/p.L;
        p.Hsill = fjord_matrix.sill_depth(i_fjord);
        p.Hgl   = fjord_matrix.gl_depth(i_fjord);
        p.H = abs(p.Hgl);
        if p.Hsill > p.Hgl || p.Hsill == -999 % if the GL sits above the sill or if the fjord has no sill
            p.sill=0;        % we ignore the sill
            p.Hsill = p.Hgl; % just so the RPM initial checks won't fail
        else
            p.sill=1;
        end

        % get subglacial discharge (separate file for better readability
        run set_subglacial_discharge.m

        % find the ocean forcing and fjord profile
        omg_data_shelf = dir([folder_ctd_casts,'/*',id_cast_shelf,'*.nc']);
        omg_data_fjord = dir([folder_ctd_casts,'/*',id_cast_fjord,'*.nc']);
        
        % add rough coordinates
        m.lon=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lon');
        m.lat=ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'lon');
        
        % get the shelf-profile forcing
        z_shelf = smoothdata(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'depth'));
        t_shelf = smoothdata(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'temperature'));
        s_shelf = smoothdata(ncread([omg_data_shelf.folder,'/',omg_data_shelf.name],'salinity'));
        is_downcast = [1; diff(z_shelf) > 0]; % removes upcast, if it exists
        
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
                zs = [z_shelf; max(z_shelf_nonans)+dz_missing];
                Ts = interp1(z_shelf_nonans,t_shelf_nonans,zs,'nearest','extrap');
                Ss = interp1(z_shelf_nonans,s_shelf_nonans,zs,'nearest','extrap');
            else
                zs = z_shelf_nonans;
                Ts = t_shelf_nonans;
                Ss = s_shelf_nonans;
            end

            % check if we need to extend the profiles upwards
            if min(zs) > 5
                zs_with_surf = [5; zs];
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
            [T_clim, S_clim, ~] = bin_forcings(f, a.H0, t);
            a.T0 = T_clim(:,1);
            a.S0 = S_clim(:,1);
    
            % obervations for reference and model fitting
            tf = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'temperature');
            sf = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'salinity');
            zf = ncread([omg_data_fjord.folder,'/',omg_data_fjord.name],'depth');
            tf_nonans = tf(~isnan(tf) & ~isnan(sf));
            sf_nonans = sf(~isnan(tf) & ~isnan(sf));
            zf_nonans = zf(~isnan(tf) & ~isnan(sf));
    
    
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
warning('on','all')
idx = arrayfun(fun,fjord_model);
fjord_model(idx)=[]; % remove the empty elements
n_fjords = length(fjord_model);
fprintf('Inputs processing complete. Total of fjords for %d: %d\n',which_year,n_fjords)

%% Running the model
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
            Tfinal=NaN([cur_fjord.p.N,length(tgt_days)]);
            Sfinal=NaN(size(Tfinal));
            Hfinal=NaN(size(Tfinal));
            for i_day=1:length(tgt_days)
                % 10-day avg centered at the target day
                tgt_day = tgt_days(i_day);
                Tfinal(:,i_day) = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)),2); 
                Sfinal(:,i_day) = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)),2);
            end
            ensemble(i_fjord,i_run).s.Tfinal = Tfinal;
            ensemble(i_fjord,i_run).s.Sfinal = Sfinal;

            zf_obs = cur_fjord.c.zf';
            Tregular=NaN([length(zf_obs),length(cur_fjord.s.t)]);
            Sregular=NaN(size(Tregular));
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
file_out = [outs_path,'rpm_GRL_fjords_n',num2str(n_runs),'_',num2str(which_year),'_',num2str(fjord_model(1).p.N),'layers_dt',num2str(dt_in_h),'h'];
save(file_out,'-v7.3','ensemble','fjord_model');
if ~isempty(fjords_crashed)
    save([file_out,'_crashed'],'-v7.3','fjords_crashed');
end
disp('Outputs saved.')

% run postprocess_plot_ensembles.m
%% Sanity-check post processing and plot
% load(file_out);

% [res_obs,res_box] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
% disp('Postprocessing ensemble done.')
% plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2);
