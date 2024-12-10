clearvars; %close all
run setup_paths
iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile
%% Initialise all needed variables


n_runs     = 400;             % number of runs per fjord
input_dt   = 30;              % time step of model input (in days)
dt_in_h    = 2.;              % time step in hours
model_dt   = dt_in_h/24.;     % time step in days for the model (e.g., 2h/24 ~ 0.083 days)
n_years    = 4;               % how many years we want to run
tgt_days   = [935,1010,1150]; % which days of the run we want vertical profiles for: year 3 peak-melt and post-melt, following pre-melt
name_days  = {'peak','post','winter'};

dir_data   = '/Volumes/leg/work/scientific_work_areas/ctd/BASproc';
file_shelf = 'SD041_ctd_034.2db'; % might want to choose a different one later (42 or 34 or 39)
file_fjord = 'SD041_ctd_035.2db'; % might want to choose a different one later (27 or 31 or 35)

%% Define parameter space
param_names = {'C0','wp','K0','A0'};

range_params = {[1e1,1e4],...    % C0 
                [5,750],...     % P0 no crashes with [10,750]
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

%% Get the data

% Subglacial discharge for 2020 (our latest available)
path_runoff = '/Users/mmeb1/Library/CloudStorage/OneDrive/data_common/greenland/runoff/Karlsson2023/';
runoff_file1 = 'SorgenfriGletsjer_CE_D189.csv';
runoff_file2 = 'SorgenfriGletsjer_CE_D193.csv';
runoff1 = readtable([path_runoff,'/',runoff_file1]);
runoff2 = readtable([path_runoff,'/',runoff_file2]);

fjord = load(sprintf('%s/%s.mat',dir_data,file_fjord),'-mat'); % our own measured ice front cast
shelf = load(sprintf('%s/%s.mat',dir_data,file_shelf),'-mat'); % our own measured shelf cast

%% Initialise model

t = 0:model_dt:365*n_years; % 4 years (we want to repeat the yearly data 4 times)

% compile metadata and default parameters
m.ID   = {'R'};
m.name = {'Ryberg'};
m.lon  = -30.545084 ;
m.lat  = 68.232404;
p      = default_parameters();
p.fixedthickness=1;
p.N=60;
p.dt=model_dt;
p.t_save = t(1):1:t(end); % save at daily resolution

% obtain CTD data
ctd.temp  = fjord.potemp1;
ctd.salt  = fjord.salin1;
ctd.depth = fjord.press;

% fjord geometry - REVIEW GL AND SILL DEPTHS
p.L     = 31e3;
p.W     = 3e3;
p.H     = 250;
p.Hsill = 175;
p.Hgl   = 130; % tried (100?), 130,160,190
p.sill  = 1;

% ocean forcing
[~,i_bottom] = max(shelf.press);
if max(shelf.press) < p.H % need to extend the profiles
    dz_missing = p.H - max(shelf.press);
    t_shelf = smoothdata(shelf.potemp1(~isnan(shelf.potemp1) & ~isnan(shelf.salin1)));
    s_shelf = smoothdata(shelf.salin1(~isnan(shelf.potemp1) & ~isnan(shelf.salin1)));
    z_shelf = shelf.press(~isnan(shelf.potemp1) & ~isnan(shelf.salin1));

    zs = [z_shelf; max(z_shelf)+dz_missing];
    Ts = interp1(z_shelf,t_shelf,zs,'nearest','extrap');
    Ss = interp1(z_shelf,s_shelf,zs,'nearest','extrap');
else
    zs = shelf.press(1:i_bottom);
    Ts = smoothdata(shelf.potemp1(1:i_bottom));
    Ss = smoothdata(shelf.salin1(1:i_bottom));
end

zs = flip(zs);
Ts = flip(Ts);
Ss = flip(Ss);

f.zs = -zs;
f.ts = [t(1),t(end)];
f.Ts = repmat(Ts,1,length(f.ts));
f.Ss = repmat(Ss,1,length(f.ts));

% obtain Qsg
total_runoff = runoff1.SurfaceMelt(end-11:end) + runoff1.BasalMelt(end-11:end)...
             + runoff2.SurfaceMelt(end-11:end) + runoff2.BasalMelt(end-11:end);
start_date  = datetime('15-Jan-2020');
cal_runoff = start_date + calmonths(0:n_years*12);
seconds_in_month = eomday(cal_runoff.Year,cal_runoff.Month)*86400; 
time_runoff = juliandate(cal_runoff)' - juliandate(start_date);
f.tsg = time_runoff';
f.Qsg = [repmat(total_runoff,n_years,1); 0]';
f.Qsg = 6.*f.Qsg./seconds_in_month;

% initial conditions
a.H0 = ones([p.N],1).*p.H/p.N;
% [a.T0, a.S0, f.Qsg] = bin_forcings(f, a.H0, t);
[T_clim, S_clim, ~] = bin_forcings(f, a.H0, t);
a.T0 = T_clim(:,1);
a.S0 = S_clim(:,1);

% idealised icebergs
p.icestatic=1; % icebergs do not evolve in MITgcm
p.A0 = 1e9;
p.nu0=25;
f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards

% create the model structure
fjord_model.t = t;
fjord_model.m = m;
fjord_model.p = p;
fjord_model.a = a;
fjord_model.f = f;
fjord_model.c.tf = fjord.potemp1(~isnan(fjord.salin1));
fjord_model.c.sf = fjord.salin1(~isnan(fjord.salin1));
fjord_model.c.zf = fjord.press(~isnan(fjord.salin1))';
fjord_model.c.ts = shelf.potemp1(~isnan(shelf.salin1));
fjord_model.c.ss = shelf.salin1(~isnan(shelf.salin1));
fjord_model.c.zs = shelf.press(~isnan(shelf.salin1))';

% fjord_model = [fjord_model; fjord_model];
%% Running the model
dims_ensemble = [length(fjord_model),size(X,1)];
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);
ensemble = struct(ensemble_fields{:});
ensemble(dims_ensemble) = struct("p",[],"t",[],"m",[],"s",[]);
fjords_crashed = {};

run_counter=0;
tic
for i_fjord=1:1%length(fjord_model)
for i_run=1:n_runs
    run_counter = run_counter+1;
    cur_fjord = fjord_model(1);
    
    for i_param=1:size(X,2)
        cur_fjord.p.(param_names{i_param}) = X(i_run,i_param);
    end
    cur_fjord.a.I0 = cur_fjord.p.A0*iceberg_fun(cur_fjord.p.nu0, abs(cur_fjord.p.Hgl), -cumsum(cur_fjord.a.H0)+cur_fjord.a.H0/2);
    
    % cur_fjord.p.plot_runtime=1;
        try
            cur_fjord.s = run_model(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
            fprintf('run %d complete. ',run_counter)
    
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
end % for i_fjord
toc
fprintf('\n')
file_out = [outs_path,'rpm_Ryberg_n',num2str(n_runs),'_',num2str(cur_fjord.p.N),'layers_dt',num2str(dt_in_h),'h_flipped_qsg6x'];
save(file_out,'-v7.3','ensemble','fjord_model');
if ~isempty(fjords_crashed)
    save([file_out,'_crashed'],'-v7.3','fjords_crashed');
end
disp('Outputs saved.')
% file_in = [outs_path,'rpm_Ryberg_n400_60layers_dt3h.mat'];
% load(file_in)
%% Post-processing
% load(file_out);
% TODO: check RMSE
[res_obs,res_box] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
disp('Postprocessing ensemble done.')

%% Plotting

plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2);
% exportgraphics(gcf,[figs_path,'profiles_temp_',num2str(2015+i_yr),'_n',num2str(n_runs),'.png'],'Resolution',300)

% plot_best_params(fjord_model,ensemble,res_box,param_names,range_params,2);
% exportgraphics(gcf,[figs_path,'best_parameters_',num2str(2015+i_yr),'_n',num2str(n_runs),'.png'],'Resolution',300)
