%% Driver file for simulating the 10 NW fjords

run setup_paths % Configuring paths

%% Initialise all needed variables

n_runs     = 400;              % number of runs per fjord
input_dt   = 30;              % time step of model input (in days)
dt_in_h    = 3.;              % time step in hours
model_dt   = dt_in_h/24.;     % time step in days for the model (e.g., 2h/24 ~ 0.083 days)
n_years    = 4;               % how many years we want to run
tgt_days   = [935,1010,1150]; % which days of the run we want vertical profiles for: year 3 peak-melt and post-melt, following pre-melt
name_days  = {'peak','post','winter'};

i_yr       = 4; % {1,2,3,4} for {2016,2017,2018,2019}
time_start = datetime(2015+i_yr,01,01);
time_end   = datetime(2015+i_yr,12,31);

letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};

fjord_names = {'Qeqertaarsuusarsuup','Nuussuup','Tuttulikassaap','Nunatakassaap','Illulip',...
               'Kakiffaat','Naajarsuit','Sermeq',...
               'Salliarsutsip','Umiammakku','Kangilliup','Kangerlussuup','SermeqSilarleq','SermeqKujalleq'} ;

fun = @(s) all(structfun(@isempty,s));              % tiny function to get rid of empty entries in array
iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile

%% Define parameter space
param_names = {'C0','wp','K0','A0'};

range_params = {[1e2,5e4],...    % C0 
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
% These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
% fjord_ids = [4,9,17,20,22,23,24,25,28,29,30,31,24,37];

run load_cowton2023_data % Reading the data from Cowton et al. (2023)

tgt_period = time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up
disp('Data compiling done.')

%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords
    if ~isnan(meta_table.Length_km_(i_fjord)) % we only want fjords for which we have all geometry data

        % compile metadata and default parameters
        m.ID     = meta_table.Var1(i_fjord);
        m.name   = meta_table.GreenlandicName(i_fjord);
        m.lon    = meta_table.Longitude__E_(i_fjord);
        m.lat    = meta_table.Latitude__N_(i_fjord);
        p        = default_parameters();
        p.N      = 60;
        t        = 0:model_dt:365*n_years; % time axis for the run: 4 years (we want to repeat the yearly data 4 times)
        p.t_save = t(1):1:t(end);          % save at daily resolution
        p.dt     = model_dt;
        
        % obtain CTD data
        ctd.temp  = ctd_data.Tdata{i_yr,i_fjord};
        ctd.salt  = ctd_data.Sdata{i_yr,i_fjord};
        ctd.depth = ctd_data.zdata{i_yr,i_fjord};
        
        % fjord geometry - from Cowton et al. (2023)
        p.L         = 1e3.*meta_table.Length_km_(i_fjord);
        p.W         = 1e3.*meta_table.Width_km_(i_fjord);
        p.Hsill     = abs(meta_table.SillDepth_m_(i_fjord));
        p.Hgl       = abs(meta_table.GroundingLineDepth_m_(i_fjord));

        % fjord geometry - our data compilation
        % zgs = NaN(size(fjords_compilation(fjord_ids(i_fjord)).glaciers));
        % for i_glacier=1:length(fjords_compilation(fjord_ids(i_fjord)).glaciers)
        %     zgs(i_glacier) = fjords_compilation(fjord_ids(i_fjord)).glaciers(i_glacier).gldepth;
        % end
        % 
        % p.L         = fjords_compilation(fjord_ids(i_fjord)).length;
        % p.W         = fjords_compilation(fjord_ids(i_fjord)).width;
        % p.Hsill = fjords_compilation(fjord_ids(i_fjord)).silldepth;
        % p.Hgl       = min(zgs);

        % fjord geometry - fjord depth based on GL/sill depths
        if p.Hsill > p.Hgl % if the GL sits above the sill
            p.sill=0; % there is no sill
        else
            p.sill=1;
        end
        p.H = abs(p.Hgl); % the fjord is always as deep as the GL
        

        % ocean forcing
        ts = ctd.temp.s;
        ss = ctd.salt.s;
        zs = -ctd.depth.s;
        if max(abs(zs)) < p.H % if our profile does not extend to the bottom of the fjord
            h_missing = p.H - max(abs(zs));
            depth_range_missing = max(abs(zs))+1:1:p.H;
            new_zs = [zs -depth_range_missing];
            ts = interp1(zs,ts,new_zs,'nearest','extrap');
            ss = interp1(zs,ss,new_zs,'nearest','extrap');
            zs = new_zs;
            clear new_zs h_missing depth_range_missing % tidying up
        end
        f.zs = zs';
        f.Ts = repmat(ts,length(t),1)';
        f.Ss = repmat(ss,length(t),1)';


        % initial conditions
        H_clim = ones([p.N,1]).*p.H/p.N;
        [T_clim,S_clim] = bin_ocean_profiles(mean(f.Ts,2),mean(f.Ss,2),f.zs,H_clim);
        a.H0 = H_clim;
        a.T0 = T_clim;
        a.S0 = S_clim;

        % flipping the ocean forcing like the model wants it
        % f.zs = flip(f.zs);
        % f.Ts = flip(f.Ts,1); 
        % f.Ss = flip(f.Ss,1);
        
        % obtain Qsg
        % interpolate from daily Qsg data to the model time axes
        taxis_qsg_julian = juliandate(taxis_qsg);
        taxis_qsg_num = taxis_qsg_julian-taxis_qsg_julian(1);
        taxis_qsg_num_dtmodel = taxis_qsg_num(1):t(2)-t(1):taxis_qsg_num(end);
        qsg_dtmodel = interp1(taxis_qsg_num,qsg_all(:,i_fjord),taxis_qsg_num_dtmodel);
        f.Qsg = repmat(qsg_dtmodel,1,11);
        f.Qsg = f.Qsg(1:length(t));
        clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up

        % idealised icebergs
        p.nu0 = 25; % how fast iceberg area decays with depth
        % p.icestatic=1; % icebergs do not evolve
        % f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
        % a.I0 = p.A0*iceberg_fun(p.nu0, abs(p.Hgl), -cumsum(a.H0)+a.H0/2);

        fjord_model(i_fjord).t = t;
        fjord_model(i_fjord).m = m;
        fjord_model(i_fjord).p = p;
        fjord_model(i_fjord).a = a;
        fjord_model(i_fjord).f = f;
        fjord_model(i_fjord).c.tf = ctd.temp.f;
        fjord_model(i_fjord).c.sf = ctd.salt.f;
        fjord_model(i_fjord).c.zf = ctd.depth.f;
        fjord_model(i_fjord).c.ts = ctd.temp.s;
        fjord_model(i_fjord).c.ss = ctd.salt.s;
        fjord_model(i_fjord).c.zs = ctd.depth.s;

        % Sanity-check plot to make sure the boxes and binned T and S make sense
        % % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
        % figure('Position',[100 100 800 500]);
        % ints=[0,-cumsum(a.H0)];
        % subplot(1,3,1); hold on; box on;
        % plot(f.Ts(:,1),f.zs);
        % scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
        % for j=1:size(ints,2)
        %     yline(ints(j));
        % end
        % yline(p.Hgl,'--k','linewidth',0.5)
        % yline(p.Hsill,'--k','linewidth',0.5)
        % ylim([-p.H 0])
        % subplot(1,3,2); hold on; box on;
        % plot(f.Ss(:,1),f.zs);
        % scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
        % for j=1:size(ints,2)
        %     yline(ints(j));
        % end
        % yline(p.Hgl,'--k','linewidth',0.5)
        % yline(p.Hsill,'--k','linewidth',0.5)
        % ylim([-p.H 0])
        % subplot(1,3,3); hold on; box on; grid on;
        % plot(t,f.Qsg);
        % plot(t,f.D);
        % legend('Subglacial discharge','Solid-ice discharge');
        % exportgraphics(gcf,[figs_path,'forcing_summary_c13_fjord_',letters{i_fjord},'_',num2str(2015+i_yr),'.png'],'Resolution',300)
        % 
        clear m p a f % bit of tidying up
    end
end

idx = arrayfun(fun,fjord_model);
fjord_model(idx)=[]; % remove the empty elements
disp('Model inputs processing done.')

%% Running the model
dims_ensemble = [length(fjord_model),size(X,1)];
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);
ensemble = struct(ensemble_fields{:});
ensemble(dims_ensemble) = struct("p",[],"t",[],"m",[],"s",[]);
fjords_crashed = {};

% we still need to add additional 'for' statements and manually add extra
% indices within the for loops for any new parameter to be added. 
% But at least now it should be easier to spot where to add
% them and what might be missing... TODO: replace these for loops for LHS
run_counter=0;
for i_fjord=1:length(fjord_model)
tic
for i_run=1:n_runs
    run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    
    for i_param=1:size(X,2)
        cur_fjord.p.(param_names{i_param}) = X(i_run,i_param);
    end
    cur_fjord.a.I0 = cur_fjord.p.A0*iceberg_fun(cur_fjord.p.nu0, abs(cur_fjord.p.Hgl), -cumsum(cur_fjord.a.H0)+cur_fjord.a.H0/2);
    
    % cur_fjord.p.plot_runtime=1;
        try
            cur_fjord.s = run_model(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
            fprintf('run %d complete. ',run_counter)
    
            ensemble(i_fjord,i_run).p = cur_fjord.p;
            ensemble(i_fjord,i_run).t = cur_fjord.t;
            ensemble(i_fjord,i_run).m = cur_fjord.m;
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
            for k=1:length(cur_fjord.s.t)
                Tregular(:,k) = interp1(z_box,cur_fjord.s.T(:,k),zf_obs,'nearest','extrap');
                Sregular(:,k) = interp1(z_box,cur_fjord.s.S(:,k),zf_obs,'nearest','extrap');
            end
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
file_out = [outs_path,'rpm_NWfjords_n',num2str(n_runs),'_',num2str(2015+i_yr),'_',num2str(cur_fjord.p.N),'layers_dt',num2str(dt_in_h),'h_v2'];
save(file_out,'-v7.3','ensemble','fjord_model');
if ~isempty(fjords_crashed)
    save([file_out,'_crashed'],'-v7.3','fjords_crashed');
end
disp('Outputs saved.')

%% Post-processing

if exist('res_obs',"var"),       clear res_obs; end
res_obs(size(fjord_model)) = struct("tf",[],"sf",[],"zf",[]);
if exist('res_box',"var"),       clear res_box; end
res_box(size(fjord_model)) = struct("tf",[],"sf",[],"zf",[],"id",[],"name",[]);
for i_fjord=1:size(ensemble,1)
    zf_obs = fjord_model(i_fjord).c.zf;
    tf_obs = fjord_model(i_fjord).c.tf;
    sf_obs = fjord_model(i_fjord).c.sf;
    n_completed = 0; % completed runs for that fjord

    tf_box_comp = NaN([length(zf_obs),n_runs,length(tgt_days)]);
    sf_box_comp = NaN(size(tf_box_comp));
    rmse_tf = NaN([n_runs,length(tgt_days)]);
    rmse_sf = NaN([n_runs,length(tgt_days)]);
    tupper_box_comp = NaN([length(cur_fjord.s.t),n_runs]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    supper_box_comp = NaN(size(tupper_box_comp));
    sinter_box_comp = NaN(size(tupper_box_comp));
    slower_box_comp = NaN(size(tupper_box_comp));

    for i_run=1:n_runs
        if ~isempty(ensemble(i_fjord,i_run).s) & ~isnan(ensemble(i_fjord,i_run).s.Tfinal(1))

            % we get a profile for each of our target days - if something
            % goes wrong here it's because of the sign/orientation of depth
            % profiles
            for i_day=1:length(tgt_days) 
                tf_box = ensemble(i_fjord,i_run).s.Tfinal(:,i_day); 
                sf_box = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                
                tf_box_comp(:,i_run,i_day) = interp1(ensemble(i_fjord,i_run).s.z,tf_box,zf_obs,'nearest','extrap');
                sf_box_comp(:,i_run,i_day) = interp1(ensemble(i_fjord,i_run).s.z,sf_box,zf_obs,'nearest','extrap');
            end

            tupper_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tupper;
            tinter_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tinter;
            tlower_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tlower;

            supper_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Supper;
            sinter_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Sinter;
            slower_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Slower;

            min_depth_rmse = 0;
            depths_rmse = zf_obs > min_depth_rmse;
            for i_day = 1:length(tgt_days)
                rmse_tf(i_run,i_day) = rmse(tf_box_comp(depths_rmse,i_run,i_day),tf_obs(depths_rmse)','omitnan')./mean(tf_obs(depths_rmse),'omitnan');
                rmse_sf(i_run,i_day) = rmse(sf_box_comp(depths_rmse,i_run,i_day),sf_obs(depths_rmse)','omitnan')./mean(sf_obs(depths_rmse),'omitnan');
            end
            n_completed = n_completed+1;
        end
    end % i_run
    res_obs(i_fjord).zf = zf_obs;
    res_obs(i_fjord).tf = tf_obs;
    res_obs(i_fjord).sf = sf_obs;
    
    res_obs(i_fjord).zs = fjord_model(i_fjord).c.zs;
    res_obs(i_fjord).ts = fjord_model(i_fjord).c.ts;
    res_obs(i_fjord).ss = fjord_model(i_fjord).c.ss;

    res_box(i_fjord).t  = cur_fjord.s.t;
    res_box(i_fjord).zf = z_box;
    res_box(i_fjord).tf = mean(tf_box_comp,2,'omitnan');
    res_box(i_fjord).tfmin = min(tf_box_comp,[],2,'omitnan');
    res_box(i_fjord).tfmax = max(tf_box_comp,[],2,'omitnan');
    res_box(i_fjord).sf = mean(sf_box_comp,2,'omitnan');
    res_box(i_fjord).sfmin = min(sf_box_comp,[],2,'omitnan');
    res_box(i_fjord).sfmax = max(sf_box_comp,[],2,'omitnan');

    res_box(i_fjord).Tupper = mean(tupper_box_comp,2,'omitnan');
    res_box(i_fjord).Tupper_min = min(tupper_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tupper_max = max(tupper_box_comp,[],2,'omitnan');

    res_box(i_fjord).Tinter = mean(tinter_box_comp,2,'omitnan');
    res_box(i_fjord).Tinter_min = min(tinter_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tinter_max = max(tinter_box_comp,[],2,'omitnan');

    res_box(i_fjord).Tlower = mean(tlower_box_comp,2,'omitnan');
    res_box(i_fjord).Tlower_min = min(tlower_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tlower_max = max(tlower_box_comp,[],2,'omitnan');

    res_box(i_fjord).Supper = mean(supper_box_comp,2,'omitnan');
    res_box(i_fjord).Sinter = mean(sinter_box_comp,2,'omitnan');
    res_box(i_fjord).Slower = mean(slower_box_comp,2,'omitnan');

    res_box(i_fjord).ensemble_tf = tf_box_comp;
    res_box(i_fjord).ensemble_sf = sf_box_comp;
    res_box(i_fjord).rmse_tf = rmse_tf;
    res_box(i_fjord).rmse_sf = rmse_sf;

    res_box(i_fjord).id = fjord_model(i_fjord).m.ID{1};
    res_box(i_fjord).name = fjord_model(i_fjord).m.name{1};
    res_box(i_fjord).n = n_completed./n_runs * 100;
end
disp('Postprocessing done.')

%% Plotting

plot_ensemble_profiles(fjord_model,ensemble,res_box,res_obs,n_runs,param_names,tgt_days(2),name_days,2,mitgcm);
% exportgraphics(gcf,[figs_path,'profiles_temp_',num2str(2015+i_yr),'_n',num2str(n_combinations),'.png'],'Resolution',300)

plot_best_params(fjord_model,ensemble,res_box,param_names,range_params,2);
% exportgraphics(gcf,[figs_path,'best_parameters_',num2str(2015+i_yr),'_n',num2str(n_combinations),'.png'],'Resolution',300)

% plot_sensitivity_to_param(res_box,res_obs,fjord_model,param_names,range_params,2);
plot_sensitivity_profiles(X,ensemble,res_box,param_names,range_params,2)