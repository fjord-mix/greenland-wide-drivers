%% Driver file for simulating the fjords simulated by Tom using MITgcm

% quick debug:
% load('/Users/mmeb1/OneDrive - University of St Andrews/data_common/greenland/FjordMIX/boxmodel/cowton2023/boxmodel_example_bad_H_negT_lim')
% [cur_fjord.s,cur_fjord.f] = zmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);

%% Initialise all needed paths and variables
run setup_paths

i_yr = 4; % we want only 2019 data for now, because we only have geometry data for glaciers with data for 2019
time_start = datetime(2019,01,01);
time_end   = datetime(2020,02,15);
input_dt   = 30;
dt_in_h    = 3.;
model_dt   = dt_in_h/24.; % time step in days (2h/24 ~ 0.083 days)
tgt_day = 300; % which day of the 400-day run we want to compare

flag_ice          = 'ideal'; %{'MIT'|'obs'|'ideal'|'no'} for choosing which iceberg treatment (MIT-derived, observation-derived, no icebergs)
fjord_ids_all_shp = [14,17,23,24,25,30,31,34,37];
fjord_ids = [14,23,24,34];

iceberg_fun = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of idealised iceberg depth profile

%% Compiling all "external" data used (Greenland fjord compilation and Cowton et al., 2023)
run compile_process_fjords % requires data-compilation repository

run load_cowton2023_data

%% Selecting which fjords from which compilation we want, and the time period for discharge data

% These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
% useful function to find them:
% plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off';

% will look for matching letters in the table to find which fjords we want
ids_cowton_fjords = zeros(size(fjord_ids_mitgcm));
for i=1:height(meta_table)
    for i1=1:length(fjord_ids_mitgcm)
        if strcmp(meta_table.Var1(i),fjord_ids_mitgcm{i1}), ids_cowton_fjords(i1)=i;end
    end
end

tgt_period = 150:250; %time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up

qsg_all    = qsg_all(:,ids_cowton_fjords);
meta_table = meta_table(ids_cowton_fjords,:);

disp('Fjord selection for ensembles done.')

%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords

    i_fjord_ctd = ids_cowton_fjords(i_fjord);
    % time axis for the run
    t = 0:model_dt:400; % 400 days (as in the MITgcm simulations)

    % compile metadata and default parameters
    m.ID   = meta_table.Var1(i_fjord);
    m.name = meta_table.GreenlandicName(i_fjord);
    m.lon  = meta_table.Longitude__E_(i_fjord);
    m.lat  = meta_table.Latitude__N_(i_fjord);
    [p,~]  = get_model_default_parameters();
    p.fixedthickness=1;
    p.N=60;
    p.dt=model_dt;
    p.t_save=t(1):1:t(end); % save output at daily resolution
    p.M0 = 5e-7;
    
    % obtain CTD data
    ctd.temp  = ctd_data.Tdata{i_yr,i_fjord_ctd};
    ctd.salt  = ctd_data.Sdata{i_yr,i_fjord_ctd};
    ctd.depth = ctd_data.zdata{i_yr,i_fjord_ctd};
    

    % fjord geometry
    p.L         = 1e3.*meta_table.Length_km_(i_fjord);
    p.W         = 1e3.*meta_table.Width_km_(i_fjord);
    p.silldepth = -1.* meta_table.SillDepth_m_(i_fjord);
    p.zgl       = -1.* meta_table.GroundingLineDepth_m_(i_fjord);
    if p.silldepth < p.zgl % if the GL sits above the sill
        p.H = abs(p.silldepth); % we consider the fjord depth to be just Hmin below the sill, 
        p.sill=0;               % so we can avoid fiddling with the number of layers
    else
        p.H = abs(p.zgl); % if the GL is below the sill, then we treat it as the same as the fjord depth
    end

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
    f.zs = zs;
    f.Ts = repmat(ts,length(t),1);
    f.Ss = repmat(ss,length(t),1);


    % initial conditions
    H_clim = ones([p.N],1).*p.H/p.N;
    [T_clim,S_clim] = bin_ocean_profiles(mean(f.Ts,1),mean(f.Ss,1),f.zs,H_clim);
    a.H0 = H_clim;
    a.T0 = T_clim;
    a.S0 = S_clim;

    % flipping the ocean forcing like the model wants it
    f.zs = flip(f.zs);
    f.Ts = flip(f.Ts,2)';
    f.Ss = flip(f.Ss,2)';
    
    % obtain Qsg
    % interpolate from daily Qsg data to the model time axes
    taxis_qsg_julian = juliandate(taxis_qsg);
    taxis_qsg_num = taxis_qsg_julian-taxis_qsg_julian(1);
    taxis_qsg_num_dtmodel = taxis_qsg_num(1):t(2)-t(1):taxis_qsg_num(end);
    qsg_dtmodel = interp1(taxis_qsg_num,qsg_all(:,i_fjord),taxis_qsg_num_dtmodel);
    f.Qsg = [zeros([1,200/model_dt]), qsg_dtmodel, zeros([1,100/model_dt])];
    clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up
    % no meltwater
    % f.Qsg=zeros(size(t));
    % p.P0=0;

    % no icebergs
    if strcmp(flag_ice,'no')
        p.M0=0; 
        f.zi = -800:10:0;
        % f.zi = f.zi';
        f.D  = zeros(1,length(t));
        f.xi = (p.nu0/abs(p.zgl))*exp(p.nu0*f.zi'/abs(p.zgl))/(1-exp(-p.nu0));
        a.I0 = 0*f.zi;
    end

    % obs-based
    if strcmp(flag_ice,'obs')
        runtime_axis = datasets.opts.time_start:model_dt:datasets.opts.time_end;
        [f.zi,ice_d,f.xi,a.I0] = get_total_iceberg_discharge(datasets,fjords_compilation(fjord_ids(i_fjord)).glaciers,runtime_axis',p);
        f.D = ice_d(1:length(t)); % we need to interp/extrapolate the 1-year discharge for 400 days
    end

    % MITgcm-based
    if strcmp(flag_ice,'MIT')
        p.icestatic=1; % icebergs do not evolve in MITgcm
        f.xi = flip(mitgcm(i_fjord).Iprofile./trapz(mitgcm(i).z,mitgcm(i).Iprofile)); % get avg. iceberg concentration in each layer and make its integral equal 1
        f.zi = -flip(mitgcm(i_fjord).z);
        f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
        D0 = sum(mitgcm(i_fjord).Ivolume(:),'omitnan');
        z0 = cumsum(a.H0);
        a.I0 = interp1(mitgcm(i_fjord).z,mitgcm(i_fjord).Iarea,z0,'pchip')'; %(D0/(p.W*p.L)).*f.xi;
        fjord_model(i_fjord).a.D0 = D0;
    end

    if strcmp(flag_ice,'ideal')
        f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
    end

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

    % % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
    % figure('Position',[100 100 800 500]);
    % ints=[0,-cumsum(a.H0)];
    % subplot(1,3,1); hold on; box on;
    % plot(f.Ts(:,1),f.zs);
    % scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,2); hold on; box on;
    % plot(f.Ss(:,1),f.zs);
    % scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,3); hold on; box on; grid on;
    % plot(t,f.Qsg);
    % yyaxis right
    % yline(D0);
    % legend('Subglacial discharge','Solid-ice discharge');
    % %
    % exportgraphics(gcf,[figs_path,'forcing_summary_MITice_fjord',fjord_ids_mitgcm{i_fjord},'.png'],'Resolution',300)IP
    % 

    clear m p a f % bit of tidying up
end
clear ctd ctd_data D0
disp('Model inputs processing done.')
%% Create parameter space

n_fjord_runs = length(fjord_model);

param_names = {'C0','P0','K0','A0'};
range_params = {[1e2,1e5],...    % C0 
                [01,40],...      % P0 no crashes with [1,40]
                [1e-4,1e-3],...  % K0
                [0,3e8]};        % A0

dims_ensemble = n_fjord_runs;
n_combinations = 1;%n_fjord_runs;
for i=1:length(range_params)
    dims_ensemble = [dims_ensemble length(range_params{i})];
    n_combinations = n_combinations.*length(range_params{i});
end

disp('Parameter space for ensemble created.')

%% Run the model
% Creating empty struct with desired fields of the size according to how
% many variables we have
ensemble_fields = {'p','t','m','s'}; % 'a','f','c','s'};
ensemble_fields{2,1} = cell(dims_ensemble);
ensemble = struct(ensemble_fields{:});
fjords_crashed = {};
% ensemble(dims_ensemble) = struct("p",[],"t",[],"m",[],"s",[]);

% we still need to add additional 'for' statements and manually add extra
% indices within the for loops for any new parameter to be added. 
% But at least now it should be easier to spot where to add
% them and what might be missing... TODO: replace these for loops for LHS
% i_fjord=3; i1=1;i2=3;i3=1;i4=3;i5=1;
% cur_fjord = fjord_model(i_fjord);
% cur_fjord.p.(param_names{1}) = range_params{1}(i1);
% cur_fjord.p.(param_names{2}) = range_params{2}(i2);
% cur_fjord.p.(param_names{3}) = range_params{3}(i3);
% cur_fjord.p.(param_names{4}) = range_params{4}(i4);
% cur_fjord.p.(param_names{5}) = range_params{5}(i5);
% [cur_fjord.s,cur_fjord.f] = zmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);

run_counter=0;
for i_fjord=1:n_fjord_runs
tic
for i1=1:dims_ensemble(2)
for i2=1:dims_ensemble(3)
for i3=1:dims_ensemble(4)
for i4=1:dims_ensemble(5)
    run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    
    cur_fjord.p.(param_names{1}) = range_params{1}(i1);
    cur_fjord.p.(param_names{2}) = range_params{2}(i2);
    cur_fjord.p.(param_names{3}) = range_params{3}(i3);
    cur_fjord.p.(param_names{4}) = range_params{4}(i4);
    if strcmp(flag_ice,'no')
        cur_fjord.p.M0=0; 
    end
    if strcmp(flag_ice,'ideal')
        cur_fjord.p.icestatic=1;
        cur_fjord.f.D = zeros(size(cur_fjord.t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
        cur_fjord.a.I0 = cur_fjord.p.A0*cur_fjord.p.if(cur_fjord.p.nu0, abs(cur_fjord.p.zgl), -cumsum(cur_fjord.a.H0)+cur_fjord.a.H0/2);
    end
    % cur_fjord.p.plot_runtime=1;
    try
        cur_fjord.s = zmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
        fprintf('run %d complete. ',run_counter)

        ensemble(i_fjord,i1,i2,i3,i4).p = cur_fjord.p;
        ensemble(i_fjord,i1,i2,i3,i4).t = cur_fjord.s.t;
        ensemble(i_fjord,i1,i2,i3,i4).m = cur_fjord.m;

        % 10-day avg centered at the target day (since MITgcm outputs a 10-day average)
        ensemble(i_fjord,i1,i2,i3,i4).s.Tfinal = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)),2); 
        ensemble(i_fjord,i1,i2,i3,i4).s.Sfinal = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)),2);
        ensemble(i_fjord,i1,i2,i3,i4).s.Hfinal = mean(cur_fjord.s.H(:,(tgt_day-5:tgt_day+5)),2);

        % "resampling" the box model results to the same z and t axes as MITgcm 
        % (regular depth levels and a smaller time axis for better memory and disk space management)
        Tregular=NaN([length(mitgcm(i_fjord).z),length(mitgcm(i_fjord).t)]);
        Sregular=NaN([length(mitgcm(i_fjord).z),length(mitgcm(i_fjord).t)]);
        temp_in_mitgcm_time = interp1(cur_fjord.s.t,cur_fjord.s.T',mitgcm(i_fjord).t,'linear','extrap')';
        salt_in_mitgcm_time = interp1(cur_fjord.s.t,cur_fjord.s.S',mitgcm(i_fjord).t,'linear','extrap')';
        h_in_mitgcm_time    = interp1(cur_fjord.s.t,cur_fjord.s.H',mitgcm(i_fjord).t,'linear','extrap')';
        
        % we get MITgcm averaged T for each of the boxmodel boxes
        tf_mitgcm_boxed = NaN([cur_fjord.p.N,length(mitgcm(i_fjord).t)]);
        sf_mitgcm_boxed = NaN(size(tf_mitgcm_boxed));
        
        for i_time=1:length(mitgcm(i_fjord).t)
            ints=[0; cumsum(h_in_mitgcm_time(:,i_time))];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            Tregular(:,i_time) = interp1(z_box,temp_in_mitgcm_time(:,i_time),mitgcm(i_fjord).z,'nearest','extrap');
            Sregular(:,i_time) = interp1(z_box,salt_in_mitgcm_time(:,i_time),mitgcm(i_fjord).z,'nearest','extrap');
        
            tf_mitgcm_boxed(:,i_time) = interp1(mitgcm(i_fjord).z,mitgcm(i_fjord).Tseries(:,i_time)',z_box,'linear','extrap');
            sf_mitgcm_boxed(:,i_time) = interp1(mitgcm(i_fjord).z,mitgcm(i_fjord).Sseries(:,i_time)',z_box,'linear','extrap');
        end

        ensemble(i_fjord,i1,i2,i3,i4).s.Tmitgcm = tf_mitgcm_boxed;
        ensemble(i_fjord,i1,i2,i3,i4).s.Smitgcm = sf_mitgcm_boxed;
        
        ensemble(i_fjord,i1,i2,i3,i4).s.Tbox   = temp_in_mitgcm_time; % we get the full (boxed) results here
        ensemble(i_fjord,i1,i2,i3,i4).s.Tupper = mean(Tregular(1:5,:),1,'omitnan');   % avg temp of the upper 50 m
        ensemble(i_fjord,i1,i2,i3,i4).s.Tinter = mean(Tregular(5:25,:),1,'omitnan');  % avg temp of the 50 - 250 m layer
        ensemble(i_fjord,i1,i2,i3,i4).s.Tlower = mean(Tregular(25:50,:),1,'omitnan'); % avg temp of the 250 - 500 m layer

        % repeat the same for salinity
        ensemble(i_fjord,i1,i2,i3,i4).s.Sbox   = salt_in_mitgcm_time; 
        ensemble(i_fjord,i1,i2,i3,i4).s.Supper = mean(Sregular(1:5,:),1,'omitnan');   
        ensemble(i_fjord,i1,i2,i3,i4).s.Sinter = mean(Sregular(5:25,:),1,'omitnan');  
        ensemble(i_fjord,i1,i2,i3,i4).s.Slower = mean(Sregular(25:50,:),1,'omitnan');

        clear temp_in_mitgcm_time salt_in_mitgcm_time % housekeeping

        fprintf('Output interpolation complete. ')
    catch ME
        fprintf('run %d failed: %s. ',run_counter,ME.message)
        cur_fjord.m.error = ME.message;
        cur_fjord.m.stack = ME.stack;
        fjords_crashed{end+1} = cur_fjord;
    end
    fprintf('\n')
end % for i4
end % for i3
end % for i2
end % for i1
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
clear i_fjord i1 i2 i3 i4

file_out = [outs_path,'boxmodel_runs_MITgcm_',flag_ice,'ice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'_60layers_dt',num2str(dt_in_h),'h'];

save(file_out,'-v7.3','ensemble','fjord_model');
save([file_out,'_crashed'],'-v7.3','fjords_crashed');
fprintf('Outputs saved in %s.\n',file_out)
% load(file_out);
%% post-processing to make things 1:1 comparable

if exist('res_obs',"var"),       clear res_obs; end
res_obs(size(ensemble)) = struct("tf",[],"sf",[],"zf",[]);
if exist('res_box',"var"),       clear res_box; end
res_box(size(ensemble)) = struct("tf",[],"sf",[],"zf",[],"ID",[],"name",[],"rmse_tf",[],"rmse_sf",[]);
for i_fjord=1:size(ensemble,1)
    zf_obs = fjord_model(i_fjord).c.zf;
    tf_obs = fjord_model(i_fjord).c.tf;
    sf_obs = fjord_model(i_fjord).c.sf;
    n_completed = 0; % completed runs for that fjord

    % n_layers = ensemble(i,1,1).p.N+ensemble(i,1,1).p.sill;
    tupper_box_comp = NaN([length(mitgcm(i_fjord).t),dims_ensemble(2:end)]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    supper_box_comp = NaN(size(tupper_box_comp));
    sinter_box_comp = NaN(size(tupper_box_comp));
    slower_box_comp = NaN(size(tupper_box_comp));

    tf_box_comp = NaN([length(zf_obs),dims_ensemble(2:end)]);
    sf_box_comp = NaN(size(tf_box_comp));
    rmse_tf = NaN(dims_ensemble(2:end));
    rmse_sf = NaN(dims_ensemble(2:end));
    tf_mitgcm_boxed = NaN([length(ensemble(i_fjord,1,1,1,1).s.Hfinal),length(mitgcm(i_fjord).t),dims_ensemble(2:end)]);
    tf_ens_mean     = NaN(size(tf_mitgcm_boxed));
    sf_mitgcm_boxed = NaN(size(tf_mitgcm_boxed));
    sf_ens_mean     = NaN(size(tf_mitgcm_boxed));
    for i1=1:size(ensemble,2)
    for i2=1:size(ensemble,3)
    for i3=1:size(ensemble,4)
    for i4=1:size(ensemble,5)
        if ~isempty(ensemble(i_fjord,i1,i2,i3,i4).s) & ~isnan(ensemble(i_fjord,i1,i2,i3,i4).s.Hfinal)
            ints=[0; cumsum(ensemble(i_fjord,i1,i2,i3,i4).s.Hfinal)];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            tf_box = ensemble(i_fjord,i1,i2,i3,i4).s.Tfinal; 
            sf_box = ensemble(i_fjord,i1,i2,i3,i4).s.Sfinal; 
            tf_mitgcm_boxed(:,:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Tmitgcm;
            sf_mitgcm_boxed(:,:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Smitgcm;
            tf_ens_mean(:,:,i1,i2,i3,i4)     = ensemble(i_fjord,i1,i2,i3,i4).s.Tbox;
            sf_ens_mean(:,:,i1,i2,i3,i4)     = ensemble(i_fjord,i1,i2,i3,i4).s.Sbox;
            
            tf_box_comp(:,i1,i2,i3,i4) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
            sf_box_comp(:,i1,i2,i3,i4) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');

            tupper_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Tupper;
            tinter_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Tinter;
            tlower_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Tlower;

            supper_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Supper;
            sinter_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Sinter;
            slower_box_comp(:,i1,i2,i3,i4) = ensemble(i_fjord,i1,i2,i3,i4).s.Slower;

            min_depth_rmse = 0;
            depths_rmse = zf_obs > min_depth_rmse;
            rmse_tf(i1,i2,i3,i4) = rmse(tf_box_comp(depths_rmse,i1,i2,i3,i4),tf_obs(depths_rmse)','omitnan')./mean(tf_obs(depths_rmse),'omitnan');
            rmse_sf(i1,i2,i3,i4) = rmse(sf_box_comp(depths_rmse,i1,i2,i3,i4),sf_obs(depths_rmse)','omitnan')./mean(sf_obs(depths_rmse),'omitnan');

            n_completed = n_completed+1;
        end
    end % i4
    end % i3
    end % i2
    end % i1
    clear i1 i2 i3 i4
    res_box(i_fjord).Tmitgcm = mean(tf_mitgcm_boxed,[3,4,5,6],'omitnan');
    res_box(i_fjord).Tbox    = mean(tf_ens_mean,[3,4,5,6],'omitnan');
    res_box(i_fjord).Smitgcm = mean(sf_mitgcm_boxed,[3,4,5,6],'omitnan');
    res_box(i_fjord).Sbox    = mean(sf_ens_mean,[3,4,5,6],'omitnan');

    res_obs(i_fjord).zf = zf_obs;
    res_obs(i_fjord).tf = tf_obs;
    res_obs(i_fjord).sf = sf_obs;
    
    res_obs(i_fjord).zs = fjord_model(i_fjord).c.zs;
    res_obs(i_fjord).ts = fjord_model(i_fjord).c.ts;
    res_obs(i_fjord).ss = fjord_model(i_fjord).c.ss;

    res_box(i_fjord).t  = mitgcm(i_fjord).t; %t(2:end); %ensemble(i,1,1,1,1).t(2:end);
    res_box(i_fjord).zf = z_box;
    res_box(i_fjord).tf = mean(tf_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).tfmin = min(tf_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).tfmax = max(tf_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).sf = mean(sf_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).sfmin = min(sf_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).sfmax = max(sf_box_comp,[],[2,3,4,5],'omitnan');


    res_box(i_fjord).Tupper = mean(tupper_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).Tupper_min = min(tupper_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).Tupper_max = max(tupper_box_comp,[],[2,3,4,5],'omitnan');

    res_box(i_fjord).Tinter = mean(tinter_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).Tinter_min = min(tinter_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).Tinter_max = max(tinter_box_comp,[],[2,3,4,5],'omitnan');

    res_box(i_fjord).Tlower = mean(tlower_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).Tlower_min = min(tlower_box_comp,[],[2,3,4,5],'omitnan');
    res_box(i_fjord).Tlower_max = max(tlower_box_comp,[],[2,3,4,5],'omitnan');

    res_box(i_fjord).Supper = mean(supper_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).Sinter = mean(sinter_box_comp,[2,3,4,5],'omitnan');
    res_box(i_fjord).Slower = mean(slower_box_comp,[2,3,4,5],'omitnan');

    res_box(i_fjord).ensemble_tf = tf_box_comp;
    res_box(i_fjord).ensemble_sf = sf_box_comp;
    res_box(i_fjord).rmse_tf = rmse_tf;
    res_box(i_fjord).rmse_sf = rmse_sf;

    res_box(i_fjord).id = fjord_model(i_fjord).m.ID{1};
    res_box(i_fjord).name = fjord_model(i_fjord).m.name{1};
    res_box(i_fjord).n = n_completed;
end % i_fjord
% housekeeping
clear tupper_box_comp tinter_box_comp tlower_box_comp supper_box_comp sinter_box_comp slower_box_comp
clear zf_obs tf_obs sf_obs tf_box_comp sf_box_comp rmse_tf rmse_sf tf_mitgcm_boxed sf_mitgcm_boxed tf_ens_mean sf_ens_mean
disp('Postprocessing done.')

%% Plotting results
plot_ensemble_results(tgt_day,res_box,res_obs,fjord_model,mitgcm,dims_ensemble,param_names,range_params,n_fjord_runs);
% exportgraphics(gcf,[figs_path,'profiles_temp_MITgcm_',flag_ice,'ice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'_',num2str(cur_fjord.p.N),'layers.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'profiles_salt_MITgcm_',flag_ice,'ice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'_',num2str(cur_fjord.p.N),'layers.png'],'Resolution',300)

plot_sensitivity_to_param(res_box,res_obs,fjord_model,param_names,range_params);
% exportgraphics(gcf,[figs_path,'sensitivity_profiles_temp_',flag_ice,'_n',num2str(n_combinations),'_day',num2str(tgt_day),'_',num2str(cur_fjord.p.N),'layers.png'],'Resolution',300)
% exportgraphics(gcf,[figs_path,'sensitivity_profiles_salt_',flag_ice,'_n',num2str(n_combinations),'_day',num2str(tgt_day),'_',num2str(cur_fjord.p.N),'layers.png'],'Resolution',300)

plot_best_params(fjord_model,res_box,dims_ensemble,param_names,range_params,tgt_day,{'peak'}); %% TODO: made this one work
% exportgraphics(gcf,[figs_path,'best_parameters_',flag_ice,'_n',num2str(n_combinations),'_day',num2str(tgt_day),'_60layers.png'],'Resolution',300)

% plot_crashed_parameters(fjords_crashed,param_names);
% plot_ens_mitgcm_boxed(res_box,fjord_model,n_fjord_runs);
% exportgraphics(gcf,[figs_path,'temp_series_box_MITgcm',num2str(n_combinations),'_3layers.png'],'Resolution',300)

%% Showcasing the different fjord geometries
% from data compilation
% for i=1:length(fjord_ids)
%     plot_compiled_fjords(fjords_map,fjords_compilation,fjord_ids(i),[figs_path,'fjords/'])
% end
