run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/pce/']; % where the model output files will be saved
figs_path = [project_path,'/figs/pce/'];                     % where the figures and animations will be saved
letters = {'a','b','c','d','e','f','g','h'};
years   = {'2016','2017','2018','2019','2020'};
fjord_names = {'NW system','Kangerlussuaq'};

file_omg_ctd = [data_path,'/greenland/obs/all_OMG_CTD_from_csv_ver_1.3.mat']; % obtained from: https://doi.org/10.1029/2021GL097081; % OMG original dataset can be accessed from: https://doi.org/10.5067/OMGEV-CTDS1 

time_start = datetime(2017,01,01);
time_end   = datetime(2020,12,31);
run compile_process_fjords % we can safely ignore the warnings, as we are using (slightly) different forcings

tgt_day = 630; % end of melt season
model_dt=0.1;
i_tgt_fjords = [9,42]; % 9 fjord system (Slorapaluup Kangerlua, Tinussarfik, Kangerluarsorujuk, Qungasissat Kangerlua, Kangerluarsuk); 42: Kangerlussuaq 
fjords_tgt   = fjords_processed(i_tgt_fjords);
ctd_data_omg = process_omg_ctd(file_omg_ctd);

%% Establishing the boundaries in space (NW fjord system; Kangerlussuaq) and time (quasi-synoptic overlap between fjord and shelf)

% NW Kanger-system
omg_bnds(1).lon_f = [-71.72, -66.];
omg_bnds(1).lon_s = [-74.05,-71.73];
omg_bnds(1).lat_f = [77.3,77.93];
omg_bnds(1).lat_s =[76.5,78.35];

% Kangerlussuaq
omg_bnds(2).lon_f = [-32.5, -31.5];
omg_bnds(2).lon_s = [-31.5,-26];
omg_bnds(2).lat_f = [68,69];
omg_bnds(2).lat_s =[67.3,68.2];

start_2016 = datetime('2016-09-15',"Format",'uuuu-MM-dd');
end_2016   = datetime('2016-09-25',"Format",'uuuu-MM-dd');
start_2017 = datetime('2017-10-10',"Format",'uuuu-MM-dd');
end_2017   = datetime('2017-10-20',"Format",'uuuu-MM-dd');
start_2018 = datetime('2018-08-20',"Format",'uuuu-MM-dd');
end_2018   = datetime('2018-08-30',"Format",'uuuu-MM-dd');
start_2019 = datetime('2019-08-10',"Format",'uuuu-MM-dd');
end_2019   = datetime('2019-08-20',"Format",'uuuu-MM-dd');
start_2019b= datetime('2019-09-01',"Format",'uuuu-MM-dd');
end_2019b  = datetime('2019-09-10',"Format",'uuuu-MM-dd');
start_2020 = datetime('2020-09-01',"Format",'uuuu-MM-dd');
end_2020   = datetime('2020-09-10',"Format",'uuuu-MM-dd');
% start_2021 = datetime('2021-08-20',"Format",'uuuu-MM-dd'); % we dont have data for 2021, so maybe not include it?
% end_2021   = datetime('2021-08-30',"Format",'uuuu-MM-dd');

fjord_or_shelf = zeros(size(ctd_data_omg)); % 1 for fjord; 2 for shelf
which_season   = zeros(size(ctd_data_omg));
which_fjord    = zeros(size(ctd_data_omg)); % 9 for NW, 42 for Kangerlussuaq
for i=1:length(ctd_data_omg)

    % check if the cast is within our region of interest
    if (ctd_data_omg(i).lat > omg_bnds(1).lat_f(1) && ctd_data_omg(i).lat < omg_bnds(1).lat_f(2)) && ...
       (ctd_data_omg(i).lon > omg_bnds(1).lon_f(1) && ctd_data_omg(i).lon < omg_bnds(1).lon_f(2))
        fjord_or_shelf(i) = 1;
        which_fjord(i) = i_tgt_fjords(1);
    elseif (ctd_data_omg(i).lat > omg_bnds(1).lat_s(1) && ctd_data_omg(i).lat < omg_bnds(1).lat_s(2)) && ...
           (ctd_data_omg(i).lon > omg_bnds(1).lon_s(1) && ctd_data_omg(i).lon < omg_bnds(1).lon_s(2))
        fjord_or_shelf(i) = 2;
        which_fjord(i) = i_tgt_fjords(1);
    end
    
    if (ctd_data_omg(i).lat > omg_bnds(2).lat_f(1) && ctd_data_omg(i).lat < omg_bnds(2).lat_f(2)) && ...
       (ctd_data_omg(i).lon > omg_bnds(2).lon_f(1) && ctd_data_omg(i).lon < omg_bnds(2).lon_f(2))
        fjord_or_shelf(i) = 1;
        which_fjord(i) = i_tgt_fjords(2);
    elseif (ctd_data_omg(i).lat > omg_bnds(2).lat_s(1) && ctd_data_omg(i).lat < omg_bnds(2).lat_s(2)) && ...
           (ctd_data_omg(i).lon > omg_bnds(2).lon_s(1) && ctd_data_omg(i).lon < omg_bnds(2).lon_s(2))
        fjord_or_shelf(i) = 2;
        which_fjord(i) = i_tgt_fjords(2);
    end

    % if it is within our region of interest, check if the cast belongs to any of the dates we are interested in
    if fjord_or_shelf(i) > 0 
        if ctd_data_omg(i).time > start_2016 & ctd_data_omg(i).time < end_2016
            which_season(i) = 2016;
        elseif ctd_data_omg(i).time > start_2017 & ctd_data_omg(i).time < end_2017
            which_season(i) = 2017;
        elseif ctd_data_omg(i).time > start_2018 & ctd_data_omg(i).time < end_2018
            which_season(i) = 2018;
        elseif (ctd_data_omg(i).time > start_2019 & ctd_data_omg(i).time < end_2019) | ...
                (ctd_data_omg(i).time > start_2019b & ctd_data_omg(i).time < end_2019b) 
            which_season(i) = 2019;
        elseif ctd_data_omg(i).time > start_2020 & ctd_data_omg(i).time < end_2020
            which_season(i) = 2020;
        % elseif ctd_data_omg(i).time > start_2021 & ctd_data_omg(i).time < end_2021
        %     which_season(i) = 6;
        end
    end
end

compiled_omg_casts = ctd_data_omg(which_season>0&fjord_or_shelf>0);
compiled_omg_years = which_season(which_season>0&fjord_or_shelf>0);
which_fjord        = which_fjord(which_season>0&fjord_or_shelf>0);
fjord_or_shelf     = fjord_or_shelf(which_season>0&fjord_or_shelf>0);

%% Get fjord attributes for each target fjord from data compilation, but swap their original (reanalysis-based) ocean forcing

if exist('fjord_model',"var"),clear fjord_model; end
fjord_model(length(years)*length(fjords_tgt)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
i_model=1;
for i_season=1:length(years)
    season_yr = str2double(years(i_season));
    for i=1:length(fjords_tgt)

        shelf_casts = compiled_omg_casts(compiled_omg_years==season_yr & fjord_or_shelf==2 & which_fjord==i_tgt_fjords(i));
        fjord_casts = compiled_omg_casts(compiled_omg_years==season_yr & fjord_or_shelf==1 & which_fjord==i_tgt_fjords(i));

        if ~isempty(shelf_casts) && ~isempty(fjord_casts)
            % replacing the ocean forcing
            t = 0:model_dt:365*2; % 2 years of simulation
            if length(shelf_casts) > 1
    
                % if there are several casts, finds deepest profile, and
                % interpolates all profiles to 1m resolution, as deep as the deepest profile
                deepest_point = 0;
                for k=1:length(shelf_casts)
                    if max(abs(shelf_casts(k).depth)) > deepest_point, deepest_point = max(abs(shelf_casts(k).depth)); end
                end
                z_common = 0:1:deepest_point;
                t_common = NaN([length(shelf_casts),length(z_common)]);
                s_common = NaN([length(shelf_casts),length(z_common)]);
                for k=1:length(shelf_casts)
                    t_common(k,:) = interp1(shelf_casts(k).depth,shelf_casts(k).temp,z_common,'nearest','extrap');
                    s_common(k,:) = interp1(shelf_casts(k).depth,shelf_casts(k).sal,z_common,'nearest','extrap');
                end
    
                % then takes the average of all (interpolated and standardised) profiles as our ocean forcing
                ts = mean(t_common,1,'omitnan');
                ss = mean(s_common,1,'omitnan');
                zs = -z_common;
            else
                ts = shelf_casts.temp;
                ss = shelf_casts.sal;
                zs = -shelf_casts.depth';
            end
            

            if max(abs(zs)) < fjords_tgt(i).p.H % if our profile does not extend to the bottom of the fjord
                h_missing = fjords_tgt(i).p.H - max(abs(zs));
                depth_range_missing = max(abs(zs))+1:1:fjords_tgt(i).p.H;
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
            H_clim = double(get_fjord_boxes_from_density(f.Ts',f.Ss',f.zs,fjords_tgt(i).p)); % negative values, top to bottom
            [T_clim,S_clim] = bin_ocean_profiles(mean(f.Ts,1),mean(f.Ss,1),f.zs,H_clim,fjords_tgt(i).p);
            a.H0 = H_clim;
            a.T0 = T_clim;
            a.S0 = S_clim;
        
            % flipping the ocean forcing like the model wants it
            f.zs = flip(f.zs);
            f.Ts = flip(f.Ts,2)'; 
            f.Ss = flip(f.Ss,2)';

            % Icebergs - MISSING
    
    
            fjord_model(i_model).p = fjords_tgt(i).p;
            fjord_model(i_model).a = fjords_tgt(i).a;
            fjord_model(i_model).f = fjords_tgt(i).f;
            fjord_model(i_model).m = fjords_tgt(i).m;

            fjord_model(i_model).f.Qsg = interp1(fjords_tgt(i).t,fjords_tgt(i).f.Qsg,t,'linear');
            fjord_model(i_model).f.D   = interp1(fjords_tgt(i).t,fjords_tgt(i).f.D,t,'linear');
            fjord_model(i_model).f.Ts = f.Ts;
            fjord_model(i_model).f.Ss = f.Ss;
            fjord_model(i_model).f.zs = f.zs;
            fjord_model(i_model).a.T0 = a.T0;
            fjord_model(i_model).a.S0 = a.S0;
            fjord_model(i_model).a.H0 = a.H0;
            fjord_model(i_model).t = t;
            fjord_model(i_model).m.name = fjord_names{i};
            fjord_model(i_model).m.date = shelf_casts(1).time;
            fjord_model(i_model).c.zs = zs;
            fjord_model(i_model).c.ts = ts;
            fjord_model(i_model).c.ss = ss;

            if length(fjord_casts) > 1
    
                % if there are several casts, finds deepest profile, and
                % interpolates all profiles to 1m resolution, as deep as the deepest profile
                deepest_point = 0;
                for k=1:length(fjord_casts)
                    if max(abs(fjord_casts(k).depth)) > deepest_point, deepest_point = max(abs(fjord_casts(k).depth)); end
                end
                z_common = 0:1:deepest_point;
                t_common = NaN([length(fjord_casts),length(z_common)]);
                s_common = NaN([length(fjord_casts),length(z_common)]);
                for k=1:length(fjord_casts)
                    t_common(k,:) = interp1(fjord_casts(k).depth,fjord_casts(k).temp,z_common,'nearest','extrap');
                    s_common(k,:) = interp1(fjord_casts(k).depth,fjord_casts(k).sal,z_common,'nearest','extrap');
                end
    
                % then takes the average of all (interpolated and standardised) profiles as our ocean forcing
                fjord_model(i_model).c.tf = mean(t_common,1,'omitnan');
                fjord_model(i_model).c.sf = mean(s_common,1,'omitnan');
                fjord_model(i_model).c.zf = -z_common;
            else
                fjord_model(i_model).c.tf = fjord_casts.temp;
                fjord_model(i_model).c.sf = fjord_casts.sal;
                fjord_model(i_model).c.zf = -fjord_casts.depth';
            end
    
            % % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
            figure('Position',[100 100 800 500]);
            ints=[0,-cumsum(a.H0)];
            subplot(2,3,1); hold on; box on;
            plot(f.Ts(:,1),f.zs);
            scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
            for j=1:size(ints,2)
                yline(ints(j));
            end
            yline(fjord_model(i_model).p.zgl,'--k','linewidth',0.5)
            yline(fjord_model(i_model).p.silldepth,'--k','linewidth',0.5)
            ylim([-fjord_model(i_model).p.H 0])
            subplot(2,3,2); hold on; box on;
            plot(fjord_model(i_model).f.Ss(:,1),fjord_model(i_model).f.zs);
            scatter(fjord_model(i_model).a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
            for j=1:size(ints,2)
                yline(ints(j));
            end
            yline(fjord_model(i_model).p.zgl,'--k','linewidth',0.5)
            yline(fjord_model(i_model).p.silldepth,'--k','linewidth',0.5)
            ylim([-fjord_model(i_model).p.H 0])
            subplot(2,3,3); hold on; box on; grid on;
            plot(fjord_model(i_model).t,fjord_model(i_model).f.Qsg);
            plot(fjord_model(i_model).t,fjord_model(i_model).f.D);
            legend('Subglacial discharge','Solid-ice discharge');

            subplot(2,3,4); hold on; box on;
            imagesc(fjord_model(i_model).t,fjord_model(i_model).f.zs,fjord_model(i_model).f.Ts)
            ylim([-fjord_model(i_model).p.H 0]);
            xlabel('Model time'); ylabel('Temperature (^oC)')
            subplot(2,3,5); hold on; box on;
            imagesc(fjord_model(i_model).t,fjord_model(i_model).f.zs,fjord_model(i_model).f.Ss)
            ylim([-fjord_model(i_model).p.H 0]);
            xlabel('Model time'); ylabel('Salinity')
            
        end
        i_model=i_model+1;
    end
    
end

%% Running the model

n_fjord_runs = length(fjord_model);
range_C0 = [1e3,5e3,1e4,5e4,1e5];
range_P0 = 5:5:30;
range_M0 = [5e-9,1e-8,2e-8,5e-8];
range_gamma = 0.4:0.05:0.6;
n_C0 = length(range_C0);
n_P0 = length(range_P0);
n_M0 = length(range_M0);
n_gamma = length(range_gamma);
n_combinations=n_gamma*n_M0*n_P0*n_C0;
if exist('ensemble',"var"),clear ensemble; end
% ensemble(n_fjord_runs, length(range_C0), length(range_P0)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[],"s",[]);
ensemble(n_fjord_runs, length(range_C0), length(range_P0),length(range_M0),length(range_gamma)) = struct("p",[],"t",[],"m",[],"s",[]);
run_counter=0;
for i_fjord=1:n_fjord_runs
tic
for i_C0=1:n_C0
for i_P0=1:n_P0
for i_M0=1:n_M0
for i_gamma=1:n_gamma
    run_counter = run_counter+1;
    cur_fjord = fjord_model(i_fjord);
    cur_fjord.p.C0 = range_C0(i_C0);
    cur_fjord.p.P0 = range_P0(i_P0);
    cur_fjord.p.M0 = range_M0(i_M0);
    cur_fjord.p.gamma = range_gamma(i_gamma);
    % cur_fjord.p.plot_runtime=1;
    try
        [cur_fjord.s,cur_fjord.f] = boxmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
        fprintf('run %d complete. ',run_counter)

        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).p = cur_fjord.p;
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).t = cur_fjord.t;
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).m = cur_fjord.m;
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tfinal = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)/model_dt),2); % 10-day avg centered at the target day
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Sfinal = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)/model_dt),2);
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Hfinal = mean(cur_fjord.s.H(:,(tgt_day-5:tgt_day+5)/model_dt),2);

        Tregular=NaN([length(cur_fjord.c.zs),length(cur_fjord.s.t)]);
        for k=1:length(cur_fjord.s.t)
            ints=[0; cumsum(cur_fjord.s.H(:,k))];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            Tregular(:,k) = interp1(z_box,cur_fjord.s.T(:,k),cur_fjord.c.zs,'nearest','extrap');
        end
        
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tupper = mean(Tregular(1:5,:),1,'omitnan');   % avg temp of the upper 50 m
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tinter = mean(Tregular(5:25,:),1,'omitnan');  % avg temp of the 50 - 250 m layer
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tlower = mean(Tregular(25:50,:),1,'omitnan'); % avg temp of the 250 - 500 m layer
        fprintf('Output interpolation complete. ')
    catch ME
        fprintf('run %d failed: %s. ',run_counter,ME.message)
    end
    fprintf('\n')
end % for i_gamma
end % for i_M0
end % for i_P0
end % for i_C0
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
save([outs_path,'boxmodel_runs_OMG_2fjords_comp_n',num2str(n_combinations),'_day',num2str(tgt_day)],'-v7.3','ensemble','fjord_model');
disp('Outputs saved.')

%% post-processing to make things 1:1 comparable


if exist('res_obs',"var"),       clear res_obs; end
res_obs(size(ensemble)) = struct("tf",[],"sf",[],"zf",[]);
if exist('res_box',"var"),       clear res_box; end
res_box(size(ensemble)) = struct("tf",[],"sf",[],"zf",[],"ID",[],"name",[]);
for i=1:size(ensemble,1)
    if isempty(fjord_model(i).f), continue; end
    zf_obs = fjord_model(i).c.zf;
    tf_obs = fjord_model(i).c.tf;
    sf_obs = fjord_model(i).c.sf;
    n_completed = 0; % completed runs for that fjord

    % n_layers = ensemble(i,1,1).p.N+ensemble(i,1,1).p.sill;
    tf_box_comp = NaN([length(zf_obs),n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    sf_box_comp = NaN(size(tf_box_comp));
    tupper_box_comp = NaN([length(fjord_model(i).t)-1,n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    for j=1:size(ensemble,2)
    for k=1:size(ensemble,3)
    for l=1:size(ensemble,4)
    for n=1:size(ensemble,5) % we skip 'm' because it is being used for the metadata structure elsewhere in the code
        if ~isempty(ensemble(i,j,k,l,n).s) & ~isnan(ensemble(i,j,k,l,n).s.Hfinal)
            ints=[0; cumsum(ensemble(i,j,k,l,n).s.Hfinal)];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            tf_box = ensemble(i,j,k,l,n).s.Tfinal; 
            sf_box = ensemble(i,j,k,l,n).s.Sfinal; 
            
            tf_box_comp(:,i,j,k,l,n) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
            sf_box_comp(:,i,j,k,l,n) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');

            tupper_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tupper;
            tinter_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tinter;
            tlower_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tlower;

            n_completed = n_completed+1;
        end
    end % n
    end % l
    end % k
    end % j
    res_obs(i).zf = zf_obs;
    res_obs(i).tf = tf_obs;
    res_obs(i).sf = sf_obs;
    res_obs(i).zs = fjord_model(i).c.zs;
    res_obs(i).ts = fjord_model(i).c.ts;
    res_obs(i).ss = fjord_model(i).c.ss;

    res_box(i).t  = ensemble(i,1,1,1,1).t(2:end);
    res_box(i).zf = z_box;
    res_box(i).tf = median(tf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).tfmin = min(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).tfmax = max(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).sf = median(sf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');


    res_box(i).Tupper = median(tupper_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tupper_min = min(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tupper_max = max(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tinter = median(tinter_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tinter_min = min(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tinter_max = max(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tlower = median(tlower_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tlower_min = min(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tlower_max = max(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    % res_box(i).id = fjord_model(i).m.ID{1};
    fjord_model(i).m.date.Format = 'dd-MM-uuuu';
    res_box(i).name = [fjord_model(i).m.name,' ',string(fjord_model(i).m.date)];
    res_box(i).n = n_completed;
end

%% Plotting results
lcolor = lines(3);

for i=1:n_fjord_runs
    if ~isempty(res_box(i)) & res_box(i).n>0
        figure('Name','Temperature comparison','Position',[40 40 1200 800]);
        tiledlayout(1,2);
        nexttile; hold on; box on; grid on
        text(0.02,1.05,sprintf("Fjord: %s",res_box(i).name),'units','normalized','fontsize',14)
        text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',14)
    
        % figure; hold on
        % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
        % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
        % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);
    
        hs = plot(res_obs(i).ts,res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
        hf = plot(res_obs(i).tf,res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
        hb = plot(res_box(i).tf,res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
        plot(res_box(i).tfmin,res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
        plot(res_box(i).tfmax,res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
        xlabel('Temperature (^oC)'); ylabel('Depth (m)');
        set(gca,'fontsize',14)
        xlim([-2 6])
        hleg = legend([hs, hf, hb], {"Shelf","Fjord","Box model"},'fontsize',14,'Location','Southeast'); 
        title(hleg,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
    
        nexttile; hold on; box on; grid on
        if length(res_box(i).t) == length(res_box(i).Tupper)
            hu = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
            hi = plot(res_box(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
            hl = plot(res_box(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
        
            % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
            % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
           
        end
        ylabel('Box model Temperature (^oC)'); xlabel('Model time (days)');
        set(gca,'fontsize',14)

        hleg = legend([hu,hi,hl],{"0-50 m","50-250 m","250-500 m"},'fontsize',14,'Location','Northeast'); 
        title(hleg,'Time series')
    end
end


% exportgraphics(gcf,[figs_path,'temp_profiles_OMG_2fjords_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'.png'],'Resolution',300)
