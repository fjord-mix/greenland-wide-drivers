%% Script to run, save, and plot outputs from the boxmodel
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and figs_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
% addpath(genpath(pwd)) % adds subfolders of local directory
% addpath(genpath(project_path))

[datasets,fjords_processed,fjords_map,msk_fjord_ids] = compile_datasets(data_path);

% this is where the model output files will be saved
fjord_model_input_folder = [data_path,'/greenland/FjordMIX/boxmodel/'];

%% Summary of data compilation 
% useful to check the link of glaciers and shelf profiles with their respective fjord and if all data was properly loaded
plt_handles     = plot_fjords_summary(datasets,fjords_map,fjords_processed);
plt_obs_handles = plot_obs(datasets);% ,datasets.obs.ctd_data.omg); % add this second argument to plot only a specific obs dataset

% just to facilitade figure browsing
%plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
%plt_obs_handles.cb1.Visible = 'off'; plt_obs_handles.cb2.Visible = 'off'; plt_obs_handles.cb3.Visible = 'off'; 

%% Running for specific "flagship fjords"
% 22 + 23 + 24 + 25: Uummannaq; Kangerlussuup Sermia drains into 25 (obs by Jackson et al., 2017)
% 34: Ilulissat (Jakobshavn)
% 35: Kangerlussuaq
% 39: Sermilik (Helheim)
% 45: Kangersuneq (Kangiata Nunaata Sermia; Nuuk's fjord?)

%Options for runtime period and time step
datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2020,01,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt = 0.05/24.; % time step in days(?)

%% Model parametres to explore
% p.P0 = 0; % no plume with entrainment efficiency being zero (default 25)
% p.C0 = 1e4; % shelf exchange efficiency (s) also limits the exchange from shelf (Default: 1e4; works with 1e3 for 2010-2017)
% p.trelax = 10*86400; % nudging relaxation time in seconds (n_days*seconds_in_a_day)



close all
fjord_ids=[35,39,45,34]; % select which fjords to run (check their IDs with plot_fjords_summary() )
fjord_names={'Kangerlussuaq','Sermilik','Kangersuneq','Ilulissat'};

name_exp = 'rlx10y_dt0p05';
output_fname=['test_fjords_',name_exp];
if exist('fjord_model','Var'), clear fjord_model; end
fjord_model(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
for id=1:2%length(fjord_ids)
    fjord_name = fjord_names{id};    
    fprintf('Running model for %s... ',fjord_name)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = fjord_name;   
    fjord_model(id) = fjord_run;
    % plot_outputs(fjord_run)
end
save([fjord_model_input_folder,output_fname],'fjord_model','-v7.3') % v7.3 allows files > 2GB
for id=1:2, plot_outputs(fjord_model(id)); end %plot_fluxes(fjord_model(id)); end

%% Loop to plot model outputs for desired fjords
fjord_model = load([fjord_model_input_folder,output_fname]).fjord_model;
figs_path='/Users/mmeb1/FjordMIX/figs';
for id=1:2%length(fjord_model)
    % fprintf('Plotting model outputs for %s... ',fjord_model(id).m.name)    
    % plot_figures(datasets,fjord_model(id),fjords_map, datasets.obs.ctd_data.omg)
    % disp('Done.')
    if ~fjord_model(id).s.status
        disp('Creating video...')
        animate_v4p1(fjord_model(id),figs_path,[output_fname,'_',fjord_model(id).m.name],50);
        disp('Done.')
    end
end
%% Running the model for sensitivity tests

% Obtains some general statistics. Might be useful in case we go for
% some type of emulator/surrogate model to carry sensitivity tests
% verbose.print = 1;
% verbose.plot  = 0;
% fjord_stats = print_fjord_statistics(fjords_processed,verbose);
% 
% fjord_model_idealised = prepare_idealised_boxmodel_input(fjord_stats,verbose);
