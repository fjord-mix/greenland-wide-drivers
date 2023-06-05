%% Script to run, save, and plot outputs from the boxmodel
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))

[datasets,fjords_processed,fjords_map,msk_fjord_ids] = compile_datasets(data_path);

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/']; % where the model output files will be saved
figs_path = [project_path,'/figs/'];                     % where the figures and animations will be saved
%% Summary of data compilation 
% useful to check the link of glaciers and shelf profiles with their respective fjord and if all data was properly loaded
plt_handles     = plot_fjords_summary(datasets,fjords_map,fjords_processed);
plt_obs_handles = plot_obs(datasets);% ,datasets.obs.ctd_data.omg); % add this second argument to plot only a specific obs dataset

% just to facilitade figure browsing
%plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
%plt_obs_handles.cb1.Visible = 'off'; plt_obs_handles.cb2.Visible = 'off'; plt_obs_handles.cb3.Visible = 'off'; 

%% Choosing model runtime
datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2020,01,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 

%% Setting up specific "flagship fjords" which might be of interest
% Their IDs were obtained visually using plot_fjords_summary(). uncomment
% the "text()" function line where IDs are plotted to see them
% 22 + 23 + 24 + 25: Uummannaq system; Kangerlussuup Sermia drains into 25
% 34: Ilulissat (Jakobshavn)
% 35: Kangerlussuaq
% 39: Sermilik (Helheim)
% 45: Kangersuneq (Kangiata Nunaata Sermia; Nuuk's fjord?)
fjord_ids=[35,39,45,34];
fjord_names={'Kangerlussuaq','Sermilik','Kangersuneq','Ilulissat'};
fjord_keys={'KF','SF','KS','IS'};

%% Exploring parametre space

run exps_parametre_space.m

%% Running for all fjords of interest

name_exp = '';
output_fname=['test_fjords_',name_exp];
run exps_flagship_fjords.m
% save([outs_path,output_fname],'fjord_model','-v7.3') % v7.3 allows files > 2GB

%% A more statistically driven approach to sensitivity tests

% Obtains some general statistics. Might be useful in case we go for
% some type of emulator/surrogate model to carry sensitivity tests
% verbose.print = 1;
% verbose.plot  = 0;
% fjord_stats = print_fjord_statistics(fjords_processed,verbose);
% 
% fjord_model_idealised = prepare_idealised_boxmodel_input(fjord_stats,verbose);
