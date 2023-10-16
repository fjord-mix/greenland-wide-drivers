%% Script to run, save, and plot outputs from the boxmodel
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/']; % where the model output files will be saved
figs_path = [project_path,'/figs/'];                     % where the figures and animations will be saved

[datasets,fjords_processed,fjords_map,msk_fjord_ids] = compile_datasets(data_path);

%% Choosing model runtime
datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,01,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 1.0; % time step in days

%% Summary of data compilation 
% useful to check the link of glaciers and shelf profiles with their respective fjord and if all data was properly loaded

% plt_handles     = plot_fjords_summary(datasets,fjords_map,fjords_processed); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
% plt_obs_handles = plot_obs(datasets);
% plt_obs_handles = plot_obs(datasets,datasets.obs.ctd_data.omg); % plots only a specific obs dataset
%plt_obs_handles.cb1.Visible = 'off'; plt_obs_handles.cb2.Visible = 'off'; plt_obs_handles.cb3.Visible = 'off'; 

%% Setting up specific "flagship fjords" which might be of interest
% Their IDs were obtained visually using plot_fjords_summary()
% Uncomment the "text()" function line where IDs are plotted to see them
% 22 + 23 + 24 + 25: Uummannaq system; Kangerlussuup Sermia drains into 25
% 34: Ilulissat (Jakobshavn)
% 35: Kangerlussuaq
% 39: Sermilik (Helheim)
% 45: Kangersuneq (Kangiata Nunaata Sermia; Nuuk's fjord?)
fjord_ids=[35,39,45,34];
fjord_names={'Kangerlussuaq','Sermilik','Kangersuneq','Ilulissat'};
fjord_keys={'KF','SF','KS','IS'};

%% Control experiment, to see how it works
id=1;
name_ctrl = sprintf('%s_ctrl',fjord_keys{id});
fjord_ctrl = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
fjord_ctrl.p.plot_runtime=1;
[fjord_ctrl.s,fjord_ctrl.f] = boxmodel(fjord_ctrl.p,fjord_ctrl.t,fjord_ctrl.f,fjord_ctrl.a); % runs and gets the results    
fjord_ctrl.o = postprocess_boxmodel(fjord_ctrl);
fjord_ctrl.m.name = name_ctrl;
% save([outs_path,'/',fjord_ctrl.m.name],'fjord_ctrl','-v7.3') % v7.3 allows files > 2GB

% Examples of different ways to plot the outputs
plot_outputs(fjord_ctrl);
plot_ts_at_depth(fjord_ctrl,[5,200,500],'nearest');
plot_fluxes(fjord_ctrl);
plot_fw_out(fjord_ctrl); % still WiP

%% Sample runs for some fjords of interest
exp_out_path='/test_benchmark_fjords/';
mkdir(exp_out_path)
run exps_benchmark_fjords.m
% run exps_benchmark_fjords_par.m % not worth running in parallel for only 4 runs
% save([exp_out_path,'example_benchmark_fjords_',num2str(p.N),'layers'],'fjord_model','-v7.3') % v7.3 allows files > 2GB

% Showing what can go wrong
run exps_bad_examples.m

%% Exploring parametre space

exp_out_path=[outs_path,'sens_model_params/'];
mkdir(exp_out_path)
run exps_parameter_space.m

exp_out_path='/sens_icebergs/';
mkdir(exp_out_path)
run exps_icebergs.m