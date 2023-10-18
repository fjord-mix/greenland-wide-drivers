%% Sets up to run the box model within the UQLab framework

%% Configuring paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

[datasets,fjords_compilation,fjords_map,~] = compile_datasets(data_path);
outs_path = [data_path,'/greenland/FjordMIX/boxmodel/pce/']; % where the model output files will be saved
figs_path = [project_path,'/figs/pce/'];                     % where the figures and animations will be saved
letters = {'a','b','c','d','e','f','g','h'};
regions = {'SW','SE','CW','CE','NW','NE','NO'};

%% Showing all input parameters first

% close all
% plot_reg_ocn_forcings(datasets,fjords_compilation)
% figure(1); exportgraphics(gcf,[figs_path,'hovmoller_Ts.png'],'Resolution',300)
% figure(2); exportgraphics(gcf,[figs_path,'hovmoller_Ss.png'],'Resolution',300)
% figure(3); exportgraphics(gcf,[figs_path,'series_discharge_hc_sc.png'],'Resolution',300)

% hs = plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
% hf = plot_distributions(datasets,fjords_compilation);
% exportgraphics(hf,[figs_path,'summary_input_probs2010-2018.png'],'Resolution',300)

%% Initialise all needed variables
n_runs    = 200;
time_step = 0.1; % in days
n_regions = length(regions);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
ohc_out_as = NaN([n_runs, n_regions]);
osc_out_as = NaN([n_runs, n_regions]);
ensemble(n_runs,n_regions) = struct("time",[],"ohc",[],"osc",[],"ohc_as",[],"osc_as",[]);
ohc_pd  = cell([1,n_regions]);
osc_pd  = cell([1,n_regions]);
ohc_ks  = cell([1, n_regions]);
osc_ks  = cell([1, n_regions]);
ohc_ks_as  = cell([1, n_regions]);
osc_ks_as  = cell([1, n_regions]);

Parameters = cell([1, n_regions]);
X          = zeros([n_runs,n_regions,10]);

uqlab % Initialise UQLab

%% sample the distributions for inputs
rng('default') % set the seed for reproducibility
for i_reg=1:n_regions
    [Parameters{i_reg},IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg,time_step); % available outputs: [Parameters,IOpts,probs,fjords_processed]
    input = uq_createInput(IOpts);
    X(:,i_reg,:) = uq_getSample(input,n_runs,'LHS'); % perform latin hypercube sampling of our parametre space
end

%% Run the model for all iterations
% load([outs_path,'hc_sc_ensemble_n',num2str(n_runs)],'ensemble') % if we have the results saved already

run model_runs_per_region_and_pdfs.m
save([outs_path,'hc_sc_ensemble_n',num2str(n_runs)],'-v7.3','ensemble') % save ensemble structure so we do not need to rerun it all the time

% save outputs so we dont have to re-run it
save([outs_path,'ohc_osc_change_runs_probs_n',num2str(n_runs)],...
      'ohc_out','osc_out','ohc_pd','osc_pd','ohc_ks','osc_ks')
%% Setting up the PCE model per region using UQLab

run compute_surrogate_and_sobol_indices.m

% ideally we would save the PCE-related variables as well, but Matlab always returns an error
% when trying to save them
% save([outs_path,'ohc_osc_pce_n50_n1e6'],'sur_model_ohc','sur_model_osc','Ysur_ohc','Ysur_osc','Yeval_ohc','Yeval_osc','sobolA_ohc','sobolA_osc')

%% Plotting the outputs

run make_figs_uqlab_outputs.m % figures are already saved inside the function