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


%% Initialise all needed variables
n_runs    = 700; % runs for producing the surrogate model
% n_valid   = floor(n_runs/10); % independent runs for surrogate model validation
n_surr    = 1e6; % sample size for the surrogate model itself
time_step = 0.1; % in days
n_regions = length(regions);
time_axis = datetime(2010,01,15):1:datetime(2018,12,15);

% initialising traiing and validation dataset structures
if exist('ensemble',"var"),       clear ensemble; end
if exist('ensemble_valid',"var"), clear ensemble_valid; end
ensemble(n_runs,n_regions)        = struct("time",[],"temp",[],"salt",[],"H",[],"ts",[],"ss",[],"zs",[],"p",[],"phi",[],"qvs",[],"qsg",[]);
% ensemble_valid(n_valid,n_regions) = struct("time",[],"temp",[],"salt",[],"H",[],"ts",[],"ss",[],"zs",[],"p",[],"phi",[],"qvs",[]);

Parameters = cell([1, n_regions]);
IOpts      = cell([1, n_regions]);
X          = zeros([n_runs,n_regions,10]);
% Xvalid     = zeros([n_valid,n_regions,10]);
Xeval     = zeros([n_surr,n_regions,10]);

uqlab % Initialise UQLab

%% sample the distributions for inputs
rng('default') % set the seed for reproducibility
for i_reg=1:n_regions
    tic
    fprintf('Generating inputs for %s...\n',regions{i_reg})
    [Parameters{i_reg},IOpts{i_reg},~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg,time_axis,time_step); % available outputs: [Parameters,IOpts,probs,fjords_processed]
    input = uq_createInput(IOpts{i_reg}); % create probability functions

    % perform latin hypercube sampling of our parametre space
    X(:,i_reg,:)      = uq_getSample(input,n_runs,'LHS');  % training dataset
    % Xvalid(:,i_reg,:) = uq_getSample(input,n_valid,'LHS'); % validation dataset
    Xeval(:,i_reg,:)  = uq_getSample(input,1e6,'LHS'); % surrogate model
    toc
end
% clear input

%% Run the model for all iterations

run run_model_compute_pdfs.m

% save outputs so we dont have to re-run it
save([outs_path,'ts_mp_ensembles_n',num2str(n_runs),''],'-v7.3','X','Xeval','Parameters','ensemble') % save ensemble structure so we do not need to rerun it all the time
save([outs_path,'ts_mp_diffs_n',num2str(n_runs)],'ohc_out','osc_out'); %,'ohc_pd','osc_pd','ohc_ks','osc_ks')
disp('Numerical model outputs saved.')
%% Setting up the PCE model per region using UQLab
% if we have the results saved already
% load([outs_path,'ts_mp_ensembles_n',num2str(n_runs),'']) 
% load([outs_path,'ts_mp_diffs_n',num2str(n_runs),''])
tic
run compute_surrogate_and_sobol_indices.m
toc
% ideally we would save the PCE-related variables as well, but the file would be too large
% save([outs_path,'pce_dt_ds_n',n_runs],'-v7.3','sur_model_ohc','sur_model_osc'),'Ysur_ohc','Ysur_osc','Yeval_ohc','Yeval_osc','sobolA_ohc','sobolA_osc','BorgonovoA_ohc','BorgonovoA_osc')

%% Plotting the outputs

run make_figs_uqlab_outputs.m % figures are already saved inside the script