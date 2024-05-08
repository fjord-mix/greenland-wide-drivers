%% Sets up to run the box model within the UQLab framework using CTD casts as the ocean forcing

%% Setting paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/emulation/']; % where the model output files will be saved
figs_path = [project_path,'/figs/emulation/'];                     % where the figures and animations will be saved
letters = {'a','b','c','d','e','f','g','h'};
dep_bnds = [50, 250, 500]; % depth boundaries for the different layers to be emulated


%% Compiling datasets to be used. 

% Note that "load_default_opts.m" should be checked for which large datastets are actually in use here
[datasets,fjords_compilation,fjords_map,~,glaciers_compilation] = compile_datasets(data_path);

% Getting all OMG casts and sorting them by their date
casts = process_omg_ctd(datasets.obs.file_omg_ctd); % gets OMG casts
T                = struct2table(casts);             % convert the struct array to a table
sortedT          = sortrows(T, 'time');
datasets.obs.ctd = table2struct(sortedT);
[~,~,datasets.obs.ctd_shelf] = add_fjord_shelf_flag(datasets.obs.ctd,fjords_map); % get only the shelf casts

clear casts T sortedT % housekeeping
%% Initialise all needed variables

n_train    = 700; % runs for training the emulator
n_emul    = 1e5; % sample size for the emulator
time_step = 0.1; % in days
n_layers = length(dep_bnds);

if datasets.opts.load_ocean && ~isfield(datasets.obs,'ctd_shelf'), n_inputs=11; else, n_inputs=10; end

% initialising training dataset structures
time_axis = datetime(2010,01,15):1:datetime(2018,12,15);
if exist('ensemble',"var"), clear ensemble; end
ensemble(n_train) = struct("time",[],"temp",[],"salt",[],"H",[],"ts",[],"ss",[],"zs",[],"p",[],"phi",[],"qvs",[],"qsg",[],"qts",[]);

X          = zeros([n_train,n_inputs]);
Xemul      = zeros([n_emul,n_inputs]);

uqlab % Initialise UQLab

%% sample the distributions for inputs
rng('default') % set the seed for reproducibility
datasets.opts.restrict_to_fjords = 0; % whether we want glacier distributions to only include 
                                      % glacier that drains into the compiled fjords 
                                      % or all in the region
tic
fprintf('Generating inputs... ')
[Parameters,IOpts,probs,~] = define_model_param_distrib(datasets,fjords_compilation,glaciers_compilation,[],time_axis,time_step); % available outputs: [Parameters,IOpts,probs,fjords_processed]
input = uq_createInput(IOpts); % create probability functions
X      = uq_getSample(input,n_train,'LHS');  % perform latin hypercube sampling of our parametre space for the training dataset
Xemul  = uq_getSample(input,n_emul,'LHS');   % same, but for the surrogate model (hence the much larger n)
disp(' Done.')
toc

% hf = plot_inputs_summary(IOpts,probs); % Plotting to see how the distributions look
%% Performing the numerical model runs (i.e., training dataset)

tic
for k_run=1:n_runs
    % if isnan(ohc_out(k_run,i_reg)) % saves time after a bug/crash fix, but requires reinitialising the variable        
    ensemble(k_run) = wrapper_boxmodel(X(k_run,:),Parameters);
    fprintf('run %d complete\n',k_run)
    % end
end
fprintf('Model training runs complete.\n')
toc    

fprintf('Computing dT for the different layers...\n')
dT_out = compute_fjord_shelf_differences(ensemble,'T');
% [ohc_out,osc_out] = compute_ensemble_metric(ensemble,length(time_axis));
fprintf('Done.\n')