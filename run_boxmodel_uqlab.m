%% Sets up to run the box model within the UQLab framework

%% Configuring paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))

[datasets,fjords_compilation,~,~] = compile_datasets(data_path);

% outs_path = [data_path,'/greenland/FjordMIX/boxmodel/pce/']; % where the model output files will be saved
% figs_path = [project_path,'/figs/pce/'];                     % where the figures and animations will be saved


%% evaluate the wrapper by itself first
[Parameters,IOpts,probs,fjords_processed] = define_model_param_distrib(datasets,fjords_compilation,1);

rng('default')
n_runs = 10;
clear ensemble fjord run
ensemble(length(n_runs)) = struct('time',[],'ohc',[],'osc',[]);
% ohc_out = NaN([n_runs, 1]);
for k=1:n_runs
try
    X = zeros(size(probs));
    for i=1:length(probs)
        X(i) = random(probs(i));
    end    
    tic
    [ensemble(k).time,ensemble(k).ohc,ensemble(k).osc,status,fjord_run(k)] = wrapper_boxmodel(X,Parameters);
    % ohc_out(k) = wrapper_boxmodel(X,Parameters);
    toc
catch ME
    fprintf('Crash on iteration #%d\n',k)
    if ~isempty(ME.stack) > 0
        fprintf('Something went wrong on %s line %d\n',ME.stack(1).name,ME.stack(1).line)
        fprintf('The error message was:\n%s\n',ME.message)
    end
end
end
figure; hold on; for k=1:n_runs, plot(ensemble(k).time,ensemble(k).ohc); end

%% Setting up the runs per region

for i_reg=1:1
    [Parameters,IOpts,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);

    % Initialise UQLab
    uqlab

    % create model
    % MOpts.mFile = 'uq_ishigami';
    ModelOpts.mFile = 'wrapper_boxmodel';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters=Parameters;
    num_model = uq_createModel(ModelOpts);
    
    % create inputs    
    input = uq_createInput(IOpts);
    X = uq_getSample(input,10,'LHS');

    % create Meta Model: specification of 14th degree LARS−based PCE     
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';    
    MetaOpts.FullModel = num_model;
    MetaOpts.Method = 'LARS';
    MetaOpts.Degree = 14;
    MetaOpts.ExpDesign.NSamples = 10;
    model=uq_createModel(MetaOpts);
    % uq_print(model)
    % uq_display(model)

    Y = uq_evalModel(model,X);
    Yvalidation = uq_evalModel(num_model,X);
    uq_plot(Yvalidation,Y,'+')
end
% 
% uq_figure
% uq_histogram(Y)