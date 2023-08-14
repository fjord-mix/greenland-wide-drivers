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
[Parameters,IOpts,probs] = define_model_param_distrib(datasets,fjords_compilation,1);

rng('default')
status=1;
while status==1
try
    X = zeros(size(probs));
    for i=1:length(probs)
        X(i) = random(probs(i));
    end    
    [ohc,osc,status] = wrapper_boxmodel(X,Parameters);
catch ME
    disp('Unsuccessful run')
end
end

% check of anomaly PDFs
for i_reg=1:7
    [Parameters,IOpts,probs] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
    figure; hold on; 
    histogram(random(probs(8),1000)); 
end

%% Setting up the runs per region - UNTESTED

% for i_reg=1:1
%     [Parameters,IOpts,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
% 
%     MOpts.mFile = 'wrapper_boxmodel';
%     MOpts.isVectorized = false;
% 
%     % Specification of 14th degree LARS−based PCE     
%     MOpts.Type = 'Metamodel';
%     MOpts.MetaType = 'PCE';    
%     MOpts.Method = 'LARS';
%     MOpts.Degree = 14;
%     MOpts.ExpDesign.NSamples = 1000;
%     model = uq_createModel(MOpts);
% 
%     input = uq_getInput(IOpts);
%     X = uq_getSample(input,500,'LHS');
% 
%     Y = uq_evalModel(model,X);
% end
% 
% uq_figure
% uq_histogram(Y)