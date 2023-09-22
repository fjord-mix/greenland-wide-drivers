%% Sets up to run the box model within the UQLab framework

%% Configuring paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

[datasets,fjords_compilation,~,~] = compile_datasets(data_path);
outs_path = [data_path,'/greenland/FjordMIX/boxmodel/pce/']; % where the model output files will be saved
figs_path = [project_path,'/figs/pce/'];                     % where the figures and animations will be saved
letters = {'a','b','c','d','e','f','g','h'};

%% Showing all input parameters first
% hf = plot_distributions(datasets,fjords_compilation);
% exportgraphics(hf,[figs_path,'summary_input_params.png'],'Resolution',300)

%% Initialise all needed variables
regions = {'SW','SE','CW','CE','NW','NE','NO'};
n_runs = 50;
n_regions = length(regions);

% X       = zeros([n_runs,10,n_regions]);
ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
ohc_pd  = cell([1,n_regions]);
osc_pd  = cell([1,n_regions]);
ohc_ks  = cell([1, n_regions]);
osc_ks  = cell([1, n_regions]);
Parameters = cell([1, n_regions]);

X       = zeros([n_runs,n_regions,10]);
% ohc_out = zeros([1,n_runs*n_regions]);
% osc_out = zeros([1,n_runs*n_regions]);

% Initialise UQLab
uqlab

%% sample the distributions for inputs
rng('default') % set the seed for reproducibility
for i_reg=1:n_regions
    [Parameters{i_reg},IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg); % available outputs: [Parameters,IOpts,probs,fjords_processed]
    input = uq_createInput(IOpts);
    X(:,i_reg,:) = uq_getSample(input,n_runs,'LHS'); % perform latin hypercube sampling
end

%% Run the model for all iterations
% num_workers=4;
% checkPool = gcp('nocreate'); % If no pool, do not create one
% if isempty(checkPool) % if there is no pool
%     parpool(num_workers);
% end

% Using a parfor inside the for loop will likely increase overhead time, 
% but will significantly reduce the amount of memory used used 
% by slicing X and Prameters before the parallel runs

for i_reg=1:n_regions
    Xreg = squeeze(X(:,i_reg,:));
    Params_reg = Parameters{i_reg};
    tic
    for k_run=1:n_runs
        try
            % [ensemble(k).time,ensemble(k).ohc,ensemble(k).osc,status,fjord_run(k)] = wrapper_boxmodel(X,Parameters);
            [ohc_out(k_run,i_reg),osc_out(k_run,i_reg)] = wrapper_boxmodel(Xreg(k_run,:),Params_reg);
        catch ME
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
        fprintf('run %d complete\n',k_run)
    end
    toc
    fprintf('region %d complete\n',i_reg)
end
% reshape input/output matrices afterwards to a more tractable format
% ohc_rsp = reshape(ohc_out,n_runs,n_regions);
% osc_rsp = reshape(osc_out,n_runs,n_regions);
% Xrsp    = reshape(X,n_runs,n_regions,10);
%% Calculate the distributions based on the numerical outputs alone

for i_reg=1:n_regions
    ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
    osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
    ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
    osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
end
% save outputs in case matlab crashes, or so we dont have to re-run it
save([outs_path,'ohc_osc_runs_probs_n30'],'ohc_out','osc_out','ohc_pd','osc_pd','ohc_ks','osc_ks')
% load([outs_path,'ohc_osc_runs_probs_n30'])
%% Setting up the PCE model per region using UQLab

Ysur_ohc      = cell([1,n_regions]); % surrogate model results evaluated at the same points as the numerical model
Ysur_osc      = cell([1,n_regions]);
Ynum_ohc      = cell([1,n_regions]); % numerical model results
Ynum_osc      = cell([1,n_regions]);
Yeval_ohc     = cell([1,n_regions]); % surrogate model results evaluated at a much larger input range
Yeval_osc     = cell([1,n_regions]);
sobolA_ohc    = cell([1,n_regions]);
sobolA_osc    = cell([1,n_regions]);
sur_model_ohc = cell([1,n_regions]);
sur_model_osc = cell([1,n_regions]);
ohc_ks_eval   = cell([1,n_regions]);
osc_ks_eval   = cell([1,n_regions]);

SobolOpts.Type             = 'Sensitivity';
SobolOpts.Method           = 'Sobol';
SobolOpts.Sobol.Order      = 1;
SobolSensOpts.SaveEvaluations = false; % to prevent excessive memory usage!
% SobolOpts.Sobol.SampleSize = 1e5; % only needed for MC, but we are using the PCE model itself

% Initialise UQLab
uqlab
for i_reg=1:n_regions
    [Params_reg,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);

    % create numerical model object
    ModelOpts.mFile = 'wrapper_boxmodel';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters=Params_reg;
    num_model = uq_createModel(ModelOpts);
    
    % create inputs object for the surrogate model
    input = uq_createInput(IOpts);
    Xeval = uq_getSample(input,1e6,'LHS');

    % create Meta (surrogate) Model: specification of 14th degree LARS−based PCE     
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';    
    MetaOpts.FullModel = num_model;
    MetaOpts.Method = 'LARS';
    MetaOpts.Degree = 14;

    % Specifying NSamples would get UQLab to perform the runs for us
    % MetaOpts.ExpDesign.NSamples = n_runs;

    % Instead we ran the simulations before to be able to ignore unstable/crashed runs
    Xreg = X(:,i_reg,:);                    
    ohc_reg = ohc_out(:,i_reg);
    osc_reg = osc_out(:,i_reg);
    Ynum_ohc{i_reg} = ohc_reg;
    Ynum_osc{i_reg} = osc_reg;

    MetaOpts.ExpDesign.X = Xreg(~isnan(ohc_reg),:);
    MetaOpts.ExpDesign.Y = ohc_reg(~isnan(ohc_reg));

    % Create the surrogate model and evaluate it
    sur_model_ohc{i_reg} = uq_createModel(MetaOpts);
    Ysur_ohc{i_reg}      = uq_evalModel(sur_model_ohc{i_reg},Xreg);  % run the surrogate model for the same inputs as the numerical model
    Yeval_ohc{i_reg}     = uq_evalModel(sur_model_ohc{i_reg},Xeval); % run the surrogate model for a much larger N
    
    sobolA_ohc{i_reg}  = uq_createAnalysis(SobolOpts);   % compute Sobol indices based on the last model run
    ohc_ks_eval{i_reg} = fitdist(Yeval_ohc{i_reg},'kernel'); % create a kernel density function from the surrogate model outputs

    % Same for salt content - but no need to change inputs because they are the same
    MetaOpts.ExpDesign.Y = osc_reg(~isnan(osc_reg)); 
    sur_model_osc{i_reg} = uq_createModel(MetaOpts);
    Ysur_osc{i_reg}      = uq_evalModel(sur_model_osc{i_reg},Xreg); 
    Yeval_osc{i_reg}     = uq_evalModel(sur_model_osc{i_reg},Xeval);
    
    sobolA_osc{i_reg}  = uq_createAnalysis(SobolOpts);
    osc_ks_eval{i_reg} = fitdist(Yeval_osc{i_reg},'kernel');
    % uq_histogram(Yeval)
end
% save([outs_path,'ohc_osc_pce_n50_n1e6'],'sur_model_ohc','sur_model_osc','Ysur_ohc','Ysur_osc','Yeval_ohc','Yeval_osc','sobolA_ohc','sobolA_osc')

%% Plotting the outputs
run make_figs_uqlab_outputs.m