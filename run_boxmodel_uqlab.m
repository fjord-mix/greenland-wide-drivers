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

%% Showing all input parameters first
% hs = plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off'; 
% hf = plot_distributions(datasets,fjords_compilation);
% exportgraphics(hf,[figs_path,'summary_input_probs2010-2018.png'],'Resolution',300)

%% Initialise all needed variables
regions = {'SW','SE','CW','CE','NW','NE','NO'};
n_runs    = 200;
time_step = 0.1; % in days
n_regions = length(regions);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
ensemble(n_runs,n_regions) = struct("time",[],"ohc",[],"osc",[]);
ohc_pd  = cell([1,n_regions]);
osc_pd  = cell([1,n_regions]);
ohc_ks  = cell([1, n_regions]);
osc_ks  = cell([1, n_regions]);


Parameters = cell([1, n_regions]);
X          = zeros([n_runs,n_regions,10]);

% Initialise UQLab
uqlab

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

% save ensemble structure so we do not need to rerun it all the time
save([outs_path,'hc_sc_ensemble_n',num2str(n_runs)],'-v7.3','ensemble')

%% Calculate the distributions based on the numerical outputs alone
% load([outs_path,'ohc_osc_runs_probs_n',num2str(n_runs)]) % if we have the results saved already

for i_reg=1:n_regions
    ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
    osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
    ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
    osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
end
% save outputs so we dont have to re-run it
save([outs_path,'ohc_osc_change_runs_probs_n',num2str(n_runs)],...
      'ohc_out','osc_out','ohc_pd','osc_pd','ohc_ks','osc_ks')
%% Setting up the PCE model per region using UQLab

Ynum_ohc      = cell([1,n_regions]); % Ynum: numerical model results
Ynum_osc      = cell([1,n_regions]);
Ysur_ohc      = cell([1,n_regions]); % Ysur: surrogate model results evaluated at the same points
Ysur_osc      = cell([1,n_regions]); % as the numerical model, used to evaluate model fit
Yeval_ohc     = cell([1,n_regions]); % Yeval: surrogate model results evaluated at a much larger input range
Yeval_osc     = cell([1,n_regions]);
sur_model_ohc = cell([1,n_regions]); % sur_model: UQLab model objects
sur_model_osc = cell([1,n_regions]);
ohc_ks_eval   = cell([1,n_regions]); % ks_eval: kernel density functions for Yeval
osc_ks_eval   = cell([1,n_regions]);

sobolA_ohc    = cell([1,n_regions]); % sobolA: Sobol indices analysis structures
sobolA_osc    = cell([1,n_regions]);
SobolOpts.Type             = 'Sensitivity';
SobolOpts.Method           = 'Sobol';
SobolOpts.Sobol.Order      = 1;
SobolSensOpts.SaveEvaluations = false; % to prevent excessive memory usage!
% SobolOpts.Sobol.SampleSize = 1e5; % only needed for MC, but we are using the PCE model itself

% ANCOVA analysis was not as successful as originally envisioned...
% ANCOVA_ohc    = cell([1,n_regions]);
% ANCOVA_osc    = cell([1,n_regions]);
% ANCOVASensOpts_ohc.Type = 'Sensitivity';
% ANCOVASensOpts_ohc.Method = 'ANCOVA';
% ANCOVASensOpts_osc = ANCOVASensOpts_ohc;
% ANCOVASensOpts.ANCOVA.SampleSize = 150;
% ANCOVAAnalysis = uq_createAnalysis(ANCOVASensOpts);


% Initialise UQLab - commented out because we initialised it earlier
% uqlab % uncomment if loading model results run in a previous Matlab session
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
    % we do not do that here because we need to exclude the unstable runs
    % from the analysis
    % MetaOpts.ExpDesign.NSamples = n_runs;

    % So instead we ran the simulations before, and just filtered the runs
    % with NaN as a result (as defined in the wrapper_boxmodel function)
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

    % ANCOVASensOpts_ohc.ANCOVA.CustomPCE = sur_model_ohc{i_reg};
    % ANCOVASensOpts_osc.ANCOVA.CustomPCE = sur_model_osc{i_reg};
    % ANCOVA_ohc{i_reg}  = uq_createAnalysis(ANCOVASensOpts_ohc);
    % ANCOVA_osc{i_reg}  = uq_createAnalysis(ANCOVASensOpts_osc);
    % uq_histogram(Yeval)
end
% ideally we would save these variables as well, but Matlab always returns an error
% when trying to save them
% save([outs_path,'ohc_osc_pce_n50_n1e6'],'sur_model_ohc','sur_model_osc','Ysur_ohc','Ysur_osc','Yeval_ohc','Yeval_osc','sobolA_ohc','sobolA_osc')

%% Plotting the outputs and an input summary
close all
plot_reg_ocn_forcings(datasets,fjords_compilation)
figure(1); exportgraphics(gcf,[figs_path,'hovmoller_Ts.png'],'Resolution',300)
figure(2); exportgraphics(gcf,[figs_path,'hovmoller_Ss.png'],'Resolution',300)
figure(3); exportgraphics(gcf,[figs_path,'series_discharge_hc_sc.png'],'Resolution',300)

run make_figs_uqlab_outputs.m