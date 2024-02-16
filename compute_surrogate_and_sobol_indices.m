Ynum_ohc      = cell([1,n_regions]); % Ynum: numerical model results for training
Ynum_osc      = cell([1,n_regions]);
Ysur_ohc      = cell([1,n_regions]); % Ysur: surrogate model results evaluated at the same points
Ysur_osc      = cell([1,n_regions]); 
Yeval_ohc     = cell([1,n_regions]); % Yeval: surrogate model results evaluated at a much larger input range
Yeval_osc     = cell([1,n_regions]);
Yboo_ohc      = cell([1,n_regions]); % Yboo: bootstrap PCEs from Ynum
Yboo_osc      = cell([1,n_regions]);
sur_model_ohc = cell([1,n_regions]); % sur_model: UQLab model objects
sur_model_osc = cell([1,n_regions]);
ohc_ks_eval   = cell([1,n_regions]); % ks_eval: kernel density functions for Yeval
osc_ks_eval   = cell([1,n_regions]);

sobolA_ohc    = cell([1,n_regions]); % sobolA: Sobol indices analysis structures
sobolA_osc    = cell([1,n_regions]);
SobolOpts.Type             = 'Sensitivity';
SobolOpts.Method           = 'Sobol';
SobolOpts.Sobol.Order      = 1;
SobolOpts.SaveEvaluations = false; % to prevent excessive memory usage!
% SobolOpts.Sobol.SampleSize = 1e5; % only needed for MC, but we are using the PCE model itself


% Initialise UQLab - commented out because we initialised it earlier
% uqlab % uncomment if loading model results run in a previous Matlab session
for i_reg=1:n_regions
    % [Params_reg,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
    Params_reg = Parameters{i_reg};
    % IOpts_reg = IOpts{i_reg};
    input = uq_createInput(IOpts{i_reg}); % create probability functions

    % create numerical model object
    ModelOpts.mFile = 'wrapper_boxmodel';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters=Params_reg;
    num_model = uq_createModel(ModelOpts);

    % create Meta (surrogate) Model object
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';    
    MetaOpts.DegreeEarlyStop = false;
    MetaOpts.FullModel = num_model;
    MetaOpts.Method = 'LARS'; % 'LARS','OMP', or 'SP'
    MetaOpts.Degree = 1:1:30; % 01-30 RMSEs: {0.11,0.11,0.09,0.14,0.13,0.10,0.16} deg. C (0.84)
                              % 10-30 RMSEs: {0.11,0.11,0.09,0.14,0.13,0.10,0.15} deg. C (0.83)
                              % 30-60 RMSEs: {0.13,0.11,0.07,0.12,0.12,0.12,0.10} deg. C (0.77)*
                              % 40-80 RMSEs: {0.14,0.10,0.09,0.14,0.15,0.08,0.14} deg. C (0.84)
    MetaOpts.TruncOptions.qNorm = 0.1:0.1:1;
    MetaOpts.TruncOptions.MaxInteraction = 2;
    MetaOpts.Bootstrap.Replications = 1e3; % for assessing the model accuracy
    % MetaOpts.LARS.LarsEarlyStop = false;
    MetaOpts.OMP.OmpEarlyStop = false;
    % MetaOpts.OMP.ModifiedLoo=0;

    % Specifying NSamples would get UQLab to perform the runs for us
    % we do not do that here because we need to exclude the unstable runs
    % from the analysis
    % MetaOpts.ExpDesign.NSamples = n_runs;

    % So instead we ran the simulations before, and just filtered the runs
    % with NaN as a result (as defined in the wrapper_boxmodel function)
    Xreg = X(:,i_reg,:);                    
    Xsur = Xeval(:,i_reg,:);
    ohc_reg = ohc_out(:,i_reg);
    osc_reg = osc_out(:,i_reg);
    Ynum_ohc{i_reg} = ohc_reg;
    Ynum_osc{i_reg} = osc_reg;
    
    MetaOpts.ExpDesign.X = single(Xreg(~isnan(ohc_reg),:));
    MetaOpts.ExpDesign.Y = single(ohc_reg(~isnan(ohc_reg)));
    
    fprintf('Creating OHC surrogate model for %s...\n',regions{i_reg})
    sur_model_ohc{i_reg} = uq_createModel(MetaOpts);                 % Create the surrogate model
    fprintf('Done. Evaluating model and storing results...\n')
    [Ysur_ohc{i_reg},~,Yboo_ohc{i_reg}]      = uq_evalModel(sur_model_ohc{i_reg},Xreg);  % run the surrogate model for the same inputs as the numerical model (incl. bootstraping results)
    % Ysur_ohc{i_reg}  = uq_evalModel(sur_model_ohc{i_reg},Xreg); % run the surrogate model for the same inputs as the numerical model
    Yeval_ohc{i_reg} = uq_evalModel(sur_model_ohc{i_reg},Xsur); % run the surrogate model for a much larger N
    
    fprintf('Done. Computing OHC Sobol indices for %s...\n',regions{i_reg})
    sobolA_ohc{i_reg}  = uq_createAnalysis(SobolOpts);   % compute Sobol indices based on the last model run
    ohc_ks_eval{i_reg} = fitdist(Yeval_ohc{i_reg},'kernel'); % create a kernel density function from the surrogate model outputs

    % Same for salt content - but no need to change inputs because they are the same
    fprintf('Creating OSC surrogate model for %s...\n',regions{i_reg})
    MetaOpts.ExpDesign.Y = osc_reg(~isnan(osc_reg)); 
    sur_model_osc{i_reg} = uq_createModel(MetaOpts);
    fprintf('Done. Evaluating model and storing results...\n')
    [Ysur_osc{i_reg},~,Yboo_osc{i_reg}]      = uq_evalModel(sur_model_osc{i_reg},Xreg);
    % Ysur_osc{i_reg}  = uq_evalModel(sur_model_osc{i_reg},Xreg); 
    Yeval_osc{i_reg} = uq_evalModel(sur_model_osc{i_reg},Xsur);
    
    fprintf('Done. Computing OSC Sobol indices for %s...\n',regions{i_reg})
    sobolA_osc{i_reg}  = uq_createAnalysis(SobolOpts);
    osc_ks_eval{i_reg} = fitdist(Yeval_osc{i_reg},'kernel');

end
clear Params_reg IOpts_reg MetaOpts ModelOpts SobolOpts num_model % we need some memory management...

fprintf('Done creating PCEs and Sobol'' indices. Computing Borgonovo indices...\n')
    
run borgonovo_analysis.m
fprintf('Done computing Borgonovo indices. Starting convergence test for n_runs...\n')

%% Checking if our choice of n_runs was enough
x_subsample=10:5:n_runs;
n_subsample=length(x_subsample);
Yconv_ohc=NaN([n_subsample,n_regions]);
Yconv_osc=NaN([n_subsample,n_regions]);
for i_reg=1:n_regions
    for i_subsample=1:n_subsample
        Xsub=X(1:x_subsample(i_subsample),i_reg,:);
        ohc_conv_sum=0;
        osc_conv_sum=0;
        for i_run=1:i_subsample
            ohc_conv_sum = ohc_conv_sum + uq_evalModel(sur_model_ohc{i_reg},Xsub(i_run,:));
            osc_conv_sum = osc_conv_sum + uq_evalModel(sur_model_osc{i_reg},Xsub(i_run,:));
        end
        Yconv_ohc(i_subsample,i_reg)=ohc_conv_sum./x_subsample(i_subsample);
        Yconv_osc(i_subsample,i_reg)=osc_conv_sum./x_subsample(i_subsample);
    end
end
clear ohc_conv_sum osc_conv_sum
fprintf('Done with convergence test for n_runs.\n')