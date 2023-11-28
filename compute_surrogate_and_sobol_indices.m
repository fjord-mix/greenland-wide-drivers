Ynum_ohc      = cell([1,n_regions]); % Ynum: numerical model results for training
Ynum_osc      = cell([1,n_regions]);
Yind_ohc      = cell([1,n_regions]); % Yind: numerical model results for validation
Yind_osc      = cell([1,n_regions]);
Ysur_ohc      = cell([1,n_regions]); % Ysur: surrogate model results evaluated at the same points
Ysur_osc      = cell([1,n_regions]); 
Yeval_ohc     = cell([1,n_regions]); % Yeval: surrogate model results evaluated at a much larger input range
Yeval_osc     = cell([1,n_regions]);
Yvld_ohc      = cell([1,n_regions]); % Yvld: surrogate model results for independent dataset (validation)
Yvld_osc      = cell([1,n_regions]);
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
SobolSensOpts.SaveEvaluations = false; % to prevent excessive memory usage!
% SobolOpts.Sobol.SampleSize = 1e5; % only needed for MC, but we are using the PCE model itself


% Initialise UQLab - commented out because we initialised it earlier
% uqlab % uncomment if loading model results run in a previous Matlab session
for i_reg=1:n_regions
    % [Params_reg,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
    Params_reg = Parameters{i_reg};
    IOpts_reg = IOpts{i_reg};

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
    MetaOpts.Method = 'OMP'; % 'LARS','OMP', or 'SP'
    MetaOpts.Degree = 40:1:80;
    MetaOpts.TruncOptions.qNorm = 0.1:0.1:1;
    MetaOpts.TruncOptions.MaxInteraction = 2;
    % MetaOpts.Bootstrap.Replications = 1e2; % for assessing the model accuracy
    % MetaOpts.LARS.LarsEarlyStop = false;
    MetaOpts.OMP.OmpEarlyStop = false;
    MetaOpts.OMP.ModifiedLoo=0;

    % Specifying NSamples would get UQLab to perform the runs for us
    % we do not do that here because we need to exclude the unstable runs
    % from the analysis
    % MetaOpts.ExpDesign.NSamples = n_runs;

    % So instead we ran the simulations before, and just filtered the runs
    % with NaN as a result (as defined in the wrapper_boxmodel function)
    Xreg = X(:,i_reg,:);                    
    Xind = Xvalid(:,i_reg,:);
    Xsur = Xeval(:,i_reg,:);
    ohc_reg = ohc_out(:,i_reg);
    osc_reg = osc_out(:,i_reg);
    ohc_ind = ohc_vld(:,i_reg);
    osc_ind = osc_vld(:,i_reg);
    Ynum_ohc{i_reg} = ohc_reg;
    Ynum_osc{i_reg} = osc_reg;
    Yind_ohc{i_reg} = ohc_ind;
    Yind_osc{i_reg} = osc_ind;

    MetaOpts.ExpDesign.X = single(Xreg(~isnan(ohc_reg),:));
    MetaOpts.ExpDesign.Y = single(ohc_reg(~isnan(ohc_reg)));
    MetaOpts.ValidationSet.X = Xind(~isnan(ohc_ind),:);
    MetaOpts.ValidationSet.Y = ohc_ind(~isnan(ohc_ind));


    fprintf('Creating OHC surrogate model for %s...\n',regions{i_reg})
    sur_model_ohc{i_reg} = uq_createModel(MetaOpts);                 % Create the surrogate model
    fprintf('Done. Evaluating model and storing results...\n')
    % [Ysur_ohc{i_reg},~,Yboo_ohc{i_reg}]      = uq_evalModel(sur_model_ohc{i_reg},Xreg);  % run the surrogate model for the same inputs as the numerical model (incl. bootstraping results)
    Yvld_ohc{i_reg}  = uq_evalModel(sur_model_ohc{i_reg},Xind);  % run the surrogate model for an independent set of inputs (validation)
    Ysur_ohc{i_reg}  = uq_evalModel(sur_model_ohc{i_reg},Xreg); % run the surrogate model for the same inputs as the numerical model
    Yeval_ohc{i_reg} = uq_evalModel(sur_model_ohc{i_reg},Xsur); % run the surrogate model for a much larger N
    
    fprintf('Done. Computing OHC Sobol indices for %s...\n',regions{i_reg})
    sobolA_ohc{i_reg}  = uq_createAnalysis(SobolOpts);   % compute Sobol indices based on the last model run
    ohc_ks_eval{i_reg} = fitdist(Yeval_ohc{i_reg},'kernel'); % create a kernel density function from the surrogate model outputs

    % Same for salt content - but no need to change inputs because they are the same
    fprintf('Creating OSC surrogate model for %s...\n',regions{i_reg})
    MetaOpts.Degree = 10:1:20;
    MetaOpts.ExpDesign.Y = osc_reg(~isnan(osc_reg)); 
    MetaOpts.ValidationSet.Y = osc_ind(~isnan(osc_ind));
    sur_model_osc{i_reg} = uq_createModel(MetaOpts);
    fprintf('Done. Evaluating model and storing results...\n')
    Ysur_osc{i_reg}  = uq_evalModel(sur_model_osc{i_reg},Xreg); 
    Yvld_osc{i_reg}  = uq_evalModel(sur_model_osc{i_reg},Xind);
    Yeval_osc{i_reg} = uq_evalModel(sur_model_osc{i_reg},Xsur);
    
    fprintf('Done. Computing OSC Sobol indices for %s...\n',regions{i_reg})
    sobolA_osc{i_reg}  = uq_createAnalysis(SobolOpts);
    osc_ks_eval{i_reg} = fitdist(Yeval_osc{i_reg},'kernel');

end
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
fprintf('Done with convergence test for n_runs.\n')