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

%% evaluate the wrapper by itself first
rng('default') % set the seed for reproducibility
regions = {'SW','SE','CW','CE','NW','NE','NO'};
n_runs = 50;
n_regions = length(regions);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
ohc_pd  = cell([1,n_regions]);
osc_pd  = cell([1,n_regions]);
ohc_ks  = cell([1, n_regions]);
osc_ks  = cell([1, n_regions]);
X       = zeros([n_runs,10,n_regions]);

% Initialise UQLab
uqlab
for i_reg=1:n_regions
    [Parameters,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg); % available outputs: [Parameters,IOpts,probs,fjords_processed]
    input = uq_createInput(IOpts);
    X(:,:,i_reg) = uq_getSample(input,n_runs,'LHS'); % perform latin hypercube sampling
    for k=1:n_runs
        try
            tic
            % [ensemble(k).time,ensemble(k).ohc,ensemble(k).osc,status,fjord_run(k)] = wrapper_boxmodel(X,Parameters);
            [ohc_out(k,i_reg),osc_out(k,i_reg)] = wrapper_boxmodel(X(k,:,i_reg),Parameters);
            toc            
        catch ME
            fprintf('Crash on iteration #%d\n',k)
            ohc_out(k,i_reg) = NaN;
            osc_out(k,i_reg) = NaN;
        end
    end
    ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
    osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
    ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
    osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
end

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
    [Parameters,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);

    % create numerical model object
    ModelOpts.mFile = 'wrapper_boxmodel';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters=Parameters;
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
    Xreg = X(:,:,i_reg);                    
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

%% Plotting the results of the numerical model alone

% we want to know how many runs were successful
ok_runs = zeros([1, n_regions]);
regions_lbl = regions;
for i_reg=1:n_regions
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions_lbl{i_reg} = [regions_lbl{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];
end

% construct the plot itself
ohc_x = linspace(0.99*min(ohc_out(:)),1.01*max(ohc_out(:)),1000);
osc_x = linspace(0.99*min(osc_out(:)),1.01*max(osc_out(:)),1000);
figure('Position',[40 40 850 300]); 
subplot(1,2,1), hold on; 
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_pd{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content (J m^{-3})',fontsize=14); ylabel('Probability',fontsize=14);  box on
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
subplot(1,2,2), hold on; 
for i_reg=1:n_regions,plot(osc_x,pdf(osc_pd{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content (g m^{-3})',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
hl = legend(regions_lbl,'fontsize',14,'Location','north');
% exportgraphics(gcf,[figs_path,'test_output_ohc_osc_pdf_n50.png'],'Resolution',300)


% figure; hold on; for k=1:n_runs, plot(ensemble(k).time,ensemble(k).ohc); end

%% Plotting surrogate vs numerical model
figure('Name','Model fit','position',[40 40 1000 400])
for i_reg=1:n_regions
    subplot(2,4,i_reg)
    uq_plot(gca,Ynum_ohc{i_reg},Ysur_ohc{i_reg},'+')
    text(0.05,0.95,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
    set(gca,'fontsize',14,'XTickLabel',[],'YTickLabel',[])
    if i_reg > 3, xlabel('Numerical model'); end
    if ismember(i_reg,[1,5]), ylabel('Surrogate model'); end
end
% exportgraphics(gcf,[figs_path,'pce_fit_indices_n50.png'],'Resolution',300)

%% Plotting the surrogate model kernel density
figure('Position',[40 40 1000 400]); hold on;
subplot(1,2,1), hold on; box on
for i_reg=1:n_regions,plot(ohc_x,pdf(ohc_ks_eval{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content (J m^{-3})',fontsize=14); ylabel('Kernel density',fontsize=14);  
set(gca,'fontsize',14)
subplot(1,2,2), hold on; box on
for i_reg=1:n_regions,plot(osc_x,pdf(osc_ks_eval{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content (g m^{-3})',fontsize=14); 
set(gca,'fontsize',14)
hl = legend(regions,'fontsize',14,'Location','northeast');
% exportgraphics(gcf,[figs_path,'test_pce_ks_n1e6.png'],'Resolution',300)


%% Plotting the Sobol indices
sobolTotal_ohc = [];
sobolFirstOrder_ohc = [];

% Gather up all regions here
for i_reg=1:7
    sobolResults_ohc  = sobolA_ohc{i_reg}.Results;
    % uq_print(sobolAnalysis)
    sobolTotal_ohc      = [sobolTotal_ohc sobolResults_ohc.Total];
    sobolFirstOrder_ohc = [sobolFirstOrder_ohc sobolResults_ohc.FirstOrder];

    sobolResults_osc  = sobolA_osc{i_reg}.Results;    
    sobolTotal_osc      = [sobolTotal_osc sobolResults_osc.Total];
    sobolFirstOrder_osc = [sobolFirstOrder_osc sobolResults_osc.FirstOrder];
end

% Plotting indices
uq_figure('Name', 'Total Sobol'' Indices','Position',[40 40 1000 500])
barWidth = 0.5;
subplot(1,2,1)
uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_ohc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('Total Sobol indices (OHC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)

subplot(1,2,2)
uq_bar(gca,1:length(IOpts.Marginals), sobolTotal_osc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('Total Sobol indices (OSC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
uq_legend(regions,'Location', 'northeastoutside')
% exportgraphics(gcf,[figs_path,'test_total_sobol_indices_n50.png'],'Resolution',300)

uq_figure('Name', 'First-Order Sobol'' Indices','Position',[40 40 1000 500])
subplot(1,2,1)
uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_ohc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OHC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_ohc.VariableNames)
subplot(1,2,2)
uq_bar(gca,1:length(IOpts.Marginals), sobolFirstOrder_osc, barWidth)
xlim([0 length(IOpts.Marginals)+1]); % ylim([0 1]); 
xlabel('Variable'); ylabel('First-Order Sobol'' Indices (OSC)')
set(gca,'XTick', 1:length(IOpts.Marginals),'XTickLabel', sobolResults_osc.VariableNames)
uq_legend(regions,'Location', 'northeastoutside')
% exportgraphics(gcf,[figs_path,'test_1st_sobol_indices_n50.png'],'Resolution',300)