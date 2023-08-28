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

%% Showing all input parameters first
% hf = plot_distributions(datasets,fjords_compilation);
% exportgraphics(hf,[figs_path,'summary_input_params.png'],'Resolution',300)

%% evaluate the wrapper by itself first
[Parameters,IOpts,probs,fjords_processed] = define_model_param_distrib(datasets,fjords_compilation,1);

rng('default')
n_runs = 30;
clear ensemble fjord_run
% ensemble(length(n_runs)) = struct('time',[],'ohc',[],'osc',[]);

ohc_pd = cell([7,1]);
ohc_out = NaN([n_runs, 7]);
osc_pd = cell([7,1]);
osc_out = NaN([n_runs, 7]);
for i_reg=1:7    
    for k=1:n_runs
        try
            X = zeros(size(probs));
            for i=1:length(probs)
                X(i) = random(probs(i));
            end    
            tic
            % [ensemble(k).time,ensemble(k).ohc,ensemble(k).osc,status,fjord_run(k)] = wrapper_boxmodel(X,Parameters);
            [ohc_out(k,i_reg),osc_out(k,i_reg)] = wrapper_boxmodel(X,Parameters);
            toc
            ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
            osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
        catch ME
            fprintf('Crash on iteration #%d\n',k)
            ohc_out(k,i_reg) = NaN;
            osc_out(k,i_reg) = NaN;
        end
    end
end

% Plotting the results
ok_runs = zeros([1, 7]);
regions = {'SW','SE','CW','CE','NW','NE','NO'};
for i_reg=1:7
    total_runs = ohc_out(:,i_reg);
    ok_runs(i_reg) = sum(~isnan(total_runs));
    regions{i_reg} = [regions{i_reg},' ( n=',num2str(ok_runs(i_reg)),')'];
end
ohc_x = linspace(400,1600,1000);
osc_x = linspace(3.4,4,1000);
figure('Position',[40 40 850 300]); 
subplot(1,2,1), hold on; for i_reg=1:7, plot(ohc_x,pdf(ohc_pd{i_reg},ohc_x),'linewidth',2); end
xlabel('Heat content (J m^{-3})',fontsize=14); ylabel('Probability',fontsize=14);  box on
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([800 1600]);
subplot(1,2,2), hold on; for i_reg=1:7, plot(osc_x,pdf(osc_pd{i_reg},osc_x),'linewidth',2); end
xlabel('Salt content (g m^{-3})',fontsize=14); box on
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
xlim([3.4 4]);
legend(regions,'fontsize',14)
exportgraphics(gcf,[figs_path,'test_output_ohc_osc_pdfs.png'],'Resolution',300)


% figure; hold on; for k=1:n_runs, plot(ensemble(k).time,ensemble(k).ohc); end

%% Setting up the runs per region

% Initialise UQLab
uqlab
for i_reg=1:1
    [Parameters,IOpts,~,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
    rng('default')
    % create numerical model object
    ModelOpts.mFile = 'wrapper_boxmodel';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters=Parameters;
    num_model = uq_createModel(ModelOpts);
    
    % create inputs object
    input = uq_createInput(IOpts);
    X = uq_getSample(input,50,'LHS');

    % create Meta (surrogate) Model: specification of 14th degree LARS−based PCE     
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';    
    MetaOpts.FullModel = num_model;
    MetaOpts.Method = 'LARS';
    MetaOpts.Degree = 14;
    MetaOpts.ExpDesign.NSamples = 50;
    sur_model=uq_createModel(MetaOpts);
    % uq_print(sur_model)
    % uq_display(sur_model)
    Y = uq_evalModel(sur_model,X);
    Yvalidation = uq_evalModel(num_model,X);
    uq_plot(Yvalidation,Y,'+')
end

% 
% uq_figure
% uq_histogram(Y)