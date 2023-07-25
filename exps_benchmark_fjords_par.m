%% Running for all fjords of interest
% starts parallel pool
num_workers=2;
if ~exist('poolobj','var')
    poolobj =  parpool(num_workers); 
    addAttachedFiles(gcp,regexp(genpath(model_path), ':', 'split'));
    addAttachedFiles(gcp,regexp(genpath(collation_path), ':', 'split'));
    addAttachedFiles(gcp,regexp(genpath(import_path), ':', 'split'));
end
if exist('fjord_model','Var'), clear fjord_model; end
fjord_model(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
% fjord_model(length(fjord_ids)) = struct;

% set up runs
for id=1:length(fjord_ids)
    fjord_name = fjord_names{id};        
    [p,~] = get_model_default_parameters();
    if id==2 || id==4 % we remove the sill for Sermilik and Ilulissat
        p.sill=0;
    end
    % p.N=2; % might be worth testing the model using 2 layers only
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p);
    fjord_run.s = [];
    fjord_run.p.plot_runtime=0;    
    fjord_run.m.name = fjord_name;
    fjord_model(id) = fjord_run;
end

% run in parallel
parfor I=1:length(fjord_model)
    fprintf('Running model for %s... \n',fjord_model(I).m.name)
    fjord_model(I).s = boxmodel(fjord_model(I).p,fjord_model(I).t,fjord_model(I).f,fjord_model(I).a); % runs and gets the results    
end
%% Loop to plot model outputs for desired fjords
% fjord_model = load([exp_out_path,'example_benchmark_fjords_',num2str(p.N),'layers']).fjord_model;
% for id=1:length(fjord_model)
%     fprintf('Plotting model outputs for %s... \n',fjord_model(id).m.name)   
%     % plot_outputs(fjord_model(id))
%     animate(fjord_model(id),[figs_path,'/animations/'],fjord_model(id).m.name,200)
% end