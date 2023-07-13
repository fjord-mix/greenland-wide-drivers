%% Running for all fjords of interest
if exist('fjord_model','Var'), clear fjord_model; end
fjord_model(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
for id=1:length(fjord_ids)
    fjord_name = fjord_names{id};    
    fprintf('Running model for %s... ',fjord_name)
    [p,~] = get_model_default_parameters();
    if id==2 || id==4 % we remove the sill for Sermilik and Ilulissat
        p.sill=0;
    end
    % p.N=2; % might be worth testing the model using 2 layers only
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p);
    fjord_run.p.plot_runtime=0;
    tic
    fjord_run.s = boxmodel(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = fjord_name;
    fjord_model(id) = fjord_run;
end
% save([outs_path,'example_benchmark_fjords_',num2str(p.N),'layers'],'fjord_model','-v7.3') % v7.3 allows files > 2GB
% for id=1:2, plot_outputs(fjord_model(id)); end %plot_fluxes(fjord_model(id)); end

%% Loop to plot model outputs for desired fjords
% fjord_model = load([outs_path,output_fname]).fjord_model;
for id=1:length(fjord_model)
    fprintf('Plotting model outputs for %s... \n',fjord_model(id).m.name)   
    % plot_outputs(fjord_model(id))
    animate(fjord_model(id),[figs_path,'/animations/'],fjord_model(id).m.name,200)
end