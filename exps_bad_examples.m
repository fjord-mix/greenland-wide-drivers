if exist('fjords_bad','Var'), clear fjords_bad; end
fjords_bad(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
for id=1:1%length(fjord_ids)
    fjord_name = fjord_names{id};    
    fprintf('Running model for %s... ',fjord_name)    
    [p,~] = get_model_default_parameters();        
    fjord_run = prepare_boxmodel_input(datasets_mod,fjords_processed(fjord_ids(id)),fjord_ids(id),p);
    fjord_run.p.trelax=5; % exacerbates numerical instability in the mixing fluxes
    fjord_run.p.plot_runtime=0;
    tic
    fjord_run.s = boxmodel(fjord_run.p,fjord_run.t,fjord_run.f,fjord_run.a);
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = fjord_name;
    fjords_bad(id) = fjord_run;
end

%% Loop to plot model outputs for desired fjords
% fjord_model = load([exp_out_path,'example_benchmark_fjords_',num2str(p.N),'layers']).fjord_model;
for id=1:length(fjords_bad)
    fprintf('Plotting model outputs for %s... \n',fjords_bad(id).m.name)   
    plot_outputs(fjord_model(id))
    % animate(fjords_bad(id),[figs_path,'/animations/bad_'],fjords_bad(id).m.name,200)
end

%% Other options to take a closer look at specific fjords
plot_TS_boxmodel(fjords_bad(3))%,fjords_map,datasets.obs.ctd_data.all);
plot_fluxes(fjords_bad(3))
plot_outputs(fjords_bad(3))