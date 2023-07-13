if exist('fjords_bad','Var'), clear fjords_bad; end
fjords_bad(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
for id=1:length(fjord_ids)
    fjord_name = fjord_names{id};    
    fprintf('Running model for %s... ',fjord_name)
    [p,~] = get_model_default_parameters();
    p.N=4; % might be worth testing the model using 2 layers only
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p);
    fjord_run.p.plot_runtime=1;
    tic
    fjord_run.s = boxmodel(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = fjord_name;
    fjords_bad(id) = fjord_run;
end