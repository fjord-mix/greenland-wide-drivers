%% Running for all fjords of interest
% name_exp = '';
% output_fname=['test_fjords_',name_exp];
if exist('fjord_model','Var'), clear fjord_model; end
fjord_model(length(fjord_ids)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);
for id=1:length(fjord_ids)
    fjord_name = fjord_names{id};    
    fprintf('Running model for %s... ',fjord_name)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = fjord_name;   
    fjord_model(id) = fjord_run;
    % plot_outputs(fjord_run)
end
% save([outs_path,output_fname],'fjord_model','-v7.3') % v7.3 allows files > 2GB
% for id=1:2, plot_outputs(fjord_model(id)); end %plot_fluxes(fjord_model(id)); end

for id=1:2, plot_outputs(fjord_model(id)); end %plot_fluxes(fjord_model(id)); end

%% Loop to plot model outputs for desired fjords
fjord_model = load([outs_path,output_fname]).fjord_model;
for id=1:length(fjord_model)
    % fprintf('Plotting model outputs for %s... ',fjord_model(id).m.name)   
    % plot_outputs(fjord_run)
    % plot_figures(datasets,fjord_model(id),fjords_map, datasets.obs.ctd_data.omg)
    % disp('Done.')
    if ~fjord_model(id).s.status
        disp('Creating video...')
        animate_v4p1(fjord_model(id),figs_path,[output_fname,'_',fjord_model(id).m.name],50);
        disp('Done.')
    end
end