%% Sensitivity to iceberg melt efficiency, i.e., varying a single parametre
id=1;
[p_def,~] = get_model_default_parameters();
if id==2 || id==4 % we remove the sill for Sermilik and Ilulissat
    p_def.sill=0;
end

M0_space = [2e-8,2e-10,2e-9,2e-7,2e-6]; % Defining the parametre space (default M0: 2e-8)

% Iterating over a loop performing one run for each combination
for i_param=1:length(M0_space)
    name_exp = sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(i_param));    
    fprintf('Running model for %s... ',name_exp)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p_def);
    fjord_run.p.M0=M0_space(i_param);
    tic
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.o = postprocess_boxmodel(fjord_run);
    fjord_run.m.name = name_exp;   
    save([exp_out_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB; MPhr stands for the ocean forcing used    
end

%% Plotting

% This is a useful way to evaluate the differences between a set of runs
name_set = sprintf('%s_M0*',fjord_keys{id});
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of icebergs','Set2')
exportgraphics(gcf,[figs_path,'/test_iceberg/sens_icebergs.png'],'ContentType','vector','BackgroundColor','none'); close gcf

% Individual plots for all simulations performed
chosen_depths=[5, 200, 500];
for i_param=1:length(M0_space)
    name_exp = sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(i_param));    
    fjord_run = load([exp_out_path,name_exp,'.mat']).fjord_run;
    hf=plot_ts_at_depth(fjord_run,chosen_depths,'nearest');
    % exportgraphics(gcf,[figs_path,'/',name_set,'.pdf'],'ContentType','vector','BackgroundColor','none'); close gcf;
    % animate(fjord_run,figs_path,name_exp,50);
end
