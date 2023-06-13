%% Experiment(s) to check whether we can fix the mixing fluxes problems
% Buoyancy goes to negative values, resulting in Re < 0 and thus complex values
dt_space = [1, 6, 12, 24]/24.; % time step - dt = [0.5,12] hours
K0_space = [5e-1, 5e-2, 5e-3, 5e-10]; % mixing coefficient (default 0.05)

id=1;
% i_dt=3;
% i_K0=1;
for i_K0=1:length(K0_space)
for i_dt=1:length(dt_space)
    datasets.opts.dt = dt_space(i_dt);
    name_exp = sprintf('%s_K%1.e_dt%0.2f',fjord_keys{id},K0_space(i_K0),dt_space(i_dt));
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    fjord_run.p.trelax = 365*86400;
    fjord_run.p.K0=K0_space(i_K0);
    fprintf('Running model for %s... ',name_exp)
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = name_exp;
    save([outs_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
    % plot_outputs([],fjord_run);
end
end
name_set_K0 = sprintf('%s_K*_dt%0.2f',fjord_keys{id},dt_space(i_dt));
plot_subset_runs([outs_path,name_set_K0,'.mat'],'Effect of K0')

name_set_dt = sprintf('%s_K%1.e_dt*',fjord_keys{id},K0_space(1));
plot_subset_runs([outs_path,name_set_dt,'.mat'],'Effect of dt')

name_set = sprintf('%s_K*_dt*',fjord_keys{id});
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of dt and K0')

%% Control run
name_ctrl = sprintf('%s_K%1.e_dt%0.2f',fjord_keys{id},K0_space(1),dt_space(4));
datasets.opts.dt = dt_space(4);
fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
fjord_run.p.trelax = 365*86400;
fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
fjord_run.m.name = name_ctrl;
plot_outputs([],fjord_run);


plot_initial_conds([],fjord_run);
plot_fluxes(fjord_run);
animate_v4p1(fjord_run,figs_path,[name_ctrl,'_mixing'],50);

%% Testing if it holds for Sermilik, because it is too good to be true (of course it did not)
id=2;
for i_K0=1:length(K0_space)
% for i_dt=1:length(dt_space)
    datasets.opts.dt = dt_space(i_dt);
    name_exp = sprintf('%s_K%1.e_dt%0.2f',fjord_keys{id},K0_space(i_K0),dt_space(i_dt));
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    fjord_run.p.trelax = 365*86400;
    fjord_run.p.K0=K0_space(i_K0);
    fprintf('Running model for %s... ',name_exp)
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = name_exp;
    save([outs_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
    % plot_outputs([],fjord_run);
% end
end
name_set_K0 = sprintf('%s_K*_dt%0.2f',fjord_keys{id},dt_space(i_dt));
plot_subset_runs([outs_path,name_set_K0,'.mat'],'Effect of K0 on Sermilik')


name_ctrl = sprintf('%s_K%1.e_dt%0.2f',fjord_keys{id},K0_space(1),dt_space(4));
datasets.opts.dt = dt_space(4);
fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
fjord_run.p.trelax = 365*86400;
fjord_run.p.K0=K0_space(2);
fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
fjord_run.m.name = name_ctrl;
% save([outs_path,name_ctrl,'','.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
plot_outputs([],fjord_run);

[~]=plot_TS_boxmodel(datasets,fjords_map,name_ctrl,fjord_run,fjord_run.s.t,datasets.obs.ctd_data.all);

plot_fluxes(fjord_run)

plot_initial_conds([],fjord_run);
plot_fluxes(fjord_run);
animate_v4p1(fjord_run,figs_path,[name_ctrl,'_mixing'],50);