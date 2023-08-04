%% Experiment(s) to check whether we can fix the initialisation problems
% 1. bad interpolation from shelf profile
% 2. nudging using runtime profile instead of fixed values
id=1;
datasets.opts.dt = 0.5;
name_ctrl = sprintf('%s_C%1.e_P%d_dt%0.2f_tau365d',fjord_keys{id},1e4,25,datasets.opts.dt);
fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
fjord_run.p.trelax = 365*86400;
fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
fjord_run.m.name = name_ctrl;
plot_outputs([],fjord_run);

plot_initial_conds([],fjord_run);
plot_fluxes(fjord_run);
animate_v4p1(fjord_run,figs_path,[name_ctrl,'_havg'],50);