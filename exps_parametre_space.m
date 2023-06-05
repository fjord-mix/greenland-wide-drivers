%% Exploring the parametre space
% 3 parametres, 4 values each - 4^3 = 64 runs
dt_space = [1, 6, 12, 24]/24.; % time step - dt = [0.5,12] hours
C0_space = [1e3, 5e3, 1e4, 5e4]; % shelf-exchange coefficient - C0 = [1e3,1e5]
P0_space=[10, 25, 50, 100]; % plume-entrainment coefficient - P0 = [10,25,50,100] m
wm_space=[3e-5,2e-5,1e-5,1e-6]; % maximum vertical mixing velocity (default: 4e-5)
K0_space=[5e-10, 5e-5, 5e-3]; % mixing coefficient (default 0.05)
% tau_space = [10, 30, 90, 365, 3650]*86400; % relaxation time - trelax = [10, 10*365] days

% test K0 and wmax?
n_dt=length(dt_space); n_C0=length(C0_space); n_P0=length(C0_space); n_wm=length(wm_space); n_K0=length(K0_space);

id=1; % we want to focus on Kangerlussuaq (id=1, KF) and Sermilik (id=2, SF)
for i_C0=1:n_C0
for i_P0=1:n_P0
    dt_run = dt_space(3);
    % name_exp = sprintf('%s_K0%1.e_wm%1.e_C%1.e_P%d_dt%0.2f',fjord_keys{id},K0_space(i_K0),wm_space(i_wm),C0_space(i_C0),P0_space(2),dt_space(1));    
    name_exp = sprintf('%s_C%1.e_P%d_dt%0.2f_nudged',fjord_keys{id},C0_space(i_C0),P0_space(i_P0),dt_run);    
    fprintf('Running model for %s... ',name_exp)
    datasets.opts.dt = dt_run;
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    fjord_run.p.P0=P0_space(i_P0);
    fjord_run.p.C0=C0_space(i_C0);
    fjord_run.p.trelax = 10*86400;
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = name_exp;   
    save([outs_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
    % plot_outputs(fjord_run)
end
end

name_ctrl = sprintf('%s_C%1.e_P%d_dt%0.2f',fjord_keys{id},C0_space(3),P0_space(2),dt_space(3));
% name_ctrl = sprintf('KF_C%1.e_P%d_dt%0.2f',C0_space(3),P0_space(2),dt_space(3));
fjord_ctrl = load([outs_path,'no_nudge/',name_ctrl]).fjord_run;
plot_outputs([],fjord_ctrl)
animate_v4p1(fjord_ctrl,figs_path,name_ctrl,50);

% Checking nudging
name_ctrl_nudge = sprintf('%s_C%1.e_P%d_dt%0.2f_nudged',fjord_keys{id},C0_space(3),P0_space(2),dt_space(3));
fjord_ctrl_nudge = load([outs_path,name_ctrl_nudge]).fjord_run;
plot_outputs([],fjord_ctrl_nudge)


name_set = sprintf('%s_C%1.e_P%d_dt*',fjord_keys{id},C0_space(3),P0_space(2));    % varying dt with default params
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of dt')

name_set = sprintf('%s_C*_P%d_dt%0.2f',fjord_keys{id},P0_space(2),dt_space(4));   % varying C0 with default params
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of C_0')

name_set = sprintf('%s_C%1.e_P*_dt%0.2f',fjord_keys{id},C0_space(3),dt_space(4)); % varying P0 with default params
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of P_0')

% Checking vertical mixing
name_set = sprintf('%s_wm*_C%1.e_P%d_dt%0.2f',fjord_keys{id},C0_space(3),P0_space(2),dt_space(1)); % varying wmax with default params
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of wmax')

name_set = sprintf('%s_K0*wm*_C%1.e_P%d_dt%0.2f',fjord_keys{id},C0_space(3),P0_space(2),0.05/24); % varying wmax and K0 with default params and low time step
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of K0-wmax')