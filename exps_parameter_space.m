%% Exploring the parametre space
id=1;
[p_def,~] = get_model_default_parameters();
% if id==2 || id==4 % we remove the sill for Sermilik and Ilulissat
%     p_def.sill=0;
% end
C0_space  = [1e4, 1e3, 5e3, 5e4];      % shelf-exchange coefficient (default: 1e4)
P0_space  = [25, 10 50, 75];           % plume-entrainment coefficient (default: 25)
K0_space  = [5e-2, 5e-10, 5e-5, 5e-3]; % mixing coefficient (default 0.05)
tau_space = [5, 10, 30, 90, 365];      % relaxation time
wmx_space = [4e-5, 1e-4, 1e-5, 1e-6];  % maximum entrainment velocity

%% Single-parametre run - relaxation time
for i_par=1:length(tau_space)
    name_exp = sprintf('%s_tau%dd',fjord_keys{id},tau_space(i_par));    
    fprintf('Running model for %s... ',name_exp)    
    fjord_run = prepare_boxmodel_input(datasets_mod,fjords_processed(fjord_ids(id)),fjord_ids(id),p_def); % arranges into the boxmodel input structure    
    fjord_run.p.trelax = tau_space(i_par);
    tic
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p,fjord_run.t,fjord_run.f,fjord_run.a); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.o = postprocess_boxmodel(fjord_run);
    fjord_run.m.name = name_exp;   
    % save([exp_out_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
end

% Plotting outputs
name_set = sprintf('%s_tau*',fjord_keys{id});
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of relaxation time','Paired')
% exportgraphics(gcf,[figs_path,'/test_rlx/fluxes_trelax.png'],'ContentType','vector','BackgroundColor','none'); close gcf

for i_par=1:length(tau_space)
    name_exp = sprintf('%s_tau%dd',fjord_keys{id},tau_space(i_par));    
    fjord_run = load([exp_out_path,name_exp,'.mat']).fjord_run;
    hf=plot_ts_at_depth(fjord_run,[5,200,500],'nearest');
    % exportgraphics(gcf,[figs_path,'/',name_exp,'.pdf'],'ContentType','vector','BackgroundColor','none'); close gcf
end

%% Single-parametre run - maximum entrainment velocity
exp_out_path=[outs_path,'/sens_mix/'];
mkdir(exp_out_path)
id=3;
for i_par=1:length(wmx_space)
    name_exp = sprintf('%s_wmax%e',fjord_keys{id},wmx_space(i_par));    
    fprintf('Running model for %s... ',name_exp)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p_def); % arranges into the boxmodel input structure    
    fjord_run.p.wmax = wmx_space(i_par);
    tic
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p,fjord_run.t,fjord_run.f,fjord_run.a); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.o = postprocess_boxmodel(fjord_run);
    fjord_run.m.name = name_exp;   
    save([exp_out_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
end

% Plotting outputs
name_set = sprintf('%s_wmax*',fjord_keys{id});
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of maximum entrainment velocity','Paired')
% exportgraphics(gcf,[figs_path,'/test_rlx/fluxes_trelax.png'],'ContentType','vector','BackgroundColor','none'); close gcf

%% Double parametre runs - tau vs. C0
exp_out_path=[outs_path,'sens_multi_params/'];
mkdir(exp_out_path)
for i_C0=1:length(C0_space)
for i_tau=1:length(tau_space)
    name_exp = sprintf('%s_C%1.e_tau%dd',fjord_keys{id},C0_space(i_C0),tau_space(i_tau));    
    fprintf('Running model for %s... ',name_exp)
    datasets_mod = datasets;
    datasets_mod.opts.dt = 0.25;
    fjord_run = prepare_boxmodel_input(datasets_mod,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure
    fjord_run.p.trelax = tau_space(i_tau);
    fjord_run.p.C0=C0_space(i_C0);
    tic
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p,fjord_run.t,fjord_run.f,fjord_run.a); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = name_exp;
    save([exp_out_path,name_exp,'_mod.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
end
end

name_set = sprintf('%s_C*tau5d',fjord_keys{id});
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of C0 on short trelax','Set3')

%% Multi-parametre runs
exp_out_path=[outs_path,'sens_multi_params/'];
mkdir(exp_out_path)
id=1; % we want to focus on Kangerlussuaq (id=1, KF) and Sermilik (id=2, SF)
for i_C0=1:length(C0_space)
for i_P0=1:length(P0_space)
for i_K0=1:length(K0_space)
    name_exp = sprintf('%s_K%1.e_C%1.e_P%d',fjord_keys{id},K0_space(i_K0),C0_space(i_C0),P0_space(i_P0));    
    fprintf('Running model for %s... ',name_exp)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id),p_def); % arranges into the boxmodel input structure
    fjord_run.p.P0=P0_space(i_P0);
    fjord_run.p.C0=C0_space(i_C0);
    fjord_run.p.K0=K0_space(i_K0);
    tic
    [fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p,fjord_run.t,fjord_run.f,fjord_run.a); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.m.name = name_exp;   
    save([exp_out_path,name_exp,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
end
end
end

% Plotting outputs
name_ctrl = sprintf('%s_K%1.e_C%1.e_P%d',fjord_keys{id},K0_space(1),C0_space(1),P0_space(1));    
fjord_ctrl = load([exp_out_path,name_ctrl,'.mat']).fjord_run;
plot_outputs(fjord_ctrl);
% exportgraphics(gcf,[figs_path,'/',name_ctrl,'.pdf'],'ContentType','vector','BackgroundColor','none');
% animate([exp_out_path,name_ctrl,'.mat'],figs_path,name_ctrl,50);

name_set = sprintf('%s_K*_C%1.e_P%d',fjord_keys{id},C0_space(1),P0_space(1)); % varying K0 with default params
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of K_0');
% exportgraphics(gcf,[figs_path,'/',name_set,'.pdf'],'ContentType','vector','BackgroundColor','none');

name_set = sprintf('%s_K%1.e_C*_P%d',fjord_keys{id},K0_space(1),P0_space(1)); % varying C0 with default params
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of C_0');
% exportgraphics(gcf,[figs_path,'/',name_set,'.pdf'],'ContentType','vector','BackgroundColor','none');

name_set = sprintf('%s_K%1.e_C%1.e_P*',fjord_keys{id},K0_space(1),C0_space(1)); % varying P0 with default params
plot_subset_runs([exp_out_path,name_set,'.mat'],'Effect of P_0');
% exportgraphics(gcf,[figs_path,'/',name_set,'.pdf'],'ContentType','vector','BackgroundColor','none');