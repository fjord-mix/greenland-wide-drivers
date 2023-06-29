%% Defining the parametre space
M0_space = [2e-8,2e-10,2e-9,2e-7,2e-6]; % default: 2e-8


%% Single-parametre run - relaxation time
for i_param=1:length(M0_space)
    name_exp = sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(i_param));    
    fprintf('Running model for %s... ',name_exp)
    fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure        
    fjord_run.p.M0=M0_space(i_param);
    tic
    fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
    if ~fjord_run.s.status, disp('Done.'); end
    toc    
    fjord_run.o = postprocess_boxmodel(fjord_run);
    fjord_run.m.name = name_exp;   
    save([outs_path,name_exp,'_MPhr.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
    % plot_outputs([],fjord_run);
end

% Control run:
fjord_run = prepare_boxmodel_input(datasets,fjords_processed(fjord_ids(id)),fjord_ids(id)); % arranges into the boxmodel input structure    
% fjord_run.p.M0=M0_space(1);
fjord_run.s = boxmodel_v4(fjord_run.p,fjord_run.f,fjord_run.a,fjord_run.t); % runs and gets the results    
fjord_run.o = postprocess_boxmodel(fjord_run);
% fjord_run.m.name = sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(1));    
fjord_run.m.name = sprintf('%s_M0off_MPhr',fjord_keys{id});
save([outs_path,fjord_run.m.name,'.mat'],'fjord_run','-v7.3') % v7.3 allows files > 2GB
% hf=plot_ts_at_depth(fjord_run,[5,50,300,500,800],'nearest');
% print([figs_path,'iceberg_test/',fjord_run.m.name],'-dpng','-r300')
% animate_v4p1(fjord_run,figs_path,[fjord_run.m.name],50);


name_set = sprintf('%s_M0*_MPhr*',fjord_keys{id});
plot_subset_runs([outs_path,name_set,'.mat'],'Effect of icebergs','Set2')

fjord_run=load([outs_path,sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(end)),'_MPhr.mat']).fjord_run;
plot_fluxes(fjord_run)

for i_param=1:length(M0_space)
    name_exp = sprintf('%s_M0%0.1e',fjord_keys{id},M0_space(i_param));    
    fjord_run = load([outs_path,name_exp,'.mat']).fjord_run;
    hf=plot_ts_at_depth(fjord_run,[5,200,500],'nearest');
    % print([figs_path,'iceberg_test/',name_exp],'-dpng','-r300')
    % close gcf;
    animate_v4p1(fjord_run,figs_path,name_exp,50);
end
