% Uncomment the code snipper below if running in parallel,
% and remember to substitute the inner for loop for a parfor!
% Obs. (1) it might be necessary to "flatten" the ensemble structure array
% Obs. (2) using a parfor inside the for loop will likely increase overhead time, 
%          but will significantly reduce the amount of memory used
%          by enabling us to slice X and Prameters before the parallel runs

% num_workers=4;
% checkPool = gcp('nocreate'); % If no pool, do not create one
% if isempty(checkPool) % if there is no pool
%     parpool(num_workers);
% end

%% Training dataset
for i_reg=1:n_regions
    Xreg = squeeze(X(:,i_reg,:));
    Params_reg = Parameters{i_reg};
    tic
    for k_run=1:n_runs
        % if isnan(ohc_out(k_run,i_reg)) % saves time after a bug/crash fix, but requires reinitialising the variable        
        ensemble(k_run,i_reg) = wrapper_boxmodel(Xreg(k_run,:),Params_reg);
        fprintf('run %d complete\n',k_run)
        % end
    end
    fprintf('region %d complete\n',i_reg)
    toc    
end
fprintf('Model training runs complete. Starting computation of heat/salt contents...\n')

[ohc_out,osc_out] = compute_ensemble_metric(ensemble,length(time_axis));
fprintf('Computation of heat/salt contents complete. Starting computation of heat/salt contents...\n')

%% Calculate the distributions based on the numerical outputs alone
% this is just for comparison with the surrogate model outputs
% load([outs_path,'ohc_osc_runs_probs_n',num2str(n_runs)]) % if we have the results saved already

% ohc_pd  = cell([1,n_regions]);
% osc_pd  = cell([1,n_regions]);
% ohc_ks  = cell([1, n_regions]);
% osc_ks  = cell([1, n_regions]);
% for i_reg=1:n_regions
%     % ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
%     % osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
%     ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
%     osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
% end

