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

for i_reg=1:n_regions
    Xreg = squeeze(X(:,i_reg,:));
    Params_reg = Parameters{i_reg};
    tic
    for k_run=1:n_runs
        try
            % [ohc_out(k_run,i_reg),osc_out(k_run,i_reg)] = wrapper_boxmodel(Xreg(k_run,:),Params_reg);
            [ensemble(k_run,i_reg).time,ensemble(k_run,i_reg).ohc,ensemble(k_run,i_reg).osc] = wrapper_boxmodel(Xreg(k_run,:),Params_reg);
            fprintf('run %d complete\n',k_run)
        catch ME
            ensemble(k_run,i_reg).ohc = NaN(size(ensemble(k_run,i_reg).time));
            ensemble(k_run,i_reg).osc = NaN(size(ensemble(k_run,i_reg).time));
            fprintf('run %d %s\n',k_run,ME.message)
        end
        
        if length(ensemble(k_run,i_reg).ohc) > 365*2
            ohc_start = mean(ensemble(k_run,i_reg).ohc(1:365)); % get the mean for the last year
            osc_start = mean(ensemble(k_run,i_reg).osc(1:365));
        
            ohc_end = mean(ensemble(k_run,i_reg).ohc(end-365:end)); % get the mean for the last year
            osc_end = mean(ensemble(k_run,i_reg).osc(end-365:end));
        
            % get the difference to account for total warming/freshening
            ohc_out(k_run,i_reg) = ohc_end-ohc_start; 
            osc_out(k_run,i_reg) = osc_end-osc_start;
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
    end
    toc
    fprintf('region %d complete\n',i_reg)
end

%% Separate tests of specific runs for debugging model instabilities
% wrapper_boxmodel(X(36,1,:),Parameters{1});
% ohc_out(36,1) = NaN;
% osc_out(36,1) = NaN;

% [ohc_out(k_run,i_reg),osc_out(k_run,i_reg)] = wrapper_boxmodel(X(50,7,:),Parameters{7});
% [ensemble(1,1).time,ensemble(1,1).ohc,ensemble(1,1).osc] = wrapper_boxmodel(X(1,1,:),Parameters{1});
% wrapper_boxmodel(X(6,1,:),Parameters{1});