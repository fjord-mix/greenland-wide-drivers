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
            [ensemble(k_run,i_reg).time,ensemble(k_run,i_reg).ohc,ensemble(k_run,i_reg).osc,...
                                        ensemble(k_run,i_reg).ohc_as,ensemble(k_run,i_reg).osc_as] = ...
             wrapper_boxmodel(Xreg(k_run,:),Params_reg);
            fprintf('run %d complete\n',k_run)
        catch ME
            ensemble(k_run,i_reg).ohc = NaN(size(ensemble(k_run,i_reg).time));
            ensemble(k_run,i_reg).osc = NaN(size(ensemble(k_run,i_reg).time));
            ensemble(k_run,i_reg).ohc_as = NaN(size(ensemble(k_run,i_reg).time));
            ensemble(k_run,i_reg).osc_as = NaN(size(ensemble(k_run,i_reg).time));
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

            % Computing for only layers above the sill
            ohc_start_as = mean(ensemble(k_run,i_reg).ohc_as(1:365)); 
            osc_start_as = mean(ensemble(k_run,i_reg).osc_as(1:365));
            ohc_end_as   = mean(ensemble(k_run,i_reg).ohc_as(end-365:end)); 
            osc_end_as   = mean(ensemble(k_run,i_reg).osc_as(end-365:end));
            ohc_out_as(k_run,i_reg) = ohc_end_as-ohc_start_as; 
            osc_out_as(k_run,i_reg) = osc_end_as-osc_start_as;
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
            ohc_out_as(k_run,i_reg) = NaN;
            osc_out_as(k_run,i_reg) = NaN;
        end
    end
    toc
    fprintf('region %d complete\n',i_reg)
end

%% Calculate the distributions based on the numerical outputs alone
% this is just for comparison with the surrogate model outputs
% load([outs_path,'ohc_osc_runs_probs_n',num2str(n_runs)]) % if we have the results saved already

for i_reg=1:n_regions
    ohc_pd{i_reg} = makedist('Normal','mu',mean(ohc_out(:,i_reg),'omitnan'),'sigma',std(ohc_out(:,i_reg),'omitnan'));
    osc_pd{i_reg} = makedist('Normal','mu',mean(osc_out(:,i_reg),'omitnan'),'sigma',std(osc_out(:,i_reg),'omitnan'));
    ohc_ks{i_reg} = fitdist(ohc_out(:,i_reg),'kernel');
    osc_ks{i_reg} = fitdist(osc_out(:,i_reg),'kernel');
    ohc_ks_as{i_reg} = fitdist(ohc_out_as(:,i_reg),'kernel');
    osc_ks_as{i_reg} = fitdist(osc_out_as(:,i_reg),'kernel');
end

%% Separate tests of specific runs for debugging model instabilities
% wrapper_boxmodel(X(36,1,:),Parameters{1});
% ohc_out(36,1) = NaN;
% osc_out(36,1) = NaN;

% [ohc_out(k_run,i_reg),osc_out(k_run,i_reg)] = wrapper_boxmodel(X(50,7,:),Parameters{7});
% [ensemble(1,1).time,ensemble(1,1).ohc,ensemble(1,1).osc] = wrapper_boxmodel(X(1,1,:),Parameters{1});
% wrapper_boxmodel(X(7,5,:),Parameters{5});