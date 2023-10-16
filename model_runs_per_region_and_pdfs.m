for i_reg=1:n_regions
    Xreg = squeeze(X(:,i_reg,:));
    Params_reg = Parameters{i_reg};
    tic
    for k_run=1:n_runs
        try
            % [ensemble(k).time,ensemble(k).ohc,ensemble(k).osc,status,fjord_run(k)] = wrapper_boxmodel(X,Parameters);
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
        else % if the time series is too short
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
    end
    toc
    fprintf('region %d complete\n',i_reg)
end

% This is the trouble boy
% wrapper_boxmodel(X(36,1,:),Parameters{1});
% ohc_out(36,1) = NaN;
% osc_out(36,1) = NaN;

% [ohc_out(k_run,i_reg),osc_out(k_run,i_reg)] = wrapper_boxmodel(X(50,7,:),Parameters{7});
% [ensemble(1,1).time,ensemble(1,1).ohc,ensemble(1,1).osc] = wrapper_boxmodel(X(1,1,:),Parameters{1});
% wrapper_boxmodel(X(6,1,:),Parameters{1});