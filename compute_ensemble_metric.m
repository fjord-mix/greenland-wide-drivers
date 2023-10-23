function [ohc_out,osc_out] = compute_ensemble_metric(ensemble,min_length)
n_runs=size(ensemble,1);
n_regions=size(ensemble,2);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
for i_reg=1:n_regions
    for k_run=1:n_runs
        if length(ensemble(k_run,i_reg).ohc) == min_length
            
            % compute quantities per unit volume
            fjord_run = ensemble(k_run,i_reg);
            heat_content = sum(fjord_run.ohc,1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);
            salt_content = sum(fjord_run.osc,1)./(fjord_run.p.L.*fjord_run.p.W.*fjord_run.p.H);

            ohc_start = mean(heat_content(1:365)); % get the mean for the first year
            osc_start = mean(salt_content(1:365));
        
            ohc_end = mean(heat_content(end-365:end)); % get the mean for the last year
            osc_end = mean(salt_content(end-365:end));
        
            % get the difference to account for total warming/freshening
            ohc_out(k_run,i_reg) = ohc_end-ohc_start; 
            osc_out(k_run,i_reg) = osc_end-osc_start;
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
    end
end

end