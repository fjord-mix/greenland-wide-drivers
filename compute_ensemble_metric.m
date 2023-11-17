function [ohc_out,osc_out] = compute_ensemble_metric(ensemble,min_length)
n_runs=size(ensemble,1);
n_regions=size(ensemble,2);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
for i_reg=1:n_regions
    for k_run=1:n_runs
        if length(ensemble(k_run,i_reg).temp) == min_length-1

            [heat_content,salt_content] = get_active_fjord_contents(ensemble(k_run,i_reg));

            ohc_start = mean(heat_content(1:365)); % get the mean for the first year
            osc_start = mean(salt_content(1:365));
        
            ohc_end = mean(heat_content(end-365:end)); % get the mean for the last year
            osc_end = mean(salt_content(end-365:end));

            p = polyfit(ensemble(k_run,i_reg).time,heat_content,1);
            ohc_trend = p(1);
            p = polyfit(ensemble(k_run,i_reg).time,salt_content,1);
            osc_trend = p(1);
        
            % get the difference to account for total warming/freshening
            % ohc_out(k_run,i_reg) = ohc_end-ohc_start; 
            % osc_out(k_run,i_reg) = osc_end-osc_start;
            ohc_out(k_run,i_reg) = ohc_trend *365; % get output "per year"
            osc_out(k_run,i_reg) = osc_trend *365;
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
    end
end

end