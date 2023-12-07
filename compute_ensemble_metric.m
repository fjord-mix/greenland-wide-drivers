function [ohc_out,osc_out] = compute_ensemble_metric(ensemble,min_length)
n_runs=size(ensemble,1);
n_regions=size(ensemble,2);

ohc_out = NaN([n_runs, n_regions]);
osc_out = NaN([n_runs, n_regions]);
for i_reg=1:n_regions
    for k_run=1:n_runs
        if length(ensemble(k_run,i_reg).temp) == min_length-1

            [heat_content,salt_content] = get_active_fjord_contents(ensemble(k_run,i_reg));
            
            zs0 = unique(sort([0,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.silldepth]));
            i_sill = find(zs0 == ensemble(k_run,i_reg).p.silldepth);
            zs0 = zs0(i_sill:end);

            Ss0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs0,'pchip','extrap');
            Ts0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs0,'pchip','extrap');

            hc_shelf = squeeze(trapz(zs0,Ts0)./abs(ensemble(k_run,i_reg).p.silldepth));
            sc_shelf = squeeze(trapz(zs0,Ss0)./abs(ensemble(k_run,i_reg).p.silldepth));
            % taxis_shelf = 1:1:size(sc_reg,1);
            % for i_reg=1:length(regions)
            %     p = polyfit(taxis_shelf,hc_reg(:,i_reg),1);
            %     tr_ohc_shelf(i_reg) = p(1)*12;
            %     % tr_ohc_shelf(i_reg) = mean(hc_reg(end-12:end,i_reg))-mean(hc_reg(1:12,i_reg));
            %     p = polyfit(taxis_shelf,sc_reg(:,i_reg),1);
            %     tr_osc_shelf(i_reg) = p(1)*12;
            %     % tr_osc_shelf(i_reg) = mean(sc_reg(end-12:end,i_reg))-mean(sc_reg(1:12,i_reg));
            % end

            % ohc_start = mean(heat_content(1:365)); % get the mean for the first year
            % osc_start = mean(salt_content(1:365));
            % 
            % ohc_end = mean(heat_content(end-365:end)); % get the mean for the last year
            % osc_end = mean(salt_content(end-365:end));
            % 
            % p = polyfit(ensemble(k_run,i_reg).time,heat_content,1);
            % ohc_trend = p(1);
            % p = polyfit(ensemble(k_run,i_reg).time,salt_content,1);
            % osc_trend = p(1);
        
            % get the difference to account for total warming/freshening
            % ohc_out(k_run,i_reg) = ohc_end-ohc_start; 
            % osc_out(k_run,i_reg) = osc_end-osc_start;
            % ohc_out(k_run,i_reg) = ohc_trend *365; % get output "per year"
            % osc_out(k_run,i_reg) = osc_trend *365;
            ohc_out(k_run,i_reg) = mean(heat_content - hc_shelf);
            osc_out(k_run,i_reg) = mean(salt_content - sc_shelf);
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(k_run,i_reg) = NaN;
            osc_out(k_run,i_reg) = NaN;
        end
    end
end

end