function [ohc_out,osc_out,avg_silldepth,avg_activedepth] = compute_ensemble_metric(ensemble,min_length)
n_runs=size(ensemble,1);
n_regions=size(ensemble,2);
avg_silldepth   = zeros(size(ensemble,2));
avg_activedepth = zeros(size(ensemble,2));

ohc_out = NaN([3,n_runs, n_regions]);
osc_out = NaN([3,n_runs, n_regions]);
for i_reg=1:n_regions
    for k_run=1:n_runs
        if (length(ensemble(k_run,i_reg).temp) == min_length-1) && all(~isnan(ensemble(k_run,i_reg).temp(:)))
            % [heat_content,salt_content,activedepth] = get_active_fjord_contents(ensemble(k_run,i_reg));
            [tf_out,sf_out,hf_out] = get_layered_fjord_properties(ensemble(k_run,i_reg));

            avg_silldepth(i_reg)   = avg_silldepth(i_reg)   +abs(ensemble(k_run,i_reg).p.silldepth)./length(ensemble(:,i_reg));
            h_active = sum(squeeze(sum(hf_out(1:2,:,:),2,'omitnan')));
            avg_activedepth(i_reg) = avg_activedepth(i_reg) +mean(h_active)./length(ensemble(:,i_reg));
            
            z_ssfc = 25;
            z_deep = abs(ensemble(k_run,i_reg).p.zgl);
            zs0 = unique(sort([0,-z_ssfc,-z_deep,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.silldepth]));
            % z_bottom = max(abs(ensemble(k_run,i_reg).p.silldepth,ensemble(k_run,i_reg).p.zgl));
            
            i_ssfc = find(abs(zs0) == abs(z_ssfc));
            i_deep = find(abs(zs0) == abs(z_deep));
            zs_top = zs0(i_ssfc:end);
            zs_int = zs0(i_deep+1:i_ssfc);
            zs_dep = zs0(1:i_deep);

            % TODO: best definition of these depth intervals or limiting Zg and Zs
            if ~isempty(zs_top)
                salt_top = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_top,'pchip','extrap');
                temp_top = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_top,'pchip','extrap');
                ts_top = squeeze(trapz(zs_top,temp_top))./max(abs(zs_top));
                ss_top = squeeze(trapz(zs_top,salt_top))./max(abs(zs_top));
                ohc_out(1,k_run,i_reg) = mean(tf_out(1,:) - ts_top);
                osc_out(1,k_run,i_reg) = mean(sf_out(1,:) - ss_top);
            end
            if ~isempty(zs_int) && ~isempty(zs_top)
                salt_int = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_int,'pchip','extrap');
                temp_int = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_int,'pchip','extrap');
                ts_int = squeeze(trapz(zs_int,temp_int))./(max(abs(zs_int)) - max(abs(zs_top)));
                ss_int = squeeze(trapz(zs_int,salt_int))./(max(abs(zs_int)) - max(abs(zs_top)));
                ohc_out(2,k_run,i_reg) = mean(tf_out(2,:) - ts_int);
                osc_out(2,k_run,i_reg) = mean(sf_out(2,:) - ss_int);
            end
            if ~isempty(zs_dep) && ~isempty(zs_int)
                salt_dep = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_dep,'pchip','extrap');
                temp_dep = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_dep,'pchip','extrap');
                ts_dep = squeeze(trapz(zs_dep,temp_dep))./(max(abs(zs_dep)) - max(abs(zs_int)));
                ss_dep = squeeze(trapz(zs_dep,salt_dep))./(max(abs(zs_dep)) - max(abs(zs_int)));
                ohc_out(3,k_run,i_reg) = mean(tf_out(3,:) - ts_dep);
                osc_out(3,k_run,i_reg) = mean(sf_out(3,:) - ss_dep);
            end    
            
        else 
            % if the time series is too short, i.e., the model crashed
            ohc_out(:,k_run,i_reg) = NaN;
            osc_out(:,k_run,i_reg) = NaN;
        end
        if sum(isinf(ohc_out(:,k_run,i_reg))) > 0 || sum(isinf(osc_out(:,k_run,i_reg))) > 0
            ohc_out(:,k_run,i_reg) = NaN;
            osc_out(:,k_run,i_reg) = NaN;
        end
    end
end

end