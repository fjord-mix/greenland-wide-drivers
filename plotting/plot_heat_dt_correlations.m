function [hf, corrs] =  plot_heat_dt_correlations(time_axis,ensemble,regions_lbl)
% which "dT" is the best predictor of heat transport into the fjord?

time_axis_plt = time_axis(2:end);
n_regions = size(ensemble,2);
n_runs    = size(ensemble,1);
heat_transp = NaN([n_runs,n_regions,length(time_axis_plt)]);
dt_reg      = NaN([3,n_runs,n_regions,length(time_axis_plt)]);    
ts_reg      = NaN(size(dt_reg));
tf_reg      = NaN(size(dt_reg));
region_line_color = lines(n_regions);
maxlag = 180;
handles_regions = [];

for i_reg=1:n_regions
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp) && (size(ensemble(k_run,i_reg).temp,2) == size(ensemble(k_run,i_reg).ts,2))
            %% Computes the different dT's
            z_ssfc = 25;
            z_deep = abs(ensemble(k_run,i_reg).p.zgl);
            % z_deep = abs(ensemble(k_run,i_reg).p.silldepth);

            [tf_reg(:,k_run,i_reg,:),~,~] = get_layered_fjord_properties(ensemble(k_run,i_reg),z_ssfc,z_deep);

            zs0 = unique(sort([0,-z_ssfc,-z_deep,ensemble(k_run,i_reg).zs]));

            i_ssfc = find(abs(zs0) == abs(z_ssfc));
            i_deep = find(abs(zs0) == abs(z_deep));
            zs_top = zs0(i_ssfc:end);
            zs_int = zs0(i_deep+1:i_ssfc);
            zs_dep = zs0(1:i_deep);

            if ~isempty(zs_top)
                temp_top = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_top,'pchip','extrap');
                ts_reg(1,k_run,i_reg,:) = squeeze(trapz(zs_top,temp_top))./max(abs(zs_top));
                dt_reg(1,k_run,i_reg,:) = tf_reg(1,k_run,i_reg,:) - ts_reg(1,k_run,i_reg,:);
            end
            if ~isempty(zs_int) && ~isempty(zs_top)
                temp_int = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_int,'pchip','extrap');
                ts_reg(2,k_run,i_reg,:) = squeeze(trapz(zs_int,temp_int))./(max(abs(zs_int)) - max(abs(zs_top)));
                dt_reg(2,k_run,i_reg,:) = tf_reg(2,k_run,i_reg,:) - ts_reg(2,k_run,i_reg,:);
            end
            if ~isempty(zs_dep) && ~isempty(zs_int)
                temp_dep = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_dep,'pchip','extrap');
                ts_reg(3,k_run,i_reg,:) = squeeze(trapz(zs_dep,temp_dep))./(max(abs(zs_dep)) - max(abs(zs_int)));
                dt_reg(3,k_run,i_reg,:) = tf_reg(3,k_run,i_reg,:) - ts_reg(3,k_run,i_reg,:);
            end
        
            %% computes heat transport
            fjord_rho = (ensemble(k_run,i_reg).p.betaS*ensemble(k_run,i_reg).salt ...
                       - ensemble(k_run,i_reg).p.betaT*ensemble(k_run,i_reg).temp);
            heat_transp_layers = (ensemble(k_run,i_reg).temp + 273.15).* fjord_rho .* ensemble(k_run,i_reg).p.cw .* ...
                                       ensemble(k_run,i_reg).qvs;
            heat_transp(k_run,i_reg,:) = sum(heat_transp_layers,1);
        end
    end
end


hf = figure('Name','Correlation between temprature and heat transport','Position',[10 10 1000 800]);
for i_reg=1:n_regions
    ht_region = detrend(squeeze(median(heat_transp(:,i_reg,:),1,'omitnan')).*1e-6);
    for i_layer=1:3
        dt_layer = detrend(squeeze(median(dt_reg(i_layer,:,i_reg,:),2,'omitnan')));
        ts_layer = detrend(squeeze(median(ts_reg(i_layer,:,i_reg,:),2,'omitnan')));
        tf_layer = detrend(squeeze(median(tf_reg(i_layer,:,i_reg,:),2,'omitnan')));

        % compute and plot correlations
        subplot(3,3,i_layer); hold on; box on;
        if i_layer==1, ylabel('dT xcorr'); end
        switch(i_layer)
            case 1, title('top');
            case 2, title('intermediate');
            case 3, title('bottom');
        end
        [xr,lags] = xcorr(ht_region-mean(ht_region),dt_layer-mean(dt_layer),maxlag,'coeff');
        plot(lags,xr,'Color',region_line_color(i_reg,:),'linewidth',2)
        ylim([-1 1])
        xlim([-maxlag maxlag])
        set(gca,'fontsize',14)

        subplot(3,3,i_layer+3); hold on; box on;
        if i_layer==1, ylabel('T_{shelf} xcorr'); end
        [xr,lags] = xcorr(ht_region-mean(ht_region),ts_layer-mean(ts_layer),maxlag,'coeff');
        plot(lags,xr,'Color',region_line_color(i_reg,:),'linewidth',2)
        ylim([-1 1])
        xlim([-maxlag maxlag])
        set(gca,'fontsize',14)

        subplot(3,3,i_layer+6); hold on; box on;
        if i_layer==1, ylabel('T_{fjord} xcorr'); end
        [xr,lags] = xcorr(ht_region-mean(ht_region),tf_layer-mean(tf_layer),maxlag,'coeff');
        plot(lags,xr,'Color',region_line_color(i_reg,:),'linewidth',2)
        ylim([-1 1])
        xlabel('Lag (days)','fontsize',14)
        xlim([-maxlag maxlag])
        set(gca,'fontsize',14)
    end
    hp = plot(lags,xr,'Color',region_line_color(i_reg,:),'linewidth',2);
    handles_regions = [handles_regions hp];
end
hl = legend(handles_regions,regions_lbl,'Location','North','fontsize',10);
hl.NumColumns = 3;


end