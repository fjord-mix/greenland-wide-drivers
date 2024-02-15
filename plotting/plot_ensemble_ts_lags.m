function hf = plot_ensemble_ts_lags(ensemble,maxlag)
if nargin < 2, maxlag = 360; end
region_line_color = lines(7);
xcor_ohc=NaN(size(ensemble));
lags_ohc=NaN(size(ensemble));
xcor_osc=NaN(size(ensemble));
lags_osc=NaN(size(ensemble));
n_regions=size(ensemble,2);
n_runs=size(ensemble,1);
for i_reg=1:n_regions
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp)
            [ohc_fjord,osc_fjord,active_depth] = get_active_fjord_contents(ensemble(k_run,i_reg));
            
            zs0 = unique(sort([0,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.zgl,ensemble(k_run,i_reg).p.silldepth]));
            % z_bottom = max(abs(ensemble(k_run,i_reg).p.silldepth,ensemble(k_run,i_reg).p.zgl));
            z_bottom = abs(ensemble(k_run,i_reg).p.zgl);
            i_bottom = find(abs(zs0) == abs(z_bottom));
            zs0 = zs0(i_bottom:end);

            Ss0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs0,'pchip','extrap');
            Ts0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs0,'pchip','extrap');

            % ohc_shelf = squeeze(trapz(zs0,Ts0)./abs(ensemble(k_run,i_reg).p.silldepth));
            % osc_shelf = squeeze(trapz(zs0,Ss0)./abs(ensemble(k_run,i_reg).p.silldepth));
            ohc_shelf = squeeze(trapz(zs0,Ts0)./abs(active_depth));
            osc_shelf = squeeze(trapz(zs0,Ss0)./abs(active_depth));

            %TODO: detrend and normalise by STD before computing
            [r_ohc,l_ohc] = xcorr(ohc_fjord-mean(ohc_fjord),ohc_shelf-mean(ohc_shelf),maxlag,'coeff');
            [xcor_ohc(k_run,i_reg),i_maxr] = max(r_ohc);
            lags_ohc(k_run,i_reg)          = l_ohc(i_maxr);
            
            [r_osc,l_osc] = xcorr(osc_fjord-mean(osc_fjord),osc_shelf-mean(osc_shelf),maxlag,'coeff');
            [xcor_osc(k_run,i_reg),i_maxr] = max(r_osc);
            lags_osc(k_run,i_reg)          = l_osc(i_maxr);
        end
    end
end

lags_xi = -maxlag:0.1:maxlag;
corr_xi = linspace(min([xcor_ohc(:); xcor_osc(:)]),max([xcor_ohc(:); xcor_osc(:)]),1e3);
hf = figure('Name','Lags between fjord and shelf','Position',[20 20 800 800]);
for i_reg=1:n_regions
    corr_ohc_kern   = fitdist(xcor_ohc(:,i_reg),'kernel');
    corr_osc_kern   = fitdist(xcor_osc(:,i_reg),'kernel');
    lags_ohc_kern   = fitdist(lags_ohc(:,i_reg),'kernel');
    lags_osc_kern   = fitdist(lags_osc(:,i_reg),'kernel');
    
    subplot(2,2,1); hold on; box on; grid on;
    plot(lags_xi,pdf(lags_ohc_kern,lags_xi),'color',region_line_color(i_reg,:),'linewidth',2)
    subplot(2,2,2); hold on; box on; grid on;
    plot(lags_xi,pdf(lags_osc_kern,lags_xi),'color',region_line_color(i_reg,:),'linewidth',2)
    subplot(2,2,3); hold on; box on; grid on;
    plot(corr_xi,pdf(corr_ohc_kern,corr_xi),'color',region_line_color(i_reg,:),'linewidth',2)
    subplot(2,2,4); hold on; box on; grid on;
    plot(corr_xi,pdf(corr_osc_kern,corr_xi),'color',region_line_color(i_reg,:),'linewidth',2)
end
subplot(2,2,1); 
ylabel('Probability function')
xlabel('Lag in temperature (days)')
set(gca,'fontsize',14)
% xlim([0 maxlag/2])
subplot(2,2,2); 
xlabel('Lag in salinity (days)')
set(gca,'fontsize',14)
% xlim([0 maxlag/2])
subplot(2,2,3); 
ylabel('Probability function')
xlabel('Cross-correlation of temperature')
set(gca,'fontsize',14)
% xlim([0.5 1])
subplot(2,2,4); 
xlabel('Cross-correlation of salinity')
set(gca,'fontsize',14)
% xlim([0.5 1])
end