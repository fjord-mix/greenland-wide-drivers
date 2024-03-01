function [hf,tt_ensemble] = plot_ensemble_dt_ds(ensemble,time_axis,regions_lbl)

% time axis for plotting the results, excluding t0
time_axis_plt = time_axis(2:end); % datetime(2010,01,15)+1:1:datetime(2018,12,15); 
% time_axis_plt = time_axis_plt';
n_regions=size(ensemble,2);
n_runs=size(ensemble,1);
region_handles = [];
region_line_color = lines(n_regions);

avg_silldepth   = NaN([1,size(ensemble,2)]);
avg_activedepth = NaN([1,size(ensemble,2)]);
tt_ensemble     = cell([1,n_regions]);
hf = figure('Name','time series model outputs','Position',[20 20 1200 800]);
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis_plt)]);
    osc_reg=NaN([n_runs,length(time_axis_plt)]);
    ohc_shf=NaN([n_runs,length(time_axis_plt)]);
    osc_shf=NaN([n_runs,length(time_axis_plt)]);
    ohc_fjd=NaN([n_runs,length(time_axis_plt)]);
    osc_fjd=NaN([n_runs,length(time_axis_plt)]);
    qsg_reg=NaN([n_runs,length(time_axis_plt)]);
    hact_fjd=NaN(size(time_axis_plt));
    zsill_fjd=NaN([1,n_runs]);

    %% Computing the mean T and S over the "active" areas of the shelf (above sill) and fjord (above sill/gl, whichever is deeper)
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp) && (size(ensemble(k_run,i_reg).temp,2) == size(ensemble(k_run,i_reg).ts,2))
            [ohc_fjd(k_run,:),osc_fjd(k_run,:),hact_fjd(k_run)] = get_active_fjord_contents(ensemble(k_run,i_reg));

            zs0 = unique(sort([0,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.zgl,ensemble(k_run,i_reg).p.silldepth]));
            % z_bottom = max(abs(ensemble(k_run,i_reg).p.silldepth,ensemble(k_run,i_reg).p.zgl));
            z_bottom = abs(ensemble(k_run,i_reg).p.zgl);
            i_bottom = find(abs(zs0) == abs(z_bottom));
            zs0 = zs0(i_bottom:end);

            Ss0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs0,'pchip','extrap');
            Ts0 = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs0,'pchip','extrap');

            % ohc_shf(k_run,:) = squeeze(trapz(zs0,Ts0)./abs(ensemble(k_run,i_reg).p.silldepth));
            % osc_shf(k_run,:) = squeeze(trapz(zs0,Ss0)./abs(ensemble(k_run,i_reg).p.silldepth));
            ohc_shf(k_run,:) = squeeze(trapz(zs0,Ts0)./abs(hact_fjd(k_run)));
            osc_shf(k_run,:) = squeeze(trapz(zs0,Ss0)./abs(hact_fjd(k_run)));

            ohc_reg(k_run,:) = ohc_fjd(k_run,:) - ohc_shf(k_run,:);
            osc_reg(k_run,:) = osc_fjd(k_run,:) - osc_shf(k_run,:);

            zsill_fjd(k_run) = abs(ensemble(k_run,i_reg).p.silldepth);
        else
            ohc_reg(k_run,:) = NaN;
            osc_reg(k_run,:) = NaN;
            qsg_reg(k_run,:) = NaN;
        end
        avg_activedepth(i_reg) = mean(hact_fjd,'omitnan');
        avg_silldepth(i_reg)   = mean(zsill_fjd,'omitnan');
    end

    %% Testing for significance of trends in d(T,S), (T,S)f, and (T,S)s
    mean_ln_temp = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],ohc_reg));
    % mean_ln_temp = median(ohc_reg,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_temp);
    fprintf("Trend in %s dT: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    mean_ln_tf = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],ohc_fjd));
    % mean_ln_tf = median(ohc_fjd,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_tf);
    fprintf("Trend in %s Tf: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    mean_ln_ts = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],ohc_shf));
    % mean_ln_ts = median(ohc_shf,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_ts);
    fprintf("Trend in %s Ts: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    mean_ln_salt = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],osc_reg));
    % mean_ln_salt = median(osc_reg,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_salt);
    fprintf("Trend in %s dS: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    mean_ln_sf = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],osc_fjd));
    % mean_ln_sf = median(osc_fjd,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_sf);
    fprintf("Trend in %s Sf: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    mean_ln_ss = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],osc_shf));
    % mean_ln_ss = median(osc_shf,1,'omitnan');
    lm = fitlm(datenum(time_axis_plt),mean_ln_ss);
    fprintf("Trend in %s Ss: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);

    %% Plotting results
    subplot(2,2,1); hold on; box on;    
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_reg);
    std_ln  = prctile(ohc_reg,[25 75],1);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    plot(time_axis_plt,mean_ln_temp,'Color',region_line_color(i_reg,:),'linewidth',1.5); 
    
    subplot(2,2,2); hold on; box on;
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],osc_reg);
    std_ln  = prctile(osc_reg,[25 75],1);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp = plot(time_axis_plt,mean_ln_salt,'Color',region_line_color(i_reg,:),'linewidth',1.5); 
    region_handles=[region_handles hp];

    subplot(2,2,3); hold on; box on;
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_fjd);
    % upper_bnd = std_ln(2,:);
    % lower_bnd = std_ln(1,:);
    % x2 = [time_axis_plt, fliplr(time_axis_plt)];
    % inBetween = [lower_bnd, fliplr(upper_bnd)];
    % fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp_fjd = plot(time_axis_plt,mean_ln_tf,'Color','k','linewidth',1.); 
    plot(time_axis_plt,mean_ln_tf,'Color',region_line_color(i_reg,:),'linewidth',1.); 
    hp_shf = plot(time_axis_plt,mean_ln_ts,'Color','k','linewidth',1,'linestyle','--'); 
    plot(time_axis_plt,mean_ln_ts,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 
    
    subplot(2,2,4); hold on; box on;
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],osc_fjd);
    % upper_bnd = std_ln(2,:);
    % lower_bnd = std_ln(1,:);
    % x2 = [time_axis_plt, fliplr(time_axis_plt)];
    % inBetween = [lower_bnd, fliplr(upper_bnd)];
    % fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp = plot(time_axis_plt,mean_ln_sf,'Color',region_line_color(i_reg,:),'linewidth',1); 
    plot(time_axis_plt,mean_ln_ss,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 


    % getting a time table with all T and S time series
    tt_ensemble{i_reg} = timetable(time_axis_plt',mean_ln_temp',mean_ln_salt',mean_ln_tf',mean_ln_ts',mean_ln_sf',mean_ln_ss',...
        'VariableNames',{'dT','dS','Tf','Ts','Sf','Ss'});

end
%% Adding labels to panels and such
subplot(2,2,1)
text(0.03,1.07,'(a)','fontsize',14,'units','normalized')
ylabel('Temperature difference (^oC)'); % xlabel('Time');
set(gca,'fontsize',14)
subplot(2,2,2)
text(0.03,1.07,'(b)','fontsize',14,'units','normalized')
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','southeast');
hl.NumColumns=3;
ylabel('Salinity difference'); % xlabel('Time'); 
set(gca,'fontsize',14)

subplot(2,2,3)
text(0.03,1.07,'(c)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Temperature (^oC)');
hl2 = legend([hp_fjd hp_shf],{'fjord','shelf'},'fontsize',10,'Location','northeast');
set(gca,'fontsize',14)
subplot(2,2,4)
text(0.03,1.07,'(d)','fontsize',14,'units','normalized')
xlabel('Time'); ylabel('Salinity');
set(gca,'fontsize',14)

end