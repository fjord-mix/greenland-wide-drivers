function [hf,tmp_ensemble] = plot_ensemble_dt_ds_layers(ensemble,time_axis,i_tgt_layer,regions_lbl)

% time axis for plotting the results, excluding t0
time_axis_plt = time_axis(2:end); % datetime(2010,01,15)+1:1:datetime(2018,12,15); 
% time_axis_plt = time_axis_plt';
n_regions=size(ensemble,2);
n_runs=size(ensemble,1);
region_handles = [];
region_line_color = lines(n_regions);

% avg_silldepth   = NaN([1,size(ensemble,2)]);
% avg_activedepth = NaN([1,size(ensemble,2)]);
tmp_ensemble     = cell([1,n_regions]);
% dT_ens_mean = NaN([3,length(time_axis_plt)]);
% tf_ens_mean = NaN([3,length(time_axis_plt)]);
% ts_ens_mean = NaN([3,length(time_axis_plt)]);
hf = figure('Name','time series model outputs','Position',[20 20 1200 800]);
hold on; box on
for i_reg=1:n_regions
    ohc_reg=NaN([n_runs,length(time_axis_plt)]);
    osc_reg=NaN([n_runs,length(time_axis_plt)]);
    ohc_shf=NaN([n_runs,length(time_axis_plt)]);
    osc_shf=NaN([n_runs,length(time_axis_plt)]);
    ohc_fjd=NaN([n_runs,length(time_axis_plt)]);
    osc_fjd=NaN([n_runs,length(time_axis_plt)]);
    % hact_fjd=NaN(size(time_axis_plt));
    % zsill_fjd=NaN([1,n_runs]);

    %% Computing the mean T and S over the "active" areas of the shelf (above sill) and fjord (above sill/gl, whichever is deeper)
    for k_run=1:n_runs
        if ~isempty(ensemble(k_run,i_reg).temp) && (size(ensemble(k_run,i_reg).temp,2) == size(ensemble(k_run,i_reg).ts,2))
            z_ssfc = 25;
            z_deep = abs(ensemble(k_run,i_reg).p.zgl);
            % z_deep = abs(ensemble(k_run,i_reg).p.silldepth);

            % [ohc_fjd(k_run,:),osc_fjd(k_run,:),hact_fjd(k_run)] = get_active_fjord_contents(ensemble(k_run,i_reg));
            [tf,sf,~] = get_layered_fjord_properties(ensemble(k_run,i_reg),z_ssfc,z_deep);
            ohc_fjd(k_run,:) = tf(i_tgt_layer,:);
            osc_fjd(k_run,:) = sf(i_tgt_layer,:);
            % avg_silldepth(i_reg)   = avg_silldepth(i_reg)   +abs(ensemble(k_run,i_reg).p.silldepth)./length(ensemble(:,i_reg));
            % h_active = sum(squeeze(sum(hf_out(1:2,:,:),2,'omitnan')));
            % avg_activedepth(i_reg) = avg_activedepth(i_reg) +mean(h_active)./length(ensemble(:,i_reg));
            
            zs0 = unique(sort([0,-z_ssfc,-z_deep,ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).p.silldepth]));
            % z_bottom = max(abs(ensemble(k_run,i_reg).p.silldepth,ensemble(k_run,i_reg).p.zgl));
            
            i_ssfc = find(abs(zs0) == abs(z_ssfc));
            i_deep = find(abs(zs0) == abs(z_deep));
            zs_top = zs0(i_ssfc:end);
            zs_int = zs0(i_deep+1:i_ssfc);
            zs_dep = zs0(1:i_deep);

            switch i_tgt_layer
                case 1
                    if ~isempty(zs_top)
                        salt_top = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_top,'pchip','extrap');
                        temp_top = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_top,'pchip','extrap');
                        ohc_shf(k_run,:) = squeeze(trapz(zs_top,temp_top))./max(abs(zs_top));
                        osc_shf(k_run,:) = squeeze(trapz(zs_top,salt_top))./max(abs(zs_top));
                        ohc_reg(k_run,:) = ohc_fjd(k_run,:) - ohc_shf(k_run,:);
                        osc_reg(k_run,:) = osc_fjd(k_run,:) - osc_shf(k_run,:);
                    end
                case 2
                    if ~isempty(zs_int) && ~isempty(zs_top)
                        salt_int = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_int,'pchip','extrap');
                        temp_int = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_int,'pchip','extrap');
                        ohc_shf(k_run,:) = squeeze(trapz(zs_int,temp_int))./(max(abs(zs_int)) - max(abs(zs_top)));
                        osc_shf(k_run,:) = squeeze(trapz(zs_int,salt_int))./(max(abs(zs_int)) - max(abs(zs_top)));
                        ohc_reg(k_run,:) = ohc_fjd(k_run,:) - ohc_shf(k_run,:);
                        osc_reg(k_run,:) = osc_fjd(k_run,:) - osc_shf(k_run,:);
                    end
                case 3
                    if ~isempty(zs_dep) && ~isempty(zs_int)
                        salt_dep = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ss,zs_dep,'pchip','extrap');
                        temp_dep = interp1(ensemble(k_run,i_reg).zs,ensemble(k_run,i_reg).ts,zs_dep,'pchip','extrap');
                        ohc_shf(k_run,:) = squeeze(trapz(zs_dep,temp_dep))./(max(abs(zs_dep)) - max(abs(zs_int)));
                        osc_shf(k_run,:) = squeeze(trapz(zs_dep,salt_dep))./(max(abs(zs_dep)) - max(abs(zs_int)));
                        ohc_reg(k_run,:) = ohc_fjd(k_run,:) - ohc_shf(k_run,:);
                        osc_reg(k_run,:) = osc_fjd(k_run,:) - osc_shf(k_run,:);
                    end
            end
        end
        % avg_activedepth(i_reg) = mean(hact_fjd,'omitnan');
        % avg_silldepth(i_reg)   = mean(zsill_fjd,'omitnan');
    end

    % Bootstrapping
    dT_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(ohc_reg)));
    tf_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(ohc_fjd)));
    ts_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(ohc_shf)));
    dS_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(osc_reg)));
    sf_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(osc_fjd)));
    ss_ens_mean = mean(bootstrp(1e3,@(x)[median(x,1,'omitnan')],squeeze(osc_shf)));

    fprintf('Layer %d:\n',i_tgt_layer)
    lm = fitlm(datenum(time_axis_plt),dT_ens_mean);
    fprintf("Trend in %s dT: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);
    
    lm = fitlm(datenum(time_axis_plt),tf_ens_mean);
    fprintf("Trend in %s Tf: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);
    
    lm = fitlm(datenum(time_axis_plt),ts_ens_mean);
    fprintf("Trend in %s Ts: %.4f/yr (p=%.2f)\n",regions_lbl{i_reg}(1:2),lm.Coefficients.Estimate(2)*365,lm.ModelFitVsNullModel.Pvalue);
   
    %% Plotting results
    subplot(2,2,1); hold on; box on;    
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_reg);
    std_ln  = prctile(ohc_reg,[25 75],1);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    hp_dt = plot(time_axis_plt,dT_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1.5); 

    subplot(2,2,3); hold on; box on;
    hp_fjd = plot(time_axis_plt,tf_ens_mean,'Color','k','linewidth',1.); 
    plot(time_axis_plt,tf_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1.); 
    hp_shf = plot(time_axis_plt,ts_ens_mean,'Color','k','linewidth',1,'linestyle','--'); 
    plot(time_axis_plt,ts_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 


    subplot(2,2,2); hold on; box on;    
    % std_ln  = bootci(1e3,@(x)[mean(x,1,'omitnan')],ohc_reg);
    std_ln  = prctile(osc_reg,[25 75],1);
    upper_bnd = std_ln(2,:);
    lower_bnd = std_ln(1,:);
    x2 = [time_axis_plt, fliplr(time_axis_plt)];
    inBetween = [lower_bnd, fliplr(upper_bnd)];
    fill(x2, inBetween, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.3);
    plot(time_axis_plt,dS_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1.5); 

    subplot(2,2,4); hold on; box on;
    plot(time_axis_plt,sf_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1.); 
    plot(time_axis_plt,ss_ens_mean,'Color',region_line_color(i_reg,:),'linewidth',1,'linestyle','--'); 
    
    region_handles = [region_handles hp_dt];

    % getting a time table with all T and S time series
    tmp_ensemble{i_reg} = timetable(time_axis_plt',dT_ens_mean',tf_ens_mean',ts_ens_mean',...
                                                   dS_ens_mean',sf_ens_mean',ss_ens_mean',...
                                    'VariableNames',{'dT','Tf','Ts','dS','Sf','Ss'});
end
    
%% Adding labels to panels and such
subplot(2,2,1)
text(0.03,1.07,'(a)','fontsize',14,'units','normalized')
ylabel('Temperature difference (^oC)'); %xlabel('Time');
hl = legend(region_handles,regions_lbl,'fontsize',10,'Location','southeast');
hl.NumColumns=2;
set(gca,'fontsize',14)

subplot(2,2,3)
text(0.03,1.07,'(b)','fontsize',14,'units','normalized')
ylabel('Temperature (^oC)'); xlabel('Time');
hl2 = legend([hp_fjd,hp_shf],{'Fjord','Shelf'},'fontsize',10,'Location','southeast');
set(gca,'fontsize',14)

subplot(2,2,2)
text(0.03,1.07,'(c)','fontsize',14,'units','normalized')
ylabel('Salinity difference'); %xlabel('Time');
set(gca,'fontsize',14)

subplot(2,2,4)
text(0.03,1.07,'(d)','fontsize',14,'units','normalized')
ylabel('Salinity'); xlabel('Time');
set(gca,'fontsize',14)
end