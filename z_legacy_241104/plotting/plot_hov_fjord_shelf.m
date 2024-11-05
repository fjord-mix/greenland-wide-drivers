function hf = plot_hov_fjord_shelf(ensemble,time_axis,i_reg,regions_lbl)

n_runs    = size(ensemble,1);
z_profile = -1500:10:0; % bottom up, negative values
depth_lim = -1000; % 1000m because that's virtually the maximum depth at which the sill or GL can be at
sfc_lim   = -25;   % because we limit Zs and Zgl at 25m depth
time_axis = time_axis(2:end); % we discard the initial conditions time step
time_axis_plt = linspace(0,length(time_axis),length(time_axis));


ts_mean = NaN([n_runs,length(z_profile),length(time_axis)]);
tf_mean = NaN(size(ts_mean));
for k_run=1:n_runs
    fjord = ensemble(k_run,i_reg);
    if ~isempty(fjord.temp) && (size(fjord.temp,2) == size(fjord.ts,2))
        
        % bottom up, negative depth values
        zshelf = fjord.zs;
        tshelf = fjord.ts;
        zfjord_series = -flip(cumsum(fjord.H,1),1); 
        tfjord = flip(fjord.temp,1);

        % Interpolating to the 10m resolution grid
        ts_mean(k_run,:,:) = interp1(zshelf,tshelf,z_profile);
        for i_time=1:length(time_axis)
            tf_mean(k_run,:,i_time) = interp1(zfjord_series(:,i_time),tfjord(:,i_time),z_profile,'nearest');%,'extrap');
        end
    end
end

% n_entries = sum(mean(~isnan(tf_mean),3),1);
n_entries = sum(~isnan(tf_mean(:,:,1)),1);
% Averaging considering only the "valid" entries
ts_mean_plt = squeeze(mean(ts_mean,1,'omitnan'));
tf_mean_plt = squeeze(mean(tf_mean,1,'omitnan'));

%% Plotting
hf = figure('Name',sprintf('%s Fjord-shelf comparison',regions_lbl{i_reg}),'Position',[20 20 1000 400]);
subplot(1,3,2); hold on
text(0.05,1.02,'(b) shelf temperature','units','normalized','fontsize',14)
imagesc(time_axis_plt,z_profile,ts_mean_plt);
xlim([0 length(time_axis_plt)]);
ylim([depth_lim sfc_lim])
xlabel('Model time (days)','fontsize',14)
set(gca,'fontsize',14,'YTick',[],'YTickLabel',[])
clims = get(gca,'CLim');

subplot(1,3,1); hold on
text(0.05,1.02,'(a) fjord temperature','units','normalized','fontsize',14)
imagesc(time_axis_plt,z_profile,tf_mean_plt);
xlim([0 length(time_axis_plt)]);
ylim([depth_lim sfc_lim])
ylabel('Depth (m)','fontsize',14); %xlabel('Model time (days)','fontsize',14)
set(gca,'fontsize',14)
clim(clims)
colormap(cmocean('thermal'))
hc = colorbar('fontsize',14);
% ylabel(hc,'Temperature (^oC)');
hc.Position(1)=hc.Position(1)+0.05;
pos_left = get(gca,'Position');

subplot(1,3,3); hold on
text(0.05,1.02,'(c) difference','units','normalized','fontsize',14)
tdiff = tf_mean_plt-ts_mean_plt;
imagesc(time_axis_plt,z_profile,tdiff);
xlim([0 length(time_axis_plt)]);
ylim([depth_lim sfc_lim])
clim([-max(abs(tdiff(:))) max(abs(tdiff(:)))])
colormap(gca,cmocean('balance'))
set(gca,'fontsize',14)

hc2 = colorbar('fontsize',14);
ylabel(hc2,'Temperature difference (^oC)');
hc2.Position(1)=hc2.Position(1)+0.05;

axes('Position',[pos_left(1)-0.12 pos_left(2) .05 pos_left(4)])
plot(n_entries,z_profile,'-k')
xlim([0 max(n_entries)])
set(gca,'fontsize',12,'YTick',[],'YTickLabel',[])

end