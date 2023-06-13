function fig_layer_handles = plot_TS_boxmodel(datasets,fjords_map,fjord_title,fjord,model_runtime,all_ctd_data)
letters = 'abcdefgh';

% bin shelf forcings according to the fjord boxes
n_layers=size(fjord.s.H,1);
n_timesteps=length(fjord.s.t);
Ts = zeros([n_layers,n_timesteps]);
Ss = zeros([n_layers,n_timesteps]);

% increasing resolution to avoid trapz integral issues
zs = fjord.f.zs;
nz_orig=1:length(zs);
nz_hr=linspace(1,length(zs),100*length(zs));
zs_hr=interp1(nz_orig,zs,nz_hr);
Ts_hr=interp1(zs,fjord.f.Ts,zs_hr);
Ss_hr=interp1(zs,fjord.f.Ss,zs_hr);
for k=1:n_timesteps    
    [Ts(:,k),Ss(:,k)] = bin_ocean_profiles(Ts_hr(:,k),Ss_hr(:,k),-zs_hr,fjord.s.H(:,k)',fjord.p);
end
% smooth out the profile for plotting purposes. The differences in salinity are O(10^-2)
Ts=smoothdata(Ts,2);
Ss=smoothdata(Ss,2);

% Ts = Ts';
% Ss = Ss';

% get fjord limits to select obs data accordingly
[Xf,Yf]    = meshgrid(fjords_map.x,fjords_map.y);
fjord_lims = [min(Xf(fjord.m.inds)), min(Yf(fjord.m.inds)); max(Xf(fjord.m.inds)), max(Yf(fjord.m.inds))];

% find which casts are in the fjord or not - refine the algorithm later!
n_in_fjord=0;
for i=1:length(all_ctd_data)
    try
    if (all_ctd_data(i).x > fjord_lims(1,1) && all_ctd_data(i).x < fjord_lims(2,1)) && ...
       (all_ctd_data(i).y > fjord_lims(1,2) && all_ctd_data(i).y < fjord_lims(2,2))
        all_ctd_data(i).inFjord = 1;
        n_in_fjord = n_in_fjord+1;
    else
        all_ctd_data(i).inFjord = 0;
    end
    catch
        disp('this is the trouble boy') % some data from Zackenberg/Nuuk might have unsolved problems here
    end
end

if n_in_fjord == 0
    disp('No casts found inside fjord!')
    casts_time_series = [];
else
    % subet casts that lie inside the fjord
    casts_in_fjord(n_in_fjord) = struct();
    i_in_fjord=1;
    for i=1:length(all_ctd_data)
        if all_ctd_data(i).inFjord
            casts_in_fjord(i_in_fjord).x     = all_ctd_data(i).x;
            casts_in_fjord(i_in_fjord).y     = all_ctd_data(i).y;
            casts_in_fjord(i_in_fjord).temp  = all_ctd_data(i).temp;
            casts_in_fjord(i_in_fjord).sal   = all_ctd_data(i).sal ;
            casts_in_fjord(i_in_fjord).depth = all_ctd_data(i).depth;
            casts_in_fjord(i_in_fjord).time  = all_ctd_data(i).time;
            i_in_fjord=i_in_fjord+1;
        end
    end
    casts_time_series = casts_in_fjord;
end
% % organise them by time
% if length(casts_in_fjord) > 1
%     T = struct2table(casts_in_fjord);
%     casts_time_series = sortrows(T,'time');
%     casts_time_series.time = datetime(casts_time_series.time,'ConvertFrom','datenum');
%     casts_time_series = table2struct(casts_time_series)';
% else
%     casts_time_series = casts_in_fjord;
%     casts_time_series.time = datetime(casts_time_series.time,'ConvertFrom','datenum');
% end


% check plot
% fig_handles = plot_obs(datasets,[],casts_in_fjord);
% print('-dpng','-r300',['/Users/mmeb1/FjordMIX/presentations/2305_StAndrews/figs/fjord_casts_',fjord_title])
%% Determining plot parameters
n_layers = size(fjord.s.H,1);
tlims = NaN([n_layers,2]);
slims = NaN([n_layers,2]);
for i=1:n_layers
    min_t = min(min(fjord.s.T(i,:)),min(Ts(i,:)));    
    tlims(i,:) = [min_t-0.5 max(max(fjord.s.T(i,:)),max(Ts(i,:)))+0.5];
    slims(i,:) = [min(min(fjord.s.S(i,:)),min(Ss(i,:)))-0.2 max(max(fjord.s.S(i,:)),max(Ss(i,:)))+0.2];
end

runtime_axis = fjord.m.time_axis;
t0 = convertTo(runtime_axis(1),'datenum');
taxis = NaT(size(model_runtime));
for i_time=1:length(taxis)
    taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
end


%% bin CTD profiles according to the fjord boxes
if ~isempty(casts_time_series)
H_obs_axis = NaN([n_layers, length(casts_time_series)]);
for i=1:n_layers
    H_obs_axis(i,:) = interp1(taxis,fjord.s.H(i,:),[casts_time_series.time],"linear","extrap");
end
t_obs = zeros([n_layers,length(casts_time_series)]);
s_obs = zeros([n_layers,length(casts_time_series)]);
if length(casts_time_series)>1
    for k=1:length(casts_time_series)        
        [t_obs(:,k),s_obs(:,k)] = bin_ocean_profiles(casts_time_series(k).temp',casts_time_series(k).sal',casts_time_series(k).depth',H_obs_axis(:,k)',fjord.p);    
    end
else
    [t_obs,s_obs] = bin_ocean_profiles(casts_time_series.temp',casts_time_series.sal',casts_time_series.depth',H_obs_axis,fjord.p);    
end
% filter out weird extrapolation
t_obs(abs(t_obs) > 1e2) = NaN;
s_obs(abs(s_obs) > 1e2) = NaN;
else
    t_obs = []; s_obs = [];
end
%% Plotting figure
% time_lims = [min(min([casts_time_series.time]),min(taxis)),max(max([casts_time_series.time]),max(taxis))];
time_lims = [min(taxis) max(taxis)];

fig_layer_handles.hf = figure();
i_plt=1;
for i_layer=1:n_layers
    subplot(4,2,i_plt); hold on; box on;
    fig_layer_handles.hfj = plot(taxis,fjord.s.T(i_layer,:),'linewidth',1.5,'color',[0.5 0.5 0.5]);
    fig_layer_handles.hsh = plot(taxis,Ts(i_layer,:),'LineStyle','--','linewidth',1.,'color',[0.5 0.5 0.5]);
    if ~isempty(t_obs)
    if size(t_obs,2) > 1
        scatter([casts_time_series.time],t_obs(i_layer,:),30,'black','filled','MarkerFaceAlpha',0.75);
    else
        scatter(casts_time_series.time,t_obs(i_layer),30,'black','filled','MarkerFaceAlpha',0.75);
    end
    end
    ylabel('T (^oC)'); 
    xlim(time_lims); ylim(tlims(i_layer,:))
    text(0.02,1.1,sprintf('Layer %d',i_layer),'Units','normalized')
    text(0.98,0.98,sprintf('(%s)',letters(i_plt)),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
    if i_layer~=n_layers
        set(gca,'xticklabels',{}); 
    else
        xlabel('Time')        
    end

    subplot(4,2,i_plt+1); hold on; box on;
    plot(taxis,fjord.s.S(i_layer,:),'linewidth',1.5,'color',[0.5 0.5 0.5])
    plot(taxis,Ss(i_layer,:),'LineStyle','--','linewidth',1.,'color',[0.5 0.5 0.5])
    if ~isempty(s_obs)
    if size(s_obs,2) > 1
        scatter([casts_time_series.time],s_obs(i_layer,:),30,'black','filled','MarkerFaceAlpha',0.75);
    else
        scatter(casts_time_series.time,s_obs(i_layer),30,'black','filled','MarkerFaceAlpha',0.75);
    end
    end
    text(0.98,0.98,sprintf('(%s)',letters(i_plt+1)),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
    ylabel('S'); 
    xlim(time_lims); ylim(slims(i_layer,:))
    if i_layer~=n_layers
        set(gca,'xticklabels',{}); 
    else
        xlabel('Time')
        fig_layer_handles.hl = legend([fig_layer_handles.hfj; fig_layer_handles.hsh],{'fjord','shelf'},'Location','southwest');
        fig_layer_handles.hl.NumColumns=2;
    end
    
    i_plt=i_plt+2;
end % end for i_plt=n_layers

htitle=axes(fig_layer_handles.hf,'visible','off'); 
htitle.Title.Visible='on';
title(htitle,fjord_title,'interpreter','none')

end