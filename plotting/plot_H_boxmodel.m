function fig_H_handles = plot_H_boxmodel(fjord_title,fjord,model_runtime)


H_series = cumsum(fjord.s.H,1); % get the actual depths of each layer


%% Determining plot parameters
n_layers = size(fjord.s.H,1);

runtime_axis = fjord.m.time_axis;
t0 = convertTo(runtime_axis(1),'datenum');
taxis = NaT(size(model_runtime));
for i_time=1:length(taxis)
    taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
end

time_lims = [min(taxis) max(taxis)];
depth_lims = [min(-H_series(:)), 0];
%% Plot figure
fig_H_handles.hf = figure();
hold on; box on;
for i_layer=1:n_layers
    plot(taxis,-1*H_series(i_layer,:),'color',[0 0 0],'linewidth',1.5);
end
scatter(taxis(1),fjord.p.zgl,50,'black','filled','o'); 
text(taxis(1),fjord.p.zgl,'  grounding line','fontsize',12,'color',[0 0 0],'HorizontalAlignment','left','VerticalAlignment','bottom');

scatter(taxis(1),-1*H_series(end-1,1),50,'black','filled','v')
text(taxis(1),-1*H_series(end-1,1),' sill','fontsize',12,'color',[0 0 0],'HorizontalAlignment','left','VerticalAlignment','top');

xlabel('Time');      xlim(time_lims);
ylabel('Depth (m)'); ylim(depth_lims);
title(fjord_title)
end