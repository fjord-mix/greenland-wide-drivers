function fig_handles = plot_ts_at_depth(fjord,depths,interp_method)

if ~isfield(fjord,'o')
    disp('No post-processed data structure found.')
    disp('Please run postprocess_boxmodel() before calling this function!')
else
    if nargin==2 || isempty(interp_method)
        interp_method='nearest';
    end
    model_runtime = fjord.s.t(1:size(fjord.s.H,2));
    runtime_axis = fjord.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis = NaT(size(model_runtime));
    for i_time=1:length(taxis)
        taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
    end

    depth_lbls = strsplit(num2str(depths));
    for i=1:length(depth_lbls)
        depth_lbls{i} = [depth_lbls{i},' m'];
    end
    
    t_series=interp1(-fjord.f.zs,fjord.o.Tf,depths,interp_method);
    s_series=interp1(-fjord.f.zs,fjord.o.Sf,depths,interp_method);
    t_series=t_series(:,1:length(model_runtime));
    s_series=s_series(:,1:length(model_runtime));

    hf = figure;
    subplot(2,1,1)    
    fig_handles.hft = plot(taxis,t_series,'linewidth',1.5);
    ylabel('Temperature (^oC)')
    title(fjord.m.name,'interpreter','none')
    subplot(2,1,2)
    fig_handles.hfs = plot(taxis,s_series,'linewidth',1.5);
    ylabel('Salinity (-)')
    xlabel('Time')
    hl=legend(depth_lbls);

end

end