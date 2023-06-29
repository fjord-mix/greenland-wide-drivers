function fig_handles = plot_fw_out(fjord)

if ~isfield(fjord,'o')
    disp('No post-processed data structure found.')
    disp('Please run postprocess_boxmodel() before calling this function!')
else    
    model_runtime = fjord.s.t(1:size(fjord.s.H,2));
    runtime_axis = fjord.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis = NaT(size(model_runtime));
    for i_time=1:length(taxis)
        taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
    end    
    avg_depth=mean(fjord.o.Ze,1);
    % TODO: integrate volume over depth intervals

    hf = figure; hold on
    plot(runtime_axis,fjord.o.Ze(1,:),':k')
    plot(runtime_axis,fjord.o.Ze(2,:),':k')
    hp=surface([runtime_axis';runtime_axis'],[avg_depth;avg_depth],...
               [zeros(size(avg_depth));zeros(size(avg_depth))],...
               [fjord.o.Ve';fjord.o.Ve'],...
                'facecol','no','edgecol','interp','linew',5);
    hc=colorbar;
    ylabel(hc,'Volume flux (m^3 s^{-1})')
    ylabel('Depth (m)')
    xlabel('Time')    

    fig_handles.hf=hf;
    fig_handles.hp=hp;
    fig_handles.hc=hc;
end

end