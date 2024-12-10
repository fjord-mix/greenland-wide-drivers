function hf = plot_crashed_parameters(fjords_crashed,param_names)

hf = figure('Name','Crashed run parameters','Position',[40 40 900 250*length(param_names)]);
ht_p = tiledlayout('flow');
i_run=1;
for i_param=1:length(param_names)
    nexttile; hold on; box on;
    for i_fjord=1:length(fjords_crashed)
        param_name = param_names{i_param};
        eval(['bad_params = fjords_crashed{i_fjord}.p.',param_name,';']);
        plot(i_run,bad_params,'-ok')
        i_run=i_run+1;
    end
    xlabel('Simulation')
    ylabel(param_names{i_param})
end

end