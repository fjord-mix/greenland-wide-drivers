function plot_reg_ocn_forcings(datasets,fjords_compilation)
letters = {'a','b','c','d','e','f','g','h'};
regions_lbl = {'SW','SE','CW','CE','NW','NE','NO'};

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

[temp_forcing, ~, ~, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
[salt_forcing, ~, ~, depths] = get_var_forcing_by_region(fjords_processed,'Ss');


% time_axis = datasets.opts.time_start:datasets.opts.dt:datasets.opts.time_end;
time_axis = fjords_processed(1).t;

%% Plotting figure
figure('Name','Temperature undisturbed','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg)
    imagesc(time_axis,-depths,temp_forcing(:,:,i_reg)');
    clim([-2.5 7]);
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('time (days)'); end
    set(gca,'fontsize',14)

    text(0.05,1.07,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
end
colormap(cmocean('thermal'))
hc = colorbar('fontsize',14);
ylabel(hc,'Temperature (^oC)');
hc.Position(1)=hc.Position(1)+0.15;


figure('Name','Salinity undisturbed','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg)
    imagesc(time_axis,-depths,salt_forcing(:,:,i_reg)');
    clim([30 36]);
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('time (days)'); end
    set(gca,'fontsize',14)

    text(0.05,1.07,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
end
colormap(cmocean('haline'))
hc = colorbar('fontsize',14);
ylabel(hc,'Salinity');
hc.Position(1)=hc.Position(1)+0.15;

end