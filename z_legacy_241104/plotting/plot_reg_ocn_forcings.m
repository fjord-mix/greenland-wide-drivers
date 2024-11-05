function plot_reg_ocn_forcings(datasets,fjords_compilation,temp_amp_seasonal)
letters = {'a','b','c','d','e','f','g','h'};
regions_lbl = {'SW','SE','CW','CE','NW','NE','NO'};

if nargin < 3 || isempty(temp_amp_seasonal)
    temp_amp_seasonal=0.5; % temperature amplitude threshold to be considered as the "seasonal layer"
end

datasets.opts.time_start = datetime(2010,01,15);
datasets.opts.time_end   = datetime(2018,12,15);
datasets.opts.time_interval = [datasets.opts.time_start,datasets.opts.time_end]; 
datasets.opts.dt            = 30.0; % time step (in days) for creating the forcings
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
end

[temp_forcing, ~, ~, depths] = get_var_forcing_by_region(fjords_processed,'Ts');
[salt_forcing, ~, ~, ~] = get_var_forcing_by_region(fjords_processed,'Ss');

% fjord_rho = (fjords_processed(1).p.betaS*salt_forcing - fjords_processed(1).p.betaT*temp_forcing);
% sc_reg = squeeze(trapz(depths,salt_forcing.* fjord_rho,2)./max(abs(depths)));
% squeeze(trapz(depths,(temp_forcing+273.15).* fjord_rho,2)./max(abs(depths)).* fjords_processed(1).p.cw);

% time_axis = datasets.opts.time_start:datasets.opts.dt:datasets.opts.time_end;
time_axis     = fjords_processed(1).t;

depth_seasonal=NaN([1,size(temp_forcing,3)]);
for i_reg=1:length(regions_lbl)    
    temp_amp_seasonal = std(detrend(temp_forcing(:,end,i_reg)));
    for i_z=length(depths):-1:1 % will go down the profile until it finds the deepest amplitude
        temp_detrend=detrend(temp_forcing(:,i_z,i_reg));
        % if (max(temp_detrend)-min(temp_detrend)) > temp_amp_seasonal
        if std(temp_detrend) > 0.15*temp_amp_seasonal
            depth_seasonal(i_reg) = depths(i_z);
        end
    end
end

%% Plotting Temperature and salinity
figure('Name','Temperature undisturbed','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg)
    imagesc(time_axis,-depths,temp_forcing(:,:,i_reg)');
    % yline(-depth_seasonal(i_reg),'--r','linewidth',2)
    clim([-2.5 7]);
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('time (days)'); end
    set(gca,'fontsize',14)
    % p = polyfit(time_axis,hc_reg(:,i_reg),1);
    % trend_ohc = p(1)*1e3;
    text(0.05,1.07,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
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
    % p = polyfit(time_axis,sc_reg(:,i_reg),1);
    % trend_osc = p(1)*1e3;

    % text(0.05,1.07,['(',letters{i_reg},') ',regions_lbl{i_reg}],'units','normalized','fontsize',14)
    text(0.05,1.07,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
colormap(cmocean('haline'))
hc = colorbar('fontsize',14);
ylabel(hc,'Salinity');
hc.Position(1)=hc.Position(1)+0.15;


end