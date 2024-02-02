function plot_reg_ts(datasets,fjords_compilation)
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
[salt_forcing, ~, ~, ~] = get_var_forcing_by_region(fjords_processed,'Ss');
fjord_stats = print_fjord_statistics(fjords_compilation);

temp_forcing = squeeze(mean(temp_forcing,1));
salt_forcing = squeeze(mean(salt_forcing,1));

figure('Name','TS diagrams','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg); hold on; box on;
    zs_prctiles = prctile(fjord_stats.Zs.total,[25, 75]);
    zg_prctiles = prctile(fjord_stats.Zg.total,[25, 75]);
    i_above_sill = depths > zs_prctiles(1);
    i_above_plume = depths > zg_prctiles(1);
    i_below_sill = depths < zs_prctiles(1);
    i_active = i_above_sill & i_above_plume;
    i_present = i_above_sill & ~i_above_plume;

    xlim([min(salt_forcing(:,i_reg))-0.1 max(salt_forcing(:,i_reg))+0.1])
    ylim([min(temp_forcing(:,i_reg))-0.1 max(temp_forcing(:,i_reg))+0.1])

    [thetai,si,dens] = get_sigma_curves(temp_forcing(:,i_reg),salt_forcing(:,i_reg));
    % [c,h]=contour(si,thetai,dens,25:0.2:28,'k');
    [c,h]=contour(si,thetai,dens,'k');
    
    hb = scatter(salt_forcing(:,i_reg),temp_forcing(:,i_reg),40,depths,'filled','^','MarkerEdgeColor','k');
    hp = scatter(salt_forcing(i_present,i_reg),temp_forcing(i_present,i_reg),60,depths(i_present),'filled','s','MarkerEdgeColor','k');
    ha = scatter(salt_forcing(i_active,i_reg),temp_forcing(i_active,i_reg),60,depths(i_active),'filled','o','MarkerEdgeColor','k');
    
    clabel(c,h,'LabelSpacing',100);
    if i_reg > 3
        xlabel('Salinity','FontSize',14);
    end
    if i_reg==1 || i_reg==5
        ylabel('Temperature (^oC)','FontSize',14)
    end
    clim([-1e3 0]);
    set(gca,'fontsize',12)
    text(0.05,1.07,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
hl = legend([hb,hp,ha],{'Below sill','Present','Active'});
colormap(flip(cmocean('deep')))
hc = colorbar('fontsize',14);
ylabel(hc,'Depth (m)');
hl.Position(1)=hl.Position(1)+0.17;
hc.Position(1)=hc.Position(1)+0.17;


function [thetai,si, dens] = get_sigma_curves(theta,s)
    %% generating background density contours
    theta=theta(:);
    s=s(:);
    smin=min(s)-0.01.*min(s);
    smax=max(s)+0.01.*max(s);
    thetamin=min(theta)-0.1*max(theta);
    thetamax=max(theta)+0.1*max(theta);
    xdim=round((smax-smin)./0.1+1);
    ydim=round((thetamax-thetamin)+1);
    dens=zeros(ydim,xdim);
    thetai=((1:ydim)-1)*1+thetamin;
    si=((1:xdim)-1)*0.1+smin;
    disp(xdim);disp(ydim);
    % something here is outputting to the console...
    for j=1:ydim
        for k=1:xdim
            dens(j,k) = gsw_sigma0(si(k),thetai(j)); 
        end
    end
end

end