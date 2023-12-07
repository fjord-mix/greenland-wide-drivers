function plot_reg_ocn_profiles(datasets,fjords_compilation)
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

% time_axis = datasets.opts.time_start:datasets.opts.dt:datasets.opts.time_end;
region_line_color = lines(7);

%% Plotting Temperature
figure('Name','Temperature profiles','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    temp_bnds = prctile(temp_forcing(:,:,i_reg),[5 95],1);
    temp_median_profile = median(temp_forcing(:,:,i_reg),1,'omitnan');
    % temp_lower_profile = temp_bnds(1,:);
    % temp_upper_profile = temp_bnds(2,:);
    temp_lower_profile = temp_median_profile-std(temp_forcing(:,:,i_reg),1,'omitnan');
    temp_upper_profile = temp_median_profile+std(temp_forcing(:,:,i_reg),1,'omitnan');
    y2 = [depths, fliplr(depths)];
    inBetween = [temp_lower_profile, fliplr(temp_upper_profile)];

    subplot(2,4,i_reg); hold on; box on; grid on;
    hp = fill(inBetween, y2, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(temp_lower_profile,depths,'color',region_line_color(i_reg,:));
    plot(temp_upper_profile,depths,'color',region_line_color(i_reg,:));
    plot(temp_median_profile,depths,'color',region_line_color(i_reg,:),'LineWidth',2);
    ylim([min(depths) 0])
    xlim([-2 6])
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('Temperature (^oC)'); end
    set(gca,'fontsize',14)
    text(0.03,1.08,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
subplot(2,4,8); hold on; box on; grid on;
hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'color','none','XAxisLocation','top','YAxisLocation','left');
hold(hAx,'on');
pdf_zg = pdf(fjord_stats.Zg.pd,depths); pdf_zg(pdf_zg==0) = NaN;
pdf_zs = pdf(fjord_stats.Zs.pd,depths); pdf_zs(pdf_zs==0) = NaN;

xline(hAx(1),0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
plot(hAx(2),pdf_zg,depths,'linewidth',1.5,'color','k');
plot(hAx(2),-pdf_zs,depths,'linewidth',1.5,'color','k');
plot(hAx(1),cdf(fjord_stats.Zg.pd,depths,'upper'),depths,'linestyle','-','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(hAx(1),-cdf(fjord_stats.Zs.pd,depths,'upper'),depths,'linestyle','-','linewidth',0.5,'color',[0.7 0.7 0.7]);
set(hAx,'Xticklabel',[])
set(hAx(2),'Yticklabel',[],'Xtick',[],'visible','off')
xlabel(hAx(1),'Cumulative probability')
text(0.02,0.2,'Z_{sill}','units','normalized','fontsize',14)
text(0.98,0.2,'Z_{gl}','units','normalized','fontsize',14,'HorizontalAlignment','right')
text(0.03,1.08,sprintf('(%s)',letters{8}),'units','normalized','fontsize',14)
set(hAx,'fontsize',14)
hAx(2).Box='off';
ylim(hAx,[min(depths) 0])
xlim(hAx(2),[-max(max(pdf(fjord_stats.Zg.pd,depths)),max(pdf(fjord_stats.Zs.pd,depths))) max(max(pdf(fjord_stats.Zg.pd,depths)),max(pdf(fjord_stats.Zs.pd,depths)))])
xlim(hAx(1),[-1 1])

%% Plotting Salinity
figure('Name','Salinity profiles','position',[40 40 1000 400])
for i_reg=1:length(regions_lbl)
    subplot(2,4,i_reg); hold on; box on; grid on;
    salt_bnds = prctile(salt_forcing(:,:,i_reg),[5 95],1);
    salt_median_profile = median(salt_forcing(:,:,i_reg),1,'omitnan');
    % salt_lower_profile = salt_bnds(1,:);
    % salt_upper_profile = salt_bnds(2,:);
    salt_lower_profile = salt_median_profile-std(salt_forcing(:,:,i_reg),1,'omitnan');
    salt_upper_profile = salt_median_profile+std(salt_forcing(:,:,i_reg),1,'omitnan');
    y2 = [depths, fliplr(depths)];
    inBetween = [salt_lower_profile, fliplr(salt_upper_profile)];
    
    subplot(2,4,i_reg); hold on; box on; grid on;
    hp = fill(inBetween, y2, region_line_color(i_reg,:),'edgecolor','none','facealpha',0.2);
    plot(salt_lower_profile,depths,'color',region_line_color(i_reg,:));
    plot(salt_upper_profile,depths,'color',region_line_color(i_reg,:));
    plot(salt_median_profile,depths,'color',region_line_color(i_reg,:));
    ylim([min(depths) 0])
    xlim([30 35.5])
    if i_reg==1 || i_reg==5, ylabel('Depth (m)'); end
    if i_reg>3, xlabel('Salinity'); end
    set(gca,'fontsize',14)
    text(0.03,1.08,sprintf('(%s) %s',letters{i_reg},regions_lbl{i_reg}),'units','normalized','fontsize',14)
end
subplot(2,4,8); hold on; box on; grid on;
hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'color','none','XAxisLocation','top','YAxisLocation','left');
hold(hAx,'on');
pdf_zg = pdf(fjord_stats.Zg.pd,depths); pdf_zg(pdf_zg==0) = NaN;
pdf_zs = pdf(fjord_stats.Zs.pd,depths); pdf_zs(pdf_zs==0) = NaN;

xline(hAx(1),0.0,'linewidth',1.5,'linestyle','--','color',[0.5 0.5 0.5]); 
plot(hAx(2),pdf_zg,depths,'linewidth',1.5,'color','k');
plot(hAx(2),-pdf_zs,depths,'linewidth',1.5,'color','k');
plot(hAx(1),cdf(fjord_stats.Zg.pd,depths,'upper'),depths,'linestyle','-','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(hAx(1),-cdf(fjord_stats.Zs.pd,depths,'upper'),depths,'linestyle','-','linewidth',0.5,'color',[0.7 0.7 0.7]);
set(hAx,'Xticklabel',[])
set(hAx(2),'Yticklabel',[],'Xtick',[],'visible','off')
xlabel(hAx(1),'Cumulative probability')
text(0.02,0.2,'Z_{sill}','units','normalized','fontsize',14)
text(0.98,0.2,'Z_{gl}','units','normalized','fontsize',14,'HorizontalAlignment','right')
text(0.03,1.08,sprintf('(%s)',letters{8}),'units','normalized','fontsize',14)
set(hAx,'fontsize',14)
hAx(2).Box='off';
ylim(hAx,[min(depths) 0])
xlim(hAx(2),[-max(max(pdf(fjord_stats.Zg.pd,depths)),max(pdf(fjord_stats.Zs.pd,depths))) max(max(pdf(fjord_stats.Zg.pd,depths)),max(pdf(fjord_stats.Zs.pd,depths)))])
xlim(hAx(1),[-1 1])
end