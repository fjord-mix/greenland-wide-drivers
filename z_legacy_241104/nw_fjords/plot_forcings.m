function hf = plot_forcings(f,p,t)

hf = figure('Position',[100 100 1000 400]);    
subplot(1,3,1); hold on; box on;
plot(f.Ts(:,1),f.zs);
scatter(mean(f.Ts(:,1)),-p.Hgl,40,'v','filled','MarkerFaceColor','black')
plot([mean(f.Ts(:,1)) mean(f.Ts(:,1))],[-p.H -p.Hsill],'-k','linewidth',2)
ylim([-p.H 0])
ylabel('Depth (m)')
xlabel('Temperature (^oC)')
set(gca,'fontsize',14)

subplot(1,3,2); hold on; box on;
plot(f.Ss(:,1),f.zs);
scatter(mean(f.Ss(:,1)),-p.Hgl,40,'v','filled','MarkerFaceColor','black')
plot([mean(f.Ss(:,1)) mean(f.Ss(:,1))],[-p.H -p.Hsill],'-k','linewidth',2)
ylim([-p.H 0])
xlabel('Salinity (^oC)')
set(gca,'fontsize',14)

subplot(1,3,3); hold on; box on; grid on;
plot(t,f.Qsg);
ylabel('Subglacial discharge (m^3s^{-1})')
xlabel('Time (days)')
set(gca,'fontsize',14)


% exportgraphics(gcf,[figs_path,'forcing_summary_c13_fjord_',letters{i_fjord},'_',num2str(2015+i_yr),'.png'],'Resolution',300)
% 
end