%% Geometry parameter distributions
figure;
subplot(2,3,1); hold on; box on
plot(fjord_stats.L.xi,fjord_stats.L.kern,'linewidth',1.5);
% plot(fjord_stats.L.xi,pdf(fjord_stats.L.pd,fjord_stats.L.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.L.xi)]);
text(0.05,0.95,'(a)','Units','normalized')
xlabel('Fjord length (m)'); ylabel('Kernel density')

subplot(2,3,2); hold on; box on
plot(fjord_stats.W.xi,fjord_stats.W.kern,'linewidth',1.5);
% plot(fjord_stats.W.xi,pdf(fjord_stats.W.pd,fjord_stats.W.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.W.xi)]);
text(0.05,0.95,'(b)','Units','normalized')
xlabel('Fjord width (m)');

subplot(2,3,3); hold on; box on
plot(fjord_stats.H.xi,fjord_stats.H.kern,'linewidth',1.5);
% plot(fjord_stats.H.xi,pdf(fjord_stats.H.pd,fjord_stats.H.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.W.xi)]);
text(0.05,0.95,'(c)','Units','normalized')
xlabel('Fjord depth (m)');

subplot(2,3,4); hold on; box on
plot(fjord_stats.Zg.xi,fjord_stats.Zg.kern,'linewidth',1.5);
% plot(fjord_stats.Zg.xi,pdf(fjord_stats.Zg.pd,fjord_stats.Zg.xi),':r','linewidth',1.5)
% xlim([min(fjord_stats.Zg.xi) 0]);
text(0.05,0.95,'(d)','Units','normalized')
xlabel('Grounding line depth (m)'); ylabel('Kernel density')

subplot(2,3,5); hold on; box on
plot(fjord_stats.Zs.xi,fjord_stats.Zs.kern,'linewidth',1.5);
% plot(fjord_stats.Zs.xi,pdf(fjord_stats.Zs.pd,fjord_stats.Zs.xi),':r','linewidth',1.5)
% xlim([min(fjord_stats.Zs.xi) 0]);
text(0.05,0.95,'(e)','Units','normalized')
xlabel('Sill depth (m)');


%% Glacier parameter distributions
x_d = linspace(min(d_anom(:,i_reg)), max(d_anom(:,i_reg)),3000);
x_q = linspace(0, 2,3000);
x_p = linspace(5,45,3000);
figure('Name','Discharge-anomalies parameter space'); hold on
subplot(1,3,1); plot(x_d,pdf(d_pd,x_d),'linewidth',1.5); box on
xlabel('D_a','fontsize',14); ylabel('Kernel density','fontsize',14); text(0.05,0.95,'(a)','Units','normalized','fontsize',14)
subplot(1,3,2); plot(x_q,pdf(q_pd,x_q),'linewidth',1.5); box on
xlabel('Q_a','fontsize',14); text(0.05,0.95,'(b)','Units','normalized','fontsize',14)
subplot(1,3,3); plot(x_p,pdf(p_pd,x_p),'linewidth',1.5); box on
xlabel('PW','fontsize',14); text(0.05,0.95,'(c)','Units','normalized','fontsize',14)

%% Ocean parameter distributions
x_t = linspace(min(tocn_anom(:,i_reg)), max(tocn_anom(:,i_reg)),1000);
x_s = linspace(min(socn_anom(:,i_reg)), max(socn_anom(:,i_reg)),1000);
x_w = linspace(0, 0.5,1000);

figure('Name','Ocean-anomalies parameter space'); hold on
subplot(1,3,1); plot(x_t,pdf(tocn_pd,x_t),'linewidth',1.5); box on
xlabel('T_a','fontsize',14); ylabel('Kernel density','fontsize',14); text(0.05,0.95,'(a)','Units','normalized','fontsize',14)
subplot(1,3,2); plot(x_s,pdf(socn_pd,x_s),'linewidth',1.5); box on
xlabel('S_a','fontsize',14); text(0.05,0.95,'(b)','Units','normalized','fontsize',14)
subplot(1,3,3); plot(x_w,pdf(omeg_pd,x_w),'linewidth',1.5); box on
xlabel('\omega','fontsize',14); text(0.05,0.95,'(c)','Units','normalized','fontsize',14)

%% Checking for interdependencies

zgl   = NaN(size(fjords_processed));
zsill = NaN(size(fjords_processed));
zfjord= NaN(size(fjords_processed));
lfjord= NaN(size(fjords_processed));
wfjord= NaN(size(fjords_processed));
for i=1:length(fjords_processed)
    zgl(i)   = fjords_processed(i).p.zgl;
    zsill(i) = fjords_processed(i).p.silldepth;
    zfjord(i)= fjords_processed(i).p.H;
    lfjord(i)= fjords_processed(i).p.L;
    wfjord(i)= fjords_processed(i).p.W;
end
[rho_zgl,p_zgl]     = corrcoef(zfjord,-zgl);
[rho_zsill,p_zsill] = corrcoef(zfjord,-zsill);
[rho_lw,p_lw] = corrcoef(lfjord,wfjord);

figure('Name','Relationship between fjord geometry parameters'); 
subplot(1,3,1); hold on; box on; grid on
scatter(zfjord,-zgl,'.k');
xlabel('Fjord depth (m)'); ylabel('Grounding line depth (m)')
text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_zgl(1,2),p_zgl(1,2)),'HorizontalAlignment','right','Units','normalized')
subplot(1,3,2); hold on; box on; grid on
scatter(zfjord,-zsill,'.k')
xlabel('Fjord depth (m)'); ylabel('Sill depth (m)')
text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_zsill(1,2),p_zsill(1,2)),'HorizontalAlignment','right','Units','normalized')
subplot(1,3,3); hold on; box on; grid on
scatter(1e-3.*lfjord,1e-3.*wfjord,'.k')
xlabel('Fjord length (km)'); ylabel('Fjord width (km)')
text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_lw(1,2),p_lw(1,2)),'HorizontalAlignment','right','Units','normalized')

%clear zgl zsill zfjord lfjord wfjord rho_zgl p_zgl rho_zsill p_zsill rho_lw p_lw