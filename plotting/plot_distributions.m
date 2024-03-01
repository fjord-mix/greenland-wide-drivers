function hf = plot_distributions(datasets,fjords_compilation)
fjord_stats = print_fjord_statistics(fjords_compilation);
figure("Name",'Probability functions','Position',[50 50 1200 700]);

%% Geometry parameter distributions

subplot(3,5,1); hold on; box on
plot(fjord_stats.L.xi,fjord_stats.L.kern,'linewidth',1.5,'color','k');
% plot(fjord_stats.L.xi,pdf(fjord_stats.L.pd,fjord_stats.L.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.L.xi)]);
text(0.95,0.95,'(a)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('L (m)','fontsize',14); ylabel('Probability function','fontsize',14)

subplot(3,5,2); hold on; box on
plot(fjord_stats.W.xi,fjord_stats.W.kern,'linewidth',1.5,'color','k');
% plot(fjord_stats.W.xi,pdf(fjord_stats.W.pd,fjord_stats.W.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.W.xi)]);
text(0.95,0.95,'(b)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('W (m)','fontsize',14);

subplot(3,5,3); hold on; box on
plot(fjord_stats.H.xi,fjord_stats.H.kern,'linewidth',1.5,'color','k');
% plot(fjord_stats.H.xi,pdf(fjord_stats.H.pd,fjord_stats.H.xi),':r','linewidth',1.5)
% xlim([0 max(fjord_stats.H.xi)]);
text(0.95,0.95,'(c)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('H (m)','fontsize',14);

subplot(3,5,4); hold on; box on
plot(fjord_stats.Zg.xi,fjord_stats.Zg.kern,'linewidth',1.5,'color','k');
% plot(fjord_stats.Zg.xi,pdf(fjord_stats.Zg.pd,fjord_stats.Zg.xi),':r','linewidth',1.5)
% xlim([min(fjord_stats.Zg.xi) 0]);
text(0.95,0.95,'(d)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('z_{gl} (m)','fontsize',14);

subplot(3,5,5); hold on; box on
plot(fjord_stats.Zs.xi,fjord_stats.Zs.kern,'linewidth',1.5,'color','k');
% plot(fjord_stats.Zs.xi,pdf(fjord_stats.Zs.pd,fjord_stats.Zs.xi),':r','linewidth',1.5)
% xlim([min(fjord_stats.Zs.xi) 0]);
text(0.95,0.95,'(e)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('z_s (m)','fontsize',14);

subplot(3,5,6); hold on; box on
plot(fjord_stats.a.xi,fjord_stats.a.kern,'linewidth',1.5,'color','k');
text(0.95,0.95,'(f)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('L/W','fontsize',14); % xlabel('W (m)','fontsize',14);
ylabel('Probability function','fontsize',14); 

subplot(3,5,7); hold on; box on
plot(fjord_stats.Zgr.xi,fjord_stats.Zgr.kern,'linewidth',1.5,'color','k');
text(0.95,0.95,'(g)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('z_{gl}/H','fontsize',14); % xlabel('W (m)','fontsize',14);

subplot(3,5,8); hold on; box on
plot(fjord_stats.Zsr.xi,fjord_stats.Zsr.kern,'linewidth',1.5,'color','k');
text(0.95,0.95,'(h)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
xlabel('z_s/H','fontsize',14); % xlabel('W (m)','fontsize',14);

for i_reg=1:7
    [~,~,probs,~] = define_model_param_distrib(datasets,fjords_compilation,i_reg);
    
    %% Ocean parameter distributions
    x_t = linspace(-5, 5,1000);
    x_s = linspace(-5, 5,1000);
    % x_w = linspace(0, 0.6,1000);
    
    % subplot(3,5,9);  hold on; box on; plot(x_w,pdf(probs(8),x_w),'linewidth',1.5,'color','k');
    subplot(3,5,11); hold on; box on; plot(x_t,pdf(probs(6),x_t),'linewidth',1.5);
    subplot(3,5,12); hold on; box on; plot(x_s,pdf(probs(7),x_s),'linewidth',1.5);

    %% Glacier parameter distributions    
    x_d = linspace(-500,500,3000);
    x_q = linspace(0, 2,3000);

    subplot(3,5,9); hold on; box on; plot(x_q,pdf(probs(8),x_q),'linewidth',1.5);      
    subplot(3,5,10); hold on; box on; plot(x_d,pdf(probs(9),x_d),'linewidth',1.5);
end

%% Model parametre distributions
x_p = linspace(0,40,3000);
x_c = linspace(1e2,1e6,10000);
subplot(3,5,13); hold on; box on; plot(x_c,pdf(probs(11),x_c),'linewidth',1.5,'color','k');
subplot(3,5,14); hold on; box on; plot(x_p,pdf(probs(10),x_p),'linewidth',1.5,'color','k');

%% Adding labels and such
% subplot(3,5,9);
% xlabel('\omega (day^{-1})','fontsize',14); text(0.95,0.95,'(i)','Units','normalized','fontsize',14,'HorizontalAlignment','right')

subplot(3,5,9);
xlabel('Q_a','fontsize',14); text(0.95,0.95,'(i)','Units','normalized','fontsize',14,'HorizontalAlignment','right');  

subplot(3,5,10);
xlabel('D_a (m^3s^{-1})','fontsize',14); text(0.95,0.95,'(j)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
ylim([0 0.02]); xlim([-300 500]);


subplot(3,5,11);
xlabel('T_a (^oC)','fontsize',14); text(0.95,0.95,'(k)','Units','normalized','fontsize',14,'HorizontalAlignment','right');
ylim([0 3]); xlim([-2 3]);
ylabel('Probability function','fontsize',14);
hl = legend('SW','SE','CW','CE','NW','NE','NO','fontsize',14); 
hl.Position=[0.8,0.12,0.06,0.19]; %TODO

subplot(3,5,12);
xlabel('S_a (-)','fontsize',14); text(0.95,0.95,'(l)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
ylim([0 6]); xlim([-1.2 1]);

subplot(3,5,13);
xlabel('C0 (s)','fontsize',14); text(0.95,0.95,'(m)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
set(gca,'XScale','log')
xlim([0.5e3 1.5e5]); %ylim([0 6]); 

subplot(3,5,14);
xlabel('PW (m)','fontsize',14); text(0.95,0.95,'(n)','Units','normalized','fontsize',14,'HorizontalAlignment','right')
lbls_p0 = get(gca,'xticklabels');
for i=1:length(lbls_p0),lbls_p0{i} = sprintf('%d',10*str2num(lbls_p0{i})); end
set(gca,'xticklabels',lbls_p0);


hf = gcf;
end
%% Checking for interdependencies

% zgl   = NaN(size(fjords_processed));
% zsill = NaN(size(fjords_processed));
% zfjord= NaN(size(fjords_processed));
% lfjord= NaN(size(fjords_processed));
% wfjord= NaN(size(fjords_processed));
% for i=1:length(fjords_processed)
%     zgl(i)   = fjords_processed(i).p.zgl;
%     zsill(i) = fjords_processed(i).p.silldepth;
%     zfjord(i)= fjords_processed(i).p.H;
%     lfjord(i)= fjords_processed(i).p.L;
%     wfjord(i)= fjords_processed(i).p.W;
% end
% [rho_zgl,p_zgl]     = corrcoef(zfjord,-zgl);
% [rho_zsill,p_zsill] = corrcoef(zfjord,-zsill);
% [rho_lw,p_lw] = corrcoef(lfjord,wfjord);
% 
% figure('Name','Relationship between fjord geometry parameters'); 
% subplot(1,3,1); hold on; box on; grid on
% scatter(zfjord,-zgl,'.k');
% xlabel('Fjord depth (m)'); ylabel('Grounding line depth (m)')
% text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_zgl(1,2),p_zgl(1,2)),'HorizontalAlignment','right','Units','normalized')
% subplot(1,3,2); hold on; box on; grid on
% scatter(zfjord,-zsill,'.k')
% xlabel('Fjord depth (m)'); ylabel('Sill depth (m)')
% text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_zsill(1,2),p_zsill(1,2)),'HorizontalAlignment','right','Units','normalized')
% subplot(1,3,3); hold on; box on; grid on
% scatter(1e-3.*lfjord,1e-3.*wfjord,'.k')
% xlabel('Fjord length (km)'); ylabel('Fjord width (km)')
% text(0.95,0.95,sprintf('R=%.2f\np=%.2f',rho_lw(1,2),p_lw(1,2)),'HorizontalAlignment','right','Units','normalized')

%clear zgl zsill zfjord lfjord wfjord rho_zgl p_zgl rho_zsill p_zsill rho_lw p_lw