figure('Position',[100 100 1000 500]);
ints=[0,-cumsum(a.H0)];
subplot(1,4,1); hold on; box on;
plot(f.Ts(:,1),f.zs); 
xlabel('T_0'); ylabel('Depth');
scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');

yline(-cumsum(a.H0),':k','linewidth',0.5); ylim([-sum(a.H0) 0])
text(0.02,(p.H+p.zgl)/p.H,'grounding line','units','normalized')
text(0.98,(p.H+p.silldepth)/p.H,'sill','units','normalized','HorizontalAlignment','right')

subplot(1,4,2); hold on; box on
plot(f.Ss(:,1),f.zs); 
scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
xlabel('S_0');
yline(-cumsum(a.H0),':k','linewidth',0.5); ylim([-sum(a.H0) 0])
text(0.02,0.08,sprintf('Fjord length: %.2f km',p.L*1e-3),'units','normalized');
text(0.02,0.05,sprintf('Fjord width: %.2f km',p.W*1e-3),'units','normalized');
% subplot(2,4,3), plot(Parameters.t,f.Ts(end,:)); xlabel('time'); ylabel('f.Ts at surface');
% subplot(2,4,4), plot(Parameters.t,f.Ss(end,:)); xlabel('time'); ylabel('f.Ss at surface');
subplot(2,4,3), imagesc(Parameters.t,-f.zs,f.Ts); xlabel('time'); ylabel('depth (m)'); title('f.Ts');
subplot(2,4,4), imagesc(Parameters.t,-f.zs,f.Ss); xlabel('time'); title('f.Ss');
subplot(2,4,7), plot(Parameters.t,f.Qsg); xlabel('time'); ylabel('Qsg');
subplot(2,4,8), plot(Parameters.t,f.D); xlabel('time'); ylabel('D');