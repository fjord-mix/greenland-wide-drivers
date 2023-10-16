function plot_forcings_summary(fjord)
    p = fjord.p;
    a = fjord.a;
    f = fjord.f;
    t = fjord.t;
    
    figure('Position',[100 100 1000 500]);
    subplot(1,4,1), plot(f.Ts(:,1),f.zs); xlabel('T_0'); ylabel('Depth');
    yline(-cumsum(a.H0),':k','linewidth',0.5); ylim([-sum(a.H0) 0])
    text(0.02,(p.H+p.zgl)/p.H,'grounding line','units','normalized')
    text(0.98,(p.H+p.silldepth)/p.H,'sill','units','normalized','HorizontalAlignment','right')
    
    subplot(1,4,2), plot(f.Ss(:,1),f.zs); xlabel('S_0');
    yline(-cumsum(a.H0),':k','linewidth',0.5); ylim([-sum(a.H0) 0])
    text(0.02,0.08,sprintf('Fjord length: %.2f km',p.L*1e-3),'units','normalized');
    text(0.02,0.05,sprintf('Fjord width: %.2f km',p.W*1e-3),'units','normalized');
    subplot(2,4,3), plot(t,f.Ts(end,:)); xlabel('time'); ylabel('f.Ts at surface');
    subplot(2,4,4), plot(t,f.Ss(end,:)); xlabel('time'); ylabel('f.Ss at surface');
    subplot(2,4,7), plot(t,f.Qsg); xlabel('time'); ylabel('Qsg');
    subplot(2,4,8), plot(t,f.D); xlabel('time'); ylabel('D');

end