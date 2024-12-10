time1=1;
time2=i+1;

figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),s.T(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('Temperature (^oC)')
xlim([t(time1) t(time2)])
end

figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),s.H(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('H (m)')
xlim([t(time1) t(time2)])
end

figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),Q.QVg(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('QVg (m^3/s)')
xlim([t(time1) t(time2)])
end

figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),Q.QVs(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('QVs (m^3/s)')
xlim([t(time1) t(time2)])
end

figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),Q.QVmi(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('QVmi (m^3/s)')
xlim([t(time1) t(time2)])
end


figure('Position',[200 200 1000 500]);
tiledlayout('flow')
for k=1:10:60
nexttile;
plot(t(time1:time2),Q.QVv(k:k+9,time1:time2))
text(0.95,0.95,sprintf("layers %d to %d",k,k+9),'Units','normalized','HorizontalAlignment','right')
xlabel('time (days)')
ylabel('QVv (m^3/s)')
xlim([t(time1) t(time2)])
end