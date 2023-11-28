function hf = plot_fjord_scatter(fjords)

qsg = NaN(size(fjords));
qsgmax = NaN(size(fjords));
vol = NaN(size(fjords));
reg = NaN(size(fjords));
for i_fjord=1:length(fjords)
    qsg(i_fjord) = mean(fjords(i_fjord).f.Qsg,'omitnan');
    qsgmax(i_fjord) = max(fjords(i_fjord).f.Qsg,[],'omitnan');
    vol(i_fjord) = fjords(i_fjord).p.L .* fjords(i_fjord).p.W .* fjords(i_fjord).p.H;
    reg(i_fjord) = fjords(i_fjord).m.regionID;
end

region_line_color = lines(7);
hf = figure('Name','Discharge versus volume','Position',[20 20 800 400]); 
subplot(1,2,1); hold on; box on; grid on;
for i=1:length(region_line_color)
    scatter(1e-9.*vol(i),qsg(i),50,'filled','MarkerFaceColor',region_line_color(i,:))
end
for i=1:length(qsg)
    scatter(1e-9.*vol(i),qsg(i),50,'filled','MarkerFaceColor',region_line_color(reg(i),:))
end
xlabel('Fjord volume (km^3)'); ylabel('Mean subglacial discharge (m^3s^{-1})')
legend('SW','SE','CW','CE','NW','NE','NO','Location','northeast')
text(0.05,0.95,'(a)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
subplot(1,2,2); hold on; box on; grid on;
for i=1:length(qsg)
    scatter(1e-9.*vol(i),qsgmax(i),50,'filled','MarkerFaceColor',region_line_color(reg(i),:))
end
xlabel('Fjord volume (km^3)'); ylabel('Maximum subglacial discharge (m^3s^{-1})')
text(0.05,0.95,'(b)','fontsize',14,'units','normalized')
set(gca,'fontsize',14)
end