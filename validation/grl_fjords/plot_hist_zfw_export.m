function [hf_zfw] = plot_hist_zfw_export(ensemble_yr,i_tgt_day)

n_years = length(ensemble_yr);
lcolor=lines(n_years);
h_hist = [];
lbl_hist = {};

hf_zfw = figure; 
tiledlayout('horizontal')
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    zfw = NaN(size(ensemble));
    znb = NaN(size(ensemble));
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s)
                zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_tgt_day);
                znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb(i_tgt_day);
            end
        end
    end
    nexttile; hold on; box on; grid on;
    title('Neutral buoyancy')
    histogram(znb(:),ensemble(1,1).p.N,'FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca, 'xdir', 'reverse');

    nexttile; hold on; box on; grid on;
    title('Maximum FW export')
    h1 = histogram(zfw(:),ensemble(1,1).p.N,'FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    
    h_sit = [h_hist h1];
    lbl_hist{end+1} = num2str(2015+i_year);
end

legend(h_hist,lbl_hist,'location','best')
ylabel('Depth (m)')
xlabel('Count')
set(gca,'fontsize',16)

end