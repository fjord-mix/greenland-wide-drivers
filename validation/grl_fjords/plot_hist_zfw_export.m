function [hf_zfw] = plot_hist_zfw_export(ensemble_yr,i_tgt_day)

n_years = length(ensemble_yr);
lcolor=lines(n_years);
h_hist = [];
lbl_hist = {};
zfw_stacked = [];
znb_stacked = [];

hf_zfw = figure; 
ht_zfw = tiledlayout('horizontal');
for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    zfw = NaN(size(ensemble));
    znb = NaN(size(ensemble));
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s)
                % depths at target day
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_tgt_day);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb(i_tgt_day);

                % % average depth over the last year
                zfw(i_fjord,i_run) = mean(ensemble(i_fjord,i_run).s.z_max_export_t(end-365:end));
                znb(i_fjord,i_run) = mean(ensemble(i_fjord,i_run).s.znb_t(end-365:end));
                % 
                % % depth of export maximum over the last year
                % [~,i_fz] = min(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export_t(end-365+i_fz);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb_t(end-365+i_fz);

                zfw_stacked = [zfw_stacked; zfw(:)];
                znb_stacked = [znb_stacked; znb(:)];
            end
        end
    end
    nexttile(1); hold on; box on; grid on;
    title('Neutral buoyancy')
    histogram(znb(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca, 'xdir', 'reverse','fontsize',16);
    ylim([-450 0])

    nexttile(2); hold on; box on; grid on;
    title('Maximum FW export')
    h1 = histogram(zfw(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca,'fontsize',16)
    ylim([-450 0])

    h_sit = [h_hist h1];
    lbl_hist{end+1} = num2str(2015+i_year);
end
nexttile(1); hold on;
histogram(znb_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);
nexttile(2); hold on;
histogram(zfw_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);

legend(h_hist,lbl_hist,'location','best')
ylabel(ht_zfw,'Depth (m)','fontsize',16)
xlabel(ht_zfw,'Probability','fontsize',16)

end