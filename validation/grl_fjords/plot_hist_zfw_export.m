function hf_zfw = plot_hist_zfw_export(ensemble_yr,res_box_yr)

rmse_threshold = 0.5;
n_years = length(ensemble_yr);
lcolor=lines(n_years);
h_hist = [];
lbl_hist = {};
zfw_stacked = [];
znb_stacked = [];
tfw_stacked = [];
t_qsg_stacked = [];

hf_zfw = figure('Position',[50 50, 1200 800]); 
ht_main = tiledlayout(2,3,'TileSpacing','Compact');

for i_year=n_years:-1:1
    ensemble = ensemble_yr{i_year};
    res_box  = res_box_yr{i_year};
    zfw = NaN(size(ensemble));
    znb = NaN(size(ensemble));
    t_fw_max = NaN(size(ensemble));
    t_qsg_max = NaN([size(ensemble,1),1]);
    for i_fjord=1:size(ensemble,1)
        for i_run=1:size(ensemble,2)
            if ~isempty(ensemble(i_fjord,i_run).s) && res_box(i_fjord).rmse_tf(i_run,2) < rmse_threshold
                % depths at target day
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export(i_tgt_day);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb(i_tgt_day);

                % % average depth over the last year
                % zfw(i_fjord,i_run) =
                % mean(ensemble(i_fjord,i_run).s.z_max_export_t(end-365:end)); % average of peak FW exp depth
                zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export;   % depth of the average peak FW exp
                znb(i_fjord,i_run) = mean(ensemble(i_fjord,i_run).s.znb_t(end-365:end));
                % 
                % % depth of export maximum over the last year
                [~,i_fz] = min(ensemble(i_fjord,i_run).s.fw_export_t(end-365:end));
                % zfw(i_fjord,i_run) = ensemble(i_fjord,i_run).s.z_max_export_t(end-365+i_fz);
                % znb(i_fjord,i_run) = ensemble(i_fjord,i_run).s.znb_t(end-365+i_fz);

                t_fw_max(i_fjord,i_run) = i_fz;
                [~,t_qsg_max(i_fjord)] = max(ensemble(i_fjord,i_run).s.Qsg(end-365:end));

                zfw_stacked = [zfw_stacked; zfw(:)];
                znb_stacked = [znb_stacked; znb(:)];
                tfw_stacked = [tfw_stacked; t_fw_max(:)];
            end
        end
    end
    t_qsg_stacked = [t_qsg_stacked; t_qsg_max(:)];

    nexttile(1,[1 1]); hold on; box on; grid on;
    histogram(t_fw_max(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','vertical','FaceAlpha',0.5);
    set(gca,'fontsize',16)
    
    nexttile(4,[1 1]); hold on; box on; grid on;
    histogram(t_qsg_max(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','vertical','FaceAlpha',0.5);
    set(gca,'fontsize',16,'ydir','reverse')

    nexttile(2,[2 1]); hold on; box on; grid on;
    histogram(znb(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca, 'xdir', 'reverse','fontsize',16);
    ylim([-450 0])

    nexttile(3,[2 1]); hold on; box on; grid on;
    h1 = histogram(zfw(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',lcolor(i_year,:),'orientation','horizontal','FaceAlpha',0.5);
    set(gca,'fontsize',16)
    ylim([-450 0])

    h_sit = [h_hist h1];
    lbl_hist{end+1} = num2str(2015+i_year);
end
nexttile(1,[1 1]); hold on;
title('(a) Maximum FW export','fontsize',14)
histogram(tfw_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','vertical','FaceAlpha',1.);
xlim([100 300])
ylabel('Probability','fontsize',16)

nexttile(4,[1 1]); hold on;
title('(b) Maximum subglacial discharge','fontsize',14)
histogram(t_qsg_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','vertical','FaceAlpha',1.);
xlim([100 300])
ylabel('Probability','fontsize',16)
xlabel('Day of year','fontsize',16)

nexttile(2,[2 1]); hold on;
title('(c) Neutral buoyancy','fontsize',14)
histogram(znb_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);
xlabel('Probability','fontsize',16)
ylabel('Depth (m)','fontsize',16)

nexttile(3,[2 1]); hold on;
title('(d) Maximum FW export','fontsize',14)
histogram(zfw_stacked(:),ensemble(1,1).p.N,'Normalization','probability','FaceColor',[0 0 0],'orientation','horizontal','FaceAlpha',1.);
xlabel('Probability','fontsize',16)
legend(h_hist,lbl_hist,'location','southeast')

end