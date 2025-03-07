function hf = plot_scatter_tempsalt(fjord_model_yr,ensemble_yr,res_box_yr,res_obs_yr,n_runs,tgt_days,i_tgt_day,which_fjords)


n_years      = length(fjord_model_yr);
w_rmse_t     = 0.5; % how much we want to weight the temperature (n)RMSE versus salinity (0.5 = 50:50; 1 = only temperature)

fsize   = 14;
lcolor  = lines(3+length(tgt_days));
letters = lower(char(65:65+length(which_fjords)*2));

fig_width = 900; 
fig_height = 300*length(which_fjords);

temp_bnds = [20, -10];
salt_bnds = [40, 0];
first_time = 1;

all_r_temp = {};
all_r_salt = {};

hf    = figure('Name','Model-obs fit plots','Position',[40 40 900 500]);
%ht     = tiledlayout(1,2,'TileSpacing','loose','padding','compact');
ht    = tiledlayout('flow');
for i_year=n_years:-1:1
    n_fjord_runs = length(fjord_model_yr{i_year});
    res_box      = res_box_yr{i_year};
    res_obs      = res_obs_yr{i_year};
    for i_fjord=1:n_fjord_runs
        
        % Getting the profiles for comparing
        [~,tf_best,sf_best] = get_best_profiles_rmse(res_box,i_fjord,n_runs,w_rmse_t,tgt_days,i_tgt_day);


        % interpolaring observations to box model depths
        tf_obs = interp1(res_obs(i_fjord).zf,res_obs(i_fjord).tf,res_box(i_fjord).zf,'linear');
        sf_obs = interp1(res_obs(i_fjord).zf,res_obs(i_fjord).sf,res_box(i_fjord).zf,'linear');

        if first_time
            nexttile(1); hold on; box on; grid on
            text(0.02,0.98,"(a)",'Units','normalized','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
            xlabel('T_{obs}') ; ylabel('T_{rpm}')
            plot([-10, 10],[-10, 10],'-r')

            nexttile(2); hold on; box on; grid on
            text(0.02,0.98,"(b)",'Units','normalized','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
            xlabel('S_{obs}') ; ylabel('S_{rpm}')
            plot([10, 36],[10, 36],'-r')

            first_time = 0;

        end

        nexttile(1); hold on; box on; grid on
        scatter(tf_obs,tf_best,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);
        
        nexttile(2); hold on; box on; grid on
        scatter(sf_obs,sf_best,40,-res_box(i_fjord).zf,'filled','markerfacealpha',0.5,'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.15);

        all_r_temp{end+1} = corr(tf_obs,tf_best,'type','pearson','rows','complete');
        all_r_salt{end+1} = corr(sf_obs,sf_best,'type','pearson','rows','complete');

        temp_bnds(1) = min(temp_bnds(1),min([tf_obs(:);tf_best(:)]));
        temp_bnds(2) = max(temp_bnds(2),max([tf_obs(:);tf_best(:)]));
        salt_bnds(1) = min(salt_bnds(1),min([sf_obs(:);sf_best(:)]));
        salt_bnds(2) = max(salt_bnds(2),max([sf_obs(:);sf_best(:)]));
        
        
    end % fjords
end % years

nexttile(1);
%xlim(temp_bnds); ylim(temp_bnds);
xlim([-2, 6]); ylim([-2, 6])
hcb = colorbar;
ylabel(hcb,'Depth (m)','FontSize',fsize)
set(gca,'fontsize',fsize)



nexttile(2);
%xlim(salt_bnds); ylim(salt_bnds);
xlim([25, 36]); ylim([25, 36])
set(gca,'fontsize',fsize)
colormap(flip(cmocean('deep')))

nexttile(3)
bin_edges = -0.1:0.05:1.0;
hold on; box on
histogram(cell2mat(all_r_temp),bin_edges,'Normalization','count','FaceAlpha',0.5);
histogram(cell2mat(all_r_salt),bin_edges,'Normalization','count','FaceAlpha',0.5);
legend('Temperature','Salinity')
xlabel('Correlation'); ylabel('Count')
set(gca,'fontsize',fsize)



end