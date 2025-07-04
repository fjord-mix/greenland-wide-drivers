function [rmse_table,tf_best,sf_best] = get_best_profiles_rmse(res_box,i_fjord,n_runs,w_rmse_t,i_tgt_day,tgt_days)

if nargin < 4 || isempty(w_rmse_t)
    w_rmse_t = 0.5;
end
if nargin < 5
    i_tgt_day=1;
end

if nargin > 5
    z_rmse_t  = normalize(res_box(i_fjord).rmse_tf(:,2),"range");
    z_rmse_s  = normalize(res_box(i_fjord).rmse_sf(:,2),"range");
    rmse_both = w_rmse_t*z_rmse_t + (1-w_rmse_t)*z_rmse_s;

    [rmse_table.tf_rpm,i_min_rmse_tf] = min(res_box(i_fjord).rmse_tf,[],'all','omitnan');
    [rmse_table.sf_rpm,i_min_rmse_sf] = min(res_box(i_fjord).rmse_sf,[],'all','omitnan');
    [rmse_table.ts_rpm,i_min_rmse]    = min(rmse_both,[],'all','omitnan');
    [rmse_table.df_rpm,i_min_rmse_df] = min(res_box(i_fjord).rmse_df,[],'all','omitnan');

    [irun_best_tf,id_best_tf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_tf);
    [irun_best_sf,id_best_sf] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_sf);
    [irun_best_df,id_best_df] = ind2sub([n_runs,length(tgt_days)],i_min_rmse_df);
    [irun_best,id_best] = ind2sub([n_runs,length(tgt_days)],i_min_rmse);

    % tf_best = res_box(i_fjord).ensemble_tf(:,irun_best_tf,id_best_tf); % using RMSE(T)
    % sf_best = res_box(i_fjord).ensemble_sf(:,irun_best_sf,id_best_sf); % using RMSE(S)

    %  % using RMSE(sigma)
    % tf_best = res_box(i_fjord).ensemble_tf(:,irun_best_df,id_best_df);
    % sf_best = res_box(i_fjord).ensemble_sf(:,irun_best_df,id_best_df);

     % using normalised RMSE(both)
    tf_best = res_box(i_fjord).ensemble_tf(:,irun_best,id_best);
    sf_best = res_box(i_fjord).ensemble_sf(:,irun_best,id_best);

    % tf_best2 = res_box(i_fjord).ensemble_tf(:,irun_best,id_best);
    % sf_best2 = res_box(i_fjord).ensemble_sf(:,irun_best,id_best);
    % 
    % inds_best_tf = [irun_best_tf,id_best_tf];
    % inds_best_sf = [irun_best_sf,id_best_sf];
    % inds_best2   = [irun_best,id_best];
else
    z_rmse_t  = normalize(res_box(i_fjord).rmse_tf(:,2),"range");
    z_rmse_s  = normalize(res_box(i_fjord).rmse_sf(:,2),"range");
    rmse_both = w_rmse_t*z_rmse_t + (1-w_rmse_t)*z_rmse_s;

    [rmse_table.tf_rpm,inds_best_tf] = min(squeeze(res_box(i_fjord).rmse_tf(:,i_tgt_day)),[],'all','omitnan');
    [rmse_table.sf_rpm,inds_best_sf] = min(squeeze(res_box(i_fjord).rmse_sf(:,i_tgt_day)),[],'all','omitnan');
    [rmse_table.ts_rpm,inds_best]    = min(squeeze(rmse_both(:,i_tgt_day)),[],'all','omitnan');
    [rmse_table.df_rpm,inds_best_df] = min(squeeze(res_box(i_fjord).rmse_df(:,i_tgt_day)),[],'all','omitnan');

    % tf_best = res_box(i_fjord).ensemble_tf(:,inds_best_tf,i_tgt_day); % using RMSE(T)
    % sf_best = res_box(i_fjord).ensemble_sf(:,inds_best_sf,i_tgt_day); % using RMSE(S)

    % using RMSE(sigma)
    % tf_best = res_box(i_fjord).ensemble_tf(:,inds_best_df,i_tgt_day); 
    % sf_best = res_box(i_fjord).ensemble_sf(:,inds_best_df,i_tgt_day); 

    % using normalised RMSE(both)
    tf_best = res_box(i_fjord).ensemble_tf(:,inds_best,i_tgt_day); 
    sf_best = res_box(i_fjord).ensemble_sf(:,inds_best,i_tgt_day); 

    % tf_best2 = res_box(i_fjord).ensemble_tf(:,inds_best2,i_tgt_day);
    % sf_best2 = res_box(i_fjord).ensemble_sf(:,inds_best2,i_tgt_day);
end

% adds Shelf-fjord RMSE to the table
rmse_table.rmse_ts = res_box(i_fjord).rmse_ts;
rmse_table.rmse_ss = res_box(i_fjord).rmse_ss;
rmse_table.rmse_ds = res_box(i_fjord).rmse_ds;

end