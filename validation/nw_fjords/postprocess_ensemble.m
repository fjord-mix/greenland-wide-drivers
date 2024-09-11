function [res_obs,res_box] = postprocess_ensemble(fjord_model,ensemble,tgt_days)


n_fjords = size(ensemble,1);
n_runs   = size(ensemble,2);

res_obs(size(fjord_model)) = struct("tf",[],"sf",[],"zf",[]);
res_box(size(fjord_model)) = struct("tf",[],"sf",[],"zf",[],"id",[],"name",[]);

for i_fjord=1:n_fjords
    zf_obs = fjord_model(i_fjord).c.zf;
    tf_obs = fjord_model(i_fjord).c.tf;
    sf_obs = fjord_model(i_fjord).c.sf;
    zs_obs = flip(fjord_model(i_fjord).c.zs);
    ts_obs = interp1(zs_obs,flip(fjord_model(i_fjord).c.ts),zf_obs,'nearest','extrap');
    ss_obs = interp1(zs_obs,flip(fjord_model(i_fjord).c.ss),zf_obs,'nearest','extrap');

    n_completed = 0; % completed runs for that fjord
    if isempty(tf_obs) || isempty(sf_obs) % we disregard this fjord because we actually do not have the data
        continue
    end

    tf_box = NaN([fjord_model(i_fjord).p.N,n_runs,length(tgt_days)]);
    sf_box = NaN(size(tf_box));
    tf_box_comp = NaN([length(zf_obs),n_runs,length(tgt_days)]);
    sf_box_comp = NaN(size(tf_box_comp));
    rmse_tf = NaN([n_runs,length(tgt_days)]);
    rmse_sf = NaN([n_runs,length(tgt_days)]);
    tupper_box_comp = NaN([length(fjord_model(i_fjord).p.t_save),n_runs]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    supper_box_comp = NaN(size(tupper_box_comp));
    sinter_box_comp = NaN(size(tupper_box_comp));
    slower_box_comp = NaN(size(tupper_box_comp));

    for i_run=1:n_runs
        if ~isempty(ensemble(i_fjord,i_run).s) & ~isnan(ensemble(i_fjord,i_run).s.Tfinal(1))

            % we get a profile for each of our target days - if something
            % goes wrong here it's because of the sign/orientation of depth
            % profiles
            for i_day=1:length(tgt_days) 
                tf_box(:,i_run,i_day) = ensemble(i_fjord,i_run).s.Tfinal(:,i_day); 
                sf_box(:,i_run,i_day) = ensemble(i_fjord,i_run).s.Sfinal(:,i_day);
                
                tf_box_comp(:,i_run,i_day) = interp1(ensemble(i_fjord,i_run).s.z,tf_box(:,i_run,i_day),zf_obs,'nearest','extrap');
                sf_box_comp(:,i_run,i_day) = interp1(ensemble(i_fjord,i_run).s.z,sf_box(:,i_run,i_day),zf_obs,'nearest','extrap');
            end

            tupper_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tupper;
            tinter_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tinter;
            tlower_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Tlower;

            supper_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Supper;
            sinter_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Sinter;
            slower_box_comp(:,i_run) = ensemble(i_fjord,i_run).s.Slower;

            min_depth_rmse = 0;
            depths_rmse = (zf_obs > min_depth_rmse) & (zf_obs < ensemble(i_fjord,i_run).p.Hgl);
            for i_day = 1:length(tgt_days)
                if size(tf_obs(depths_rmse),1) ~= size(tf_box_comp(depths_rmse,i_run,i_day),1)
                    tf_obs_ref = tf_obs(depths_rmse)';
                    sf_obs_ref = sf_obs(depths_rmse)';
                else
                    tf_obs_ref = tf_obs(depths_rmse);
                    sf_obs_ref = sf_obs(depths_rmse);
                end

                rmse_tf(i_run,i_day) = rmse(tf_box_comp(depths_rmse,i_run,i_day),tf_obs_ref,'omitnan');%./mean(tf_obs(depths_rmse),'omitnan');
                rmse_sf(i_run,i_day) = rmse(sf_box_comp(depths_rmse,i_run,i_day),sf_obs_ref,'omitnan');%./mean(sf_obs(depths_rmse),'omitnan');
                % rmse_tf(i_run,i_day) = rmse(tf_box_comp(depths_rmse,i_run,i_day),tf_obs(depths_rmse),'omitnan')./std(tf_obs(depths_rmse),'omitnan');
                % rmse_sf(i_run,i_day) = rmse(sf_box_comp(depths_rmse,i_run,i_day),sf_obs(depths_rmse),'omitnan')./std(sf_obs(depths_rmse),'omitnan');
                % rmse_tf(i_run,i_day) = rmse(tf_box_comp(depths_rmse,i_run,i_day),tf_obs(depths_rmse),'omitnan')./(max(tf_obs(depths_rmse),[],'omitnan')-min(tf_obs(depths_rmse),[],'omitnan'));
                % rmse_sf(i_run,i_day) = rmse(sf_box_comp(depths_rmse,i_run,i_day),sf_obs(depths_rmse),'omitnan')./(max(sf_obs(depths_rmse),[],'omitnan')-min(sf_obs(depths_rmse),[],'omitnan'));
            end
            n_completed = n_completed+1;
            zf_box = ensemble(i_fjord,i_run).s.z;
        end
    end % i_run
    if size(ts_obs(depths_rmse),1) ~= size(tf_obs_ref,1)
        ts_obs_ref = ts_obs(depths_rmse)';
        ss_obs_ref = ss_obs(depths_rmse)';
    else
        ts_obs_ref = ts_obs(depths_rmse);
        ss_obs_ref = ss_obs(depths_rmse);
    end

    res_obs(i_fjord).zf = zf_obs;
    res_obs(i_fjord).tf = tf_obs;
    res_obs(i_fjord).sf = sf_obs;
    
    res_obs(i_fjord).zs = fjord_model(i_fjord).c.zs;
    res_obs(i_fjord).ts = fjord_model(i_fjord).c.ts;
    res_obs(i_fjord).ss = fjord_model(i_fjord).c.ss;

    res_box(i_fjord).t  = fjord_model(i_fjord).p.t_save;
    res_box(i_fjord).zf = zf_box;
    res_box(i_fjord).tf    = squeeze(mean(tf_box,2,'omitnan'));
    res_box(i_fjord).tfmin = squeeze(min(tf_box,[],2,'omitnan'));
    res_box(i_fjord).tfmax = squeeze(max(tf_box,[],2,'omitnan'));
    res_box(i_fjord).tf1sl = res_box(i_fjord).tf - squeeze(std(tf_box,[],2,'omitnan'));
    res_box(i_fjord).tf1su = res_box(i_fjord).tf + squeeze(std(tf_box,[],2,'omitnan'));

    res_box(i_fjord).sf    = squeeze(mean(sf_box,2,'omitnan'));
    res_box(i_fjord).sfmin = squeeze(min(sf_box,[],2,'omitnan'));
    res_box(i_fjord).sfmax = squeeze(max(sf_box,[],2,'omitnan'));
    res_box(i_fjord).sf1sl = res_box(i_fjord).sf - squeeze(std(sf_box,[],2,'omitnan'));
    res_box(i_fjord).sf1su = res_box(i_fjord).sf + squeeze(std(sf_box,[],2,'omitnan'));

    res_box(i_fjord).Tupper = mean(tupper_box_comp,2,'omitnan');
    res_box(i_fjord).Tupper_min = min(tupper_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tupper_max = max(tupper_box_comp,[],2,'omitnan');

    res_box(i_fjord).Tinter = mean(tinter_box_comp,2,'omitnan');
    res_box(i_fjord).Tinter_min = min(tinter_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tinter_max = max(tinter_box_comp,[],2,'omitnan');

    res_box(i_fjord).Tlower = mean(tlower_box_comp,2,'omitnan');
    res_box(i_fjord).Tlower_min = min(tlower_box_comp,[],2,'omitnan');
    res_box(i_fjord).Tlower_max = max(tlower_box_comp,[],2,'omitnan');

    res_box(i_fjord).Supper = mean(supper_box_comp,2,'omitnan');
    res_box(i_fjord).Sinter = mean(sinter_box_comp,2,'omitnan');
    res_box(i_fjord).Slower = mean(slower_box_comp,2,'omitnan');

    res_box(i_fjord).ensemble_tf = tf_box;
    res_box(i_fjord).ensemble_sf = sf_box;
    res_box(i_fjord).rmse_tf = rmse_tf;
    res_box(i_fjord).rmse_sf = rmse_sf;

    res_box(i_fjord).rmse_ts = rmse(ts_obs_ref,tf_obs_ref,'omitnan');%./mean(tf_obs_ref,'omitnan');
    res_box(i_fjord).rmse_ss = rmse(ss_obs_ref,sf_obs_ref,'omitnan');%./mean(sf_obs_ref,'omitnan');

    res_box(i_fjord).id = fjord_model(i_fjord).m.ID;
    res_box(i_fjord).name = fjord_model(i_fjord).m.name{1};
    res_box(i_fjord).n = n_completed./n_runs * 100;
end

end