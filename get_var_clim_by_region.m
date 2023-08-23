function [var_clim_forcing, var_anom, decay_func, depths] = get_var_clim_by_region(fjords_processed,var)
    % Get one average T,S profile per region
    depths = fjords_processed(1).f.zs;
    dims_var = size(fjords_processed(1).f.(var));
    if dims_var(1) > 1
        dims = [length(fjords_processed(1).t), length(depths),7];
    else
        dims = [length(fjords_processed(1).t),7];
    end
    
    var_series = zeros(dims);
    n_reg = zeros([7,1]);

    % assign ocean profiles to the fjord's respective region
    for i_fjord=1:length(fjords_processed)
        switch fjords_processed(i_fjord).m.region
            case 'SW'
                i_reg=1;
            case 'SE'
                i_reg=2;
            case 'CW'
                i_reg=3;
            case 'CE'
                i_reg=4;
            case 'NW'
                i_reg=5;
            case 'NE'
                i_reg=6;
            case 'NO'
                i_reg=7;
        end

        % we need to add back the mean because it is removed by the detrend function
        var_mean = mean(fjords_processed(i_fjord).f.(var),2);        
        var_detrend = detrend(fjords_processed(i_fjord).f.(var)') + var_mean'; 

        if dims_var(1) > 1
            var_series(:,:,i_reg) = var_series(:,:,i_reg) + var_detrend;
        else
            var_series(:,i_reg) = var_series(:,i_reg) + var_detrend;
        end
        n_reg(i_reg) = n_reg(i_reg)+1;

        % will take the longest depth profile (likely not necessary)
        if length(depths) < length(fjords_processed(i_fjord).f.zs)
            depths = fjords_processed(i_fjord).f.zs;
        end
    end
    if sum(n_reg) ~= length(fjords_processed)
        disp('Some fjords have a missing region?!')
        return
    end
    for i_reg=1:length(n_reg)
        if dims_var(1) > 1
            var_series(:,:,i_reg) = var_series(:,:,i_reg)./n_reg(i_reg);
        else
            var_series(:,i_reg) = var_series(:,i_reg)./n_reg(i_reg);
        end
    end
    
    % Compute climatology in [month,depth,region] dimensions
    if dims_var(1) > 1
        var_reshape = reshape(var_series,12,[],length(fjords_processed(1).f.zs),7);
    else
        var_reshape = reshape(var_series,12,[],7);
    end
    var_clim    = squeeze(mean(var_reshape,2));
    
    % Compute anomalies anom = var - clim
    var_anom = NaN(size(var_reshape));

    for y=1:size(var_reshape,2)
        if dims_var(1) > 1
            var_anom(:,y,:,:) = squeeze(var_reshape(:,y,:,:)) - var_clim;
        else
            var_anom(:,y,:) = squeeze(var_reshape(:,y,:)) - var_clim;
        end
    end
    
    % get the anomalies in [time,region] dimensions, over the desired depth
    % range
    var_anom = reshape(var_anom,dims);

    if dims_var(1) > 1
        % Our "decay function" is computed as the normalised standard
        % deviation of the time series. Could not use multi-variate
        % regression b/c our detrended data is not "positive definite"
        var_std = squeeze(std(var_series,1));
        std_norm = (var_std - min(var_std)) ./ (max(var_std) - min(var_std));
        % valid_range = std_norm>=0.75; % used so we know over which depth range to compute the anomaly from
        valid_range = std_norm~=0; % used so we know over which depth range to compute the anomaly from

        var_anom_mean = NaN([dims(1),dims(3)]);
        for i_reg=1:dims(3)
            var_anom_mean(:,i_reg) = squeeze(mean(var_anom(:,valid_range(:,i_reg),i_reg),2));
        end
        decay_func = NaN(dims);
        for i=1:dims_var(2)
            decay_func(i,:,:) = std_norm;
        end
        var_anom = var_anom_mean;
        var_clim_forcing = repmat(var_clim,size(var_anom,1)/size(var_clim,1),1,1);
        % var_clim_forcing = var_clim_forcing(:,:,:);
    else        
        decay_func=[];
        var_clim_forcing = repmat(var_clim,size(var_anom,1)/size(var_clim,1),1,1);
        % var_clim_forcing = var_clim_forcing(:,:);
    end
    

    %% Sanity check to see if what we did actually makes sense    
    % time_axis = linspace(1,size(var_anom,1),size(var_anom,1));
    % if dims_var(1) > 1
    %     var_check = NaN(size(var_series));
    %     for i_reg=1:size(var_clim_forcing,3)
    %         var_check(:,:,i_reg) = var_clim_forcing(:,:,i_reg)+decay_func(:,:,i_reg).*var_anom_mean(:,i_reg);
    %     end
    %     residuals = (var_check-var_series);%./mean(var_series,1)+1e-10; % add a regularisation term for when Ts->0
    %     err = squeeze(rmse(var_check,var_series));%./mean(var_series,1))+1e-10;%(max(var_series,1)-min(var_series,1)); % compute RMSE to have a metric for our error
    %     err_mean = mean(err(std_norm~=0),1);
    %     fprintf('Mean error for %s: %.4f (%.4f - %.4f)\n',var,mean(err_mean),min(err_mean),max(err_mean));
    % 
    %     %Plotting the comparison
    %     cmap = cmocean('balance');
    %     for i_reg=1:7
    %         figure; 
    %         subplot(1,3,1), imagesc(depths,time_axis,var_series(:,:,i_reg)); title(['original ',var]); colorbar; clim([min(var_series(:)),max(var_series(:))]);
    %         subplot(1,3,2), imagesc(depths,time_axis,var_check(:,:,i_reg)); title('reconstructed'); colorbar;    clim([min(var_series(:)),max(var_series(:))]);
    %         subplot(1,3,3), imagesc(depths,time_axis,residuals(:,:,i_reg)); title('Residuals'); colorbar        
    %         colormap(gca,cmap); clim(gca,[-max(max(abs(residuals(:)))),max(max(abs(residuals(:))))]); 
    %     end
    %     regs=linspace(1,7,7);%{'SW','SE','CW','CE','NW','NE'};
    %     figure; imagesc(depths,regs,err'); colorbar; colormap(gca,'hot'); title('mean RMSE')
    %     % figure; imagesc(depths,regs,std_norm'); colorbar; colormap(gca,cmap); clim(gca,[0,1]); title('Normalised STDev')
    %     % showcasing how bad the "worst reconstructions" are
    %     figure; hold on;         
    %     plot(time_axis,var_series(:,end,3),'-b'); plot(time_axis,var_check(:,end,3),'--b'); 
    %     plot(time_axis,var_series(:,end,5),'-r'); plot(time_axis,var_check(:,end,5),'--r'); 
    %     plot(time_axis,var_series(:,end,6),'-m'); plot(time_axis,var_check(:,end,6),'--m'); 
    %     legend('original CW','reconstructed CW','original NW','reconstructed NW','original NE','reconstructed NE')
    % 
    % else
    %     var_check = var_clim_forcing+var_anom;
    %     residuals = var_check-var_series;
    %     %err = rmse(var_check,var_series)./(max(var_series,1)-min(var_series,1)); % compute RMSE to have a metric for our error
    %     for i_reg=1:7
    %         figure; 
    %         subplot(1,2,1); hold on; 
    %         plot(var_series(:,i_reg)); plot(var_check(:,i_reg));
    %         legend(['original ',var],'reconstructed');
    %         subplot(1,2,2), plot(residuals(:,i_reg)); title('difference');
    %     end
    % end
    

    
end