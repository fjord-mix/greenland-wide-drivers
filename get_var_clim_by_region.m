function [var_clim, var_anom, decay_func, depths] = get_var_clim_by_region(fjords_processed,var)
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
        % deviation of the time series
        var_std = std(var_detrend,1);
        std_norm = (var_std - min(var_std)) / ( max(var_std) - min(var_std) );
        valid_range = std_norm~=0; % used so we know over which depth range to compute the anomaly from

        var_anom = squeeze(mean(var_anom(:,valid_range,:),2));
        decay_func = repmat(std_norm,60,1,7);    
    else        
        decay_func=[];
    end

    % Sanity check to see if what we did actually makes sense
    % var_check = repmat(var_clim(:,:,i_reg),5,1,1)+decay_func(:,:,i_reg).*var_anom(:,i_reg);
    % figure; 
    % subplot(1,3,1), imagesc(var_series(:,:,i_reg)); title('original'); colorbar
    % subplot(1,3,2), imagesc(var_check); title('reconstructed'); colorbar
    % subplot(1,3,3), imagesc(var_series(:,:,i_reg)-var_check); title('difference'); colorbar

    
end