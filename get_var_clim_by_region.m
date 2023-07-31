function [var_clim, var_anom, depths] = get_var_clim_by_region(fjords_processed,var)
    % Get one average T,S profile per region
    depths = fjords_processed(1).f.zs;
    dims_var = size(fjords_processed(1).f.(var));
    if size(dims_var,1) > 1
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
        var_detrend = detrend(fjords_processed(i_fjord).f.(var)');
        if size(dims_var,1) > 1
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
        if size(dims_var,1) > 1
            var_series(:,:,i_reg) = var_series(:,:,i_reg)./n_reg(i_reg);
        else
            var_series(:,i_reg) = var_series(:,i_reg)./n_reg(i_reg);
        end
    end
    
    % Compute climatology in [month,depth,region] dimensions
    if size(dims_var,1) > 1
        var_reshape = reshape(var_series,12,[],length(fjords_processed(1).f.zs),7);
    else
        var_reshape = reshape(var_series,12,[],7);
    end
    var_clim    = squeeze(mean(var_reshape,2));
    
    % Compute anomalies anom = var - clim
    var_anom = NaN(size(var_reshape));

    for y=1:size(var_reshape,2)
        if length(dims_var) > 1
            var_anom(:,y,:,:) = squeeze(var_reshape(:,y,:,:)) - var_clim;
        else
            var_anom(:,y,:) = squeeze(var_reshape(:,y,:)) - var_clim;
        end
    end
    
    % get the anomalies in [time,depth,region] dimensions
    var_anom = reshape(var_anom,dims);
    
end