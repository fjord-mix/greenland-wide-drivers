function [tocn_clim, socn_clim, tocn_anom, socn_anom, depths] = get_ocean_clim_by_region(fjords_processed)
    % Get one average T,S profile per region
    depths = fjords_processed(1).f.zs;
    dims_ocn = [length(fjords_processed(1).t), length(depths),7];
    
    tocn_series = zeros(dims_ocn);
    socn_series = zeros(dims_ocn);
    n_ocn = zeros([7,1]);

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
        tocn_detrend = detrend(fjords_processed(i_fjord).f.Ts');
        socn_detrend = detrend(fjords_processed(i_fjord).f.Ss');
        tocn_series(:,:,i_reg) = tocn_series(:,:,i_reg) + tocn_detrend;
        socn_series(:,:,i_reg) = socn_series(:,:,i_reg) + socn_detrend;
        n_ocn(i_reg) = n_ocn(i_reg)+1;

        % will take the longest depth profile (likely not necessary)
        if length(depths) < length(fjords_processed(i_fjord).f.zs)
            depths = fjords_processed(i_fjord).f.zs;
        end
    end
    if sum(n_ocn) ~= length(fjords_processed)
        disp('Some fjords have a missing region?!')
        return
    end
    for i_reg=1:length(n_ocn)
        tocn_series(:,:,i_reg) = tocn_series(:,:,i_reg)./n_ocn(i_reg);
        socn_series(:,:,i_reg) = socn_series(:,:,i_reg)./n_ocn(i_reg);
    end
    
    % Compute T and S climatologies (Tclim, Sclim) in [month,depth,region] dimensions
    tocn_reshape = reshape(tocn_series,12,[],length(fjords_processed(1).f.zs),7);
    socn_reshape = reshape(socn_series,12,[],length(fjords_processed(1).f.zs),7);
    tocn_clim    = squeeze(mean(tocn_reshape,2));
    socn_clim    = squeeze(mean(socn_reshape,2));
    
    % Compute anomalies (Ta,Sa) = (T,S) - (Tclim,Sclim)
    tocn_anom = NaN(size(tocn_reshape));
    socn_anom = NaN(size(socn_reshape));
    for y=1:size(tocn_reshape,2)
        tocn_anom(:,y,:,:) = squeeze(tocn_reshape(:,y,:,:)) - tocn_clim;
        socn_anom(:,y,:,:) = squeeze(socn_reshape(:,y,:,:)) - socn_clim;    
    end
    
    % get the anomalies in [time,depth,region] dimensions
    tocn_anom = reshape(tocn_anom,dims_ocn);
    socn_anom = reshape(socn_anom,dims_ocn);
    
end