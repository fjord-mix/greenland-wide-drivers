function [temp_out,salt_out] = heave_profiles(Tz,Sz,z_in,depression)
warning('off','all')
try
% default is to heave the isopycnal by 20 m ("depression" amplitude of 40)
if nargin < 4, depression = ones(size(Sz,2))*40; end

% Make the vertical profiles at regular depth intervals (1m)

% nz_orig = 1:np;
z_reg = min(z_in):1:-0.5;
temp_reg = interp1(z_in,Tz,z_reg);
salt_reg = interp1(z_in,Sz,z_reg);
np = length(z_reg); % number of points in profile

% The algorithm below is adapted from Cowton et al. (2016; Journal of Glaciology)
temp_stretch=NaN(size(temp_reg));
salt_stretch=NaN(size(salt_reg));
for i=1:size(salt_reg,2)
    % we only apply the isopycnal heaving if the change is larger than 1m because of the profile vertical resolution
    if abs(depression(i)) > 1
        % find where the pycnocline is
        sigma_profile = gsw_sigma0(salt_reg(:,i),temp_reg(:,i)); % depth in m is roughly equivalent to pressure in dbar

        [~,indx] = min(abs(sigma_profile-27.3));
        oldstrat = indx(end);

        % [~,i_z] = findpeaks(-diff(sigma_profile),'NPeaks',3);
        % pyc=sigma_profile(i_z(1));
        % % oldstrat = findnearest(pyc,sigma_profile);
        % % oldstrat = oldstrat(1);
        % oldstrat = i_z(end); 
        newstrat = abs(oldstrat+floor(depression(i))); % determine where it should end up at
        
        % stretch salinity
        var = salt_reg(:,i);
        varnew = NaN([np,1]);
        if newstrat < np
            varnew(1:newstrat) = interp1(1:oldstrat,var(1:oldstrat),oldstrat/newstrat:oldstrat/newstrat:oldstrat,'linear',var(oldstrat));
            varnew(newstrat+1:np) = interp1(oldstrat+1:np,var(oldstrat+1:np),oldstrat+1:((np-oldstrat-1)/(np-newstrat-1)):np,'linear',var(end));
        else % the isopycnal is reaching the surface
            varnew(1:np) = interp1(1:oldstrat,var(1:oldstrat),1:np,'linear',var(oldstrat));
        end
        X = ~isnan(varnew); % get rid of potential gaps
        Y = cumsum(X-diff([1;X])/2);
        varnew = interp1(1:nnz(X),varnew(X),Y,'linear',varnew(end));
        salt_stretch(:,i) = varnew;
        
        % stretch temperature
        var = temp_reg(:,i);
        varnew = NaN([np,1]);
       if newstrat < np
            varnew(1:newstrat) = interp1(1:oldstrat,var(1:oldstrat),oldstrat/newstrat:oldstrat/newstrat:oldstrat,'linear',var(oldstrat));
            varnew(newstrat+1:np) = interp1(oldstrat+1:np,var(oldstrat+1:np),oldstrat+1:((np-oldstrat-1)/(np-newstrat-1)):np,'linear',var(end));
        else % the isopycnal is reaching the surface
            varnew(1:np) = interp1(1:oldstrat,var(1:oldstrat),1:np,'linear',var(oldstrat));
        end
        X = ~isnan(varnew);
        Y = cumsum(X-diff([1;X])/2);
        varnew = interp1(1:nnz(X),varnew(X),Y,'linear',varnew(end));
        temp_stretch(:,i) = varnew;
    else
        salt_stretch(:,i) = salt_reg(:,i);
        temp_stretch(:,i) = temp_reg(:,i);
    end
end

salt_out = interp1(z_reg,salt_stretch,z_in);
temp_out = interp1(z_reg,temp_stretch,z_in);
for i=1:size(salt_out,2)
    for k=1:size(salt_out,1)
        if isnan(salt_out(k,i))
            salt_out(k,i) = salt_out(k-1,i);
            temp_out(k,i) = temp_out(k-1,i);
        end
    end
end

% % sanity check for the isopycnal heaving
% for i=1:100:size(Sz,2)
%     figure(99);
%     subplot(1,2,1)
%     plot(Sz(:,i),z_in); hold on;
%     plot(salt_stretch(:,i),z_reg); 
%     hold off;
%     xlim([33 36])
%     subplot(1,2,2)
%     plot(Tz(:,i),z_in); hold on;
%     plot(temp_stretch(:,i),z_reg);
%     hold off;
%     xlim([-0.5 5])
% end

warning('on','all')
catch ME
    fprintf('heave did not work here: %s\n',ME.message)
end
end