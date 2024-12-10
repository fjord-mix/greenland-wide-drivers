function [tf_out,sf_out,hf_out] = get_layered_fjord_properties(fjord_run,z_ssfc,z_deep)

if nargin < 2 || isempty(z_ssfc)
    z_ssfc = 25; % 25 metres is what we define to separate surface from subsurface (remember we filter out grounding lines shallower than 25 m)
end
if nargin < 3 || isempty(z_deep)
    z_deep = abs(fjord_run.p.zgl); % uses fjord grounding line as the lower integration boundary
end

length_run = size(fjord_run.H,2);
z_fjrd = fjord_run.p.H;

ints = cumsum(fjord_run.H,1); % gets depths of all layers for every time step

% will store the layer thicknesses we want to integrate over
h_sfc=NaN([fjord_run.p.N+fjord_run.p.sill,length_run]); % surface
h_int=NaN([fjord_run.p.N+fjord_run.p.sill,length_run]); % intermediate
h_dep=NaN([fjord_run.p.N+fjord_run.p.sill,length_run]); % deep
for i_time=1:length_run
    
    % finds lower-boundary index of surface layer
    kssfc = find(ints(:,i_time)>=abs(z_ssfc)-1e-6,1);      
    if kssfc==1
        h_sfc(1,i_time) = abs(z_ssfc); % if z=25m is in the first box, we take 25m as 'h_sfc'
    else
        for k=1:kssfc
            if k==kssfc
                h_sfc(k,i_time) = abs(z_ssfc-ints(k-1,i_time));
            else
                h_sfc(k,i_time) = fjord_run.H(k,i_time);
            end
        end
    end
    % finds lower-boundary index of intermediate layer
    kint = find(ints(:,i_time)>=abs(z_deep)-1e-6,1);      
    if kint==kssfc
        % if GL is in the same box as the 25m line, we take 'h_int' as Zgl-25m.
        % This should also take care of the case when they are both in box 1.
        h_int(kint,i_time) = abs(z_deep-z_ssfc);
    else
        for k=kssfc:kint
            if k==kssfc
                h_int(k,i_time) = abs(ints(k,i_time)-z_ssfc);
            elseif k==kint
                h_int(k,i_time) = abs(z_deep-ints(k-1,i_time));
            else
                h_int(k,i_time) = fjord_run.H(k,i_time);
            end
        end
    end
    % lower-boundary index of bottom layer is the fjord depth
    kdeep = fjord_run.p.N+fjord_run.p.sill;
    if kdeep==kint
        % if GL is in the bottom box, we take 'h_dep' as H-Zgl.
        h_dep(kdeep,i_time) = abs(z_fjrd-z_deep);
    else
        for k=kint:kdeep
            if k==kint
                h_dep(k,i_time) = abs(ints(k,i_time)-z_deep);
            else
                h_dep(k,i_time) = fjord_run.H(k,i_time);
            end
        end
    end
end

% weighted averages of T and S by how much of each boxmodel layer falls
% within the desired intervals
tf_sfc = sum((fjord_run.temp).*h_sfc,'omitnan')./sum(h_sfc,'omitnan');
tf_int = sum((fjord_run.temp).*h_int,'omitnan')./sum(h_int,'omitnan');
tf_dep = sum((fjord_run.temp).*h_dep,'omitnan')./sum(h_dep,'omitnan');

sf_sfc = sum((fjord_run.salt).*h_sfc,'omitnan')./sum(h_sfc,'omitnan');
sf_int = sum((fjord_run.salt).*h_int,'omitnan')./sum(h_int,'omitnan');
sf_dep = sum((fjord_run.salt).*h_dep,'omitnan')./sum(h_dep,'omitnan');

tf_out = [tf_sfc; tf_int; tf_dep];
sf_out = [sf_sfc; sf_int; sf_dep];

hf_out(1,:,:) = h_sfc;
hf_out(2,:,:) = h_int;
hf_out(3,:,:) = h_dep;

%% Plot to (sanity) check relative contribution of each boxmodel layer to the surf-intermediate-deep structure
% figure; hold on;
% taxis=1:3256;
% for i=1:3
%     subplot(3,1,i)
%     area(taxis,-squeeze(hf_out(i,:,:))')
% end
end