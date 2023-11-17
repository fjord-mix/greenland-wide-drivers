function [heat_content,salt_content] = get_active_fjord_contents(fjord_run)

% integration from grounding line/sill (whichever is deeper) to the surface
length_run = size(fjord_run.H,2);
zbottom = max(abs(fjord_run.p.zgl),abs(fjord_run.p.silldepth)); % finds which one is deeper, sill or grounding line            
ints = cumsum(fjord_run.H,1);                   % gets depths of all layers for every time step            
h_int=NaN([fjord_run.p.N+fjord_run.p.sill,length_run]);  % will store the layer thicknesses we want to integrate over

for i_time=1:length_run
    kbottom = find(ints(:,i_time)>=abs(zbottom)-1e-6,1);      % finds index of "bottom" layer for every time step
    for k=1:kbottom
        if k==kbottom
            h_int(k,i_time) = abs(zbottom-ints(k-1,i_time));
        else
            h_int(k,i_time) = fjord_run.H(k,i_time);
        end
    end
end
% fjord_dens = (fjord_run.p.betaS*fjord_run.salt - fjord_run.p.betaT*fjord_run.temp);
% sc = fjord_run.salt           .*                   fjord_dens  .* fjord_run.p.L.*fjord_run.p.W.*h_int;
% hc = ((fjord_run.temp+273.15) .* fjord_run.p.cw .* fjord_dens) .* fjord_run.p.L.*fjord_run.p.W.*h_int;
% % gets output quantities in "per unit volume"
% heat_content = sum(hc./(fjord_run.p.L.*fjord_run.p.W.*h_int),1,'omitnan');
% salt_content = sum(sc./(fjord_run.p.L.*fjord_run.p.W.*h_int),1,'omitnan');

heat_content = sum((fjord_run.temp+273.15).*h_int,'omitnan')./sum(h_int,'omitnan');
salt_content = sum(fjord_run.salt.*h_int,'omitnan')./sum(h_int,'omitnan');
end