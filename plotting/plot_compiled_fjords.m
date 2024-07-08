function plot_compiled_fjords(fjords_map,fjords_compilation,i_fjord,path_figs)

fjord = fjords_compilation(i_fjord);

%% Plot settings
cmap_bed = cptcmap('GMT_globe');
climits=[-1.5e3 1.5e3];
pad = 20e3;
lspace = 0.01;
rspace = 0.01;
tspace = 0.01;
bspace = 0.3;
pw = 1-lspace-rspace;
ph = 1-bspace-tspace;
x = fjords_map.x;
y = fjords_map.y;
[X,Y] = meshgrid(x,y);


% plot dimensions
minx = min(X(fjord.inds));
maxx = max(X(fjord.inds));
miny = min(Y(fjord.inds));
maxy = max(Y(fjord.inds));
%TODO: ocean "coupling" check if cast is inside these bounds
if isfield(fjord,"ocean")
    xp = fjord.ocean.x;
    yp = fjord.ocean.y;
    
    if xp>minx && xp<maxx && yp>miny && yp<maxy
        disp('Cast is inside plot');
    else
        disp('Cast is outside plot: fixing...');
        if xp<minx
            minx = xp;
            if yp<miny, miny = yp;
            elseif yp>maxy, maxy = yp;
            end
        elseif xp>maxx
            maxx = xp;
            if yp<miny, miny = yp;
            elseif yp>maxy, maxy = yp;
            end
        elseif yp<miny
            miny = yp;
        elseif yp>maxy
            maxy = yp;
        end
    end
end
dx = maxx-minx;
dy = maxy-miny;
[dmax,coord] = max([dx,dy]);
if coord==1
    xlims=[minx-pad,maxx+pad];
    ylims=0.5*(miny+maxy)+(dmax/2+pad)*[-1,1];
else
    xlims=0.5*(minx+maxx)+(dmax/2+pad)*[-1,1];
    ylims=[miny-pad,maxy+pad];
end
fb2 = fjords_map.fb;
fb2(fjord.inds)=2;

%% Plotting
% if nargin > 3, figure('Visible','off'); else, figure; end
figure;

i_entry=1;
a1 = axes('position',[lspace,bspace,pw,ph]); hold on;
% a1 = axes; hold on;
hp1 = pcolor(x,y,fjords_map.bed); axis xy;
set(hp1,'EdgeColor','none')
colormap(a1,cmap_bed); clim(climits);
xlim(xlims); ylim(ylims);
set(a1,'box','on','xtick',[],'ytick',[]);
a2 = axes('position',get(a1,'position')); hold on;
hp2 = pcolor(x,y,fjords_map.ice_mask); axis xy;
set(hp2,'EdgeColor','none')
colormap(a2,0.75*[1,1,1]);
% q(i_entry)=plot(xb,yb,'k--','linewidth',1); % we do not have the ice-ocean boundary coords at the moment
% leg_entry{i_entry} = 'fjords-shelf boundary';
% i_entry=i_entry+1; 

q(i_entry)=plot(X(fjord.b_inds),Y(fjord.b_inds),'k^','markerfacecolor','r','markersize',5);
leg_entry{i_entry}='fjord mouth';
i_entry=i_entry+1;

contour(x,y,fb2,[2,2],'r','linewidth',1);
q(i_entry)=plot(NaN,NaN,'r','linewidth',1);
leg_entry{i_entry}='fjord extent';
i_entry=i_entry+1; 

for i_glacier=1:length(fjord.glaciers)
    q(i_entry)=plot([fjord.glaciers.x],[fjord.glaciers.y],'kp','markerfacecolor','y','markersize',20);
    % q(i_entry)=plot(fjord.glaciers(i_glacier).x,fjord.glaciers(i_glacier).y,'kp','markerfacecolor','r','markersize',13);
end
leg_entry{i_entry}='glaciers';
i_entry=i_entry+1;

if isfield(fjord,"ocean")
    q(i_entry)=plot(fjord.ocean.x,fjord.ocean.y,'ko','markerfacecolor','r','markersize',15);
    leg_entry{i_entry}='Ocean profile';
    % i_entry=i_entry+1;
end

set(a2,'box','off','xtick',[],'ytick',[],'visible','off','color','none');
xlim(xlims); ylim(ylims);
legend(q,leg_entry);
% colorbar(a1,'Westoutside')
% lp = l.Position;
% l.Position = [0.75,bspace-0.13,0,0];

%% Fjord stats
% dx = 0.02;
% dy = 0.02;
% d1 = 0.01;

% annotation('textbox',[lspace+dx,bspace-0*dy-1*d1,1,0],'string',['Glacier: ',twglaciers(i).name,' (',num2str(0.01*round(100*twglaciers(i).lat)),char(176),'N, ',num2str(0.01*round(100*abs(twglaciers(i).lon))),char(176),'W)'],'fontsize',12,'interpreter','tex','edgecolor','none');
% 
% annotation('textbox',[lspace+dx,bspace-1*dy-3*d1,1,0],'string',['Area = ',num2str(round(twglaciers(i).fjord.area)),' km$^2$'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-2*dy-3*d1,1,0],'string',['Volume = ',num2str(round(twglaciers(i).fjord.vol)),' km$^3$'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-3*dy-3*d1,1,0],'string',['Length = ',num2str(round(twglaciers(i).fjord.length)),' km'],'fontsize',12,'interpreter','latex','edgecolor','none');
% 
% annotation('textbox',[lspace+dx,bspace-4*dy-4*d1,1,0],'string',['Mean depth = ',num2str(round(twglaciers(i).fjord.meandepth)),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-5*dy-4*d1,1,0],'string',['Max depth = ',num2str(round(twglaciers(i).fjord.maxdepth)),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-6*dy-4*d1,1,0],'string',['GL depth = ',num2str(round(abs(twglaciers(i).gldepth))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-7*dy-4*d1,1,0],'string',['Effective depth = ',num2str(round(abs(twglaciers(i).fjord.effdepth))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-8*dy-4*d1,1,0],'string',['Plume NB depth = ',num2str(round(abs(twglaciers(i).plume(1).znb))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
% 
% annotation('textbox',[lspace+dx,bspace-9*dy-5*d1,1,0],'string',['Mean SGD flux = ',num2str(round(twglaciers(i).plume(1).Qgl)),' m$^3$/s = ',num2str(0.01*round(100*twglaciers(i).plume(1).Qgl/1000)),' mSv'],'fontsize',12,'interpreter','latex','edgecolor','none');
% annotation('textbox',[lspace+dx,bspace-10*dy-5*d1,1,0],'string',['Upwelled flux = ',num2str(round(twglaciers(i).plume(1).Qnb)),' m$^3$/s = ',num2str(0.01*round(100*twglaciers(i).plume(1).Qnb/1000)),' mSv'],'fontsize',12,'interpreter','latex','edgecolor','none');

%% Saving
if nargin > 3
    savenum = sprintf('%04d',i_fjord);
    fw = 14;
    fh = fw*pw/ph;
    % exportgraphics(gcf,[path_figs,'fjord_bnds',savenum,'.png'],'Resolution',300);
    saveplot(fw,fh,300,[path_figs,'fjord_bnds',savenum,'.png']);
    close all;
end

end