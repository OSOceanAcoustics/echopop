function bubble_plot_NASC_number_biomass(lat,lon,var,var_name,fac)
%% bubble plots of NASC, Abundance (number), or Biomass
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

if isstruct(var)
    biomass=var.biomass;
    var=var.func;
end
clr_max=64;
min_markersize=3;
markersize_factor=0.2;
figure(12)
clf
m_proj('Mercator','lat',[34 60],'long',[-136 -120]);
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
cmap=colormap('jet'); colormap(cmap);
hold on
m_plot(lon,lat,'.','markersize',min_markersize)
xlabel('Longitude (deg)','fontsize',20,'fontweight','bold')
ylabel('Latitude (deg)','fontsize',20,'fontweight','bold')
hold on
ind=find(var == 0);
var(ind)=[];
lon(ind)=[];
lat(ind)=[];
n=length(var);
clr=max(0,floor(clr_max*(var-min(var))/(fac*max(var)-min(var))))+1;
clr=min(clr,clr_max);
biomass_sz=clr_max*(biomass-min(biomass))/(max(biomass)-min(biomass));
markersize=markersize_factor*biomass_sz+min_markersize;
cmap=colormap;
for i=1:n
   m_plot(lon(i),lat(i),'o','color',cmap(clr(i),:),'markersize',markersize(i));   
end
m_grid('linest','none','linewidth',2,'tickdir','in','xaxisloc','bottom');
set(gca,'dataaspectratio',[1 2 1]);
grid
ntick=4;
hold off
set(gca,'fontsize',14)

% %% colorbar
pos = get(gca,'Position'); 
stripe = 0.075; edge = 0.02; 
[az,el] = view;
if all([az,el]==[0 90]), space = 0.05; else space = .1; end
set(gca,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
ax = axes('Position', rect);
image([0 1],[min(var) max(var)],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')'); 
set(ax,'Ydir','normal')
set(ax,'YAxisLocation','right')
set(ax,'xtick',[])
title(var_name,'fontsize',12,'fontweight','bold')


