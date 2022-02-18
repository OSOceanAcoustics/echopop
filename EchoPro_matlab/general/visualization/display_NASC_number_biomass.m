function display_NASC_number_biomass(lat,lon,var,var_name,fac)
% plot NASC/Biomass/Abundance as a function of lat & lon
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

clr_max=64;
figure(3)
clf
m_proj('Mercator','lat',[34 60],'long',[-136 -120]);
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
hold on

m_plot(lon(1:end),lat(1:end),'.','markersize',3)
xlabel('Longitude (deg)','fontsize',20,'fontweight','bold')
ylabel('Latitude (deg)','fontsize',20,'fontweight','bold')
hold on
ind=find(var == 0 | isnan(var) == 1);
var(ind)=[];
lon(ind)=[];
lat(ind)=[];
n=length(var);
[cnt,var_bin]=hist(var,1000);
cnt_acc=cumsum(cnt);
cnt_ratio=cnt_acc/max(cnt_acc);
low_threshold=0.025;
high_threshold=0.975;
[v, ind_low]=min(abs(cnt_ratio-low_threshold));
[v, ind_high]=min(abs(cnt_ratio-high_threshold));

vmin=var_bin(ind_low);
vmax=var_bin(ind_high);
var(var>vmax)=vmax;
var(var<vmin)=nan;
% figure(2);plot(var_bin, cnt_ratio, [vmin vmin], [0 1],'r', [vmax vmax], [0 1],'r');
figure(3)
clr=max(0,floor(clr_max*(abs(var)-min(abs(var)))/(fac*max(abs(var))-min(abs(var)))))+1;
indx_nan=find(isnan(abs(var)) == 1 );
clr(isnan(var) == 1)=nan;
clr=min(clr,clr_max);
cmap=colormap('jet');

for i=n:-1:1
  m_plot(lon(i),lat(i),'o','color',cmap(clr(i),:),'markersize',2);   
end
hold off
set(gca,'fontsize',12)
m_grid('linest','none','linewidth',2,'tickdir','in','xaxisloc','bottom');
set(gca,'dataaspectratio',[1 2 1]);

%% colorbar
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


