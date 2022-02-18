function   SD_var = SD_map(CV_var, var1, data, para)
% map of standard deviation

SD_var = CV_var.*var1;     % standard deviation of biomass (metric tons)

sd_max=0.3*max(SD_var);
SD_disp=SD_var;
SD_disp(SD_disp> sd_max)=sd_max;
% SD_disp(end)=cv_max;
clr_max=length(colormap);
clr=max(0,floor(clr_max*(SD_disp-min(SD_disp))/(max(SD_disp)-min(SD_disp))))+1;
indx_nan=find(isnan(SD_disp) == 1);
SD_disp(indx_nan)=0;
clr=min(clr,clr_max);
cmap=colormap;

figure(10)  % color-coded kriged CV map
clf
m_proj('Mercator','lat',[34 60],'long',[-140 -120]);
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
hold on
m_plot(data.out.krig.lon,data.out.krig.lat,'.','markersize',1)
i0=1;
ncv=length(SD_var);
for i=1:ncv
    if isempty(indx_nan) | min(abs(indx_nan-i)) ~= 0
        m_plot(data.out.krig.lon(i),data.out.krig.lat(i),'.','color',cmap(clr(i),:),'markersize',2);
        if i == i0;hold on,end
    else
        i0=i0+1;
    end
end
xlabel('LONGITUDE','fontsize',16,'fontweight','bold')
ylabel('LATITUDE','fontsize',16,'fontweight','bold')
title(['COEFFICIENT OF VARIANCE (' para.survey_year ' Hake Survey)'],'fontsize',14,'fontweight','bold');
m_grid('linest','none','linewidth',2,'tickdir','in','xaxisloc','bottom');
set(gca,'dataaspectratio',[1 2 1]);

colormap('jet')
pos = get(gca,'Position');
stripe = 0.075; edge = 0.02;
[az,el] = view;
if all([az,el]==[0 90]), space = 0.05; else space = .1; end
set(gca,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
ax= axes('Position', rect);
image([0 1],[min(SD_disp) max(SD_disp)+0.03],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')');
set(ax,'Ydir','normal')
set(ax,'YAxisLocation','right')
set(ax,'xtick',[])
title(ax,'SD (mt)','fontsize',14,'fontweight','bold')
drawnow