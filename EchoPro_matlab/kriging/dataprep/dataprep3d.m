%% Load data file and specify necessary specification
% function		dataprep3d(status)
%					status = 1		load input file
%						   = 2		do not load input file
%						   = 3		do not load input file, just perform date transformation
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

function		dataprep3d(status)

global data hdl para

clr_max=length(colormap);

if status > 1 & isempty(data)
    krig_message(1,'Data have not been loaded yet !!!');
    return
end
if status == 1
    %%%%%%%%%%%%%%%%%%%%%%%% load data into memory %%%%%%%%%%%%%%%%%%
    loaddatfile(1);				% load data file from data preparation window
end

if status <= 2
    %%%%%%%%  delete nans, filtering and data reduction %%%%%%%%%
    krig_datachk(1);					% call from dataprep window
end

%% Data Transformation
data.in.tvar=datatransform(1,data.in.var);		% forward transformation on the original data
data.in.tv=datatransform(1,data.in.v);			% forward transformation on the normalized data

%%%%%%%%%%%%%%%%%  plot the data  %%%%%%%%%%%%%%%%%
if para.dataprep.data_disptype == 1
    if para.dataprep.data_disptype0 == 1
        hv1=hdl.dataprep.axes1;
    elseif para.dataprep.data_disptype0 == 2
        hv1=hdl.dataprep.axes2;
    end
    delete(hdl.dataprep.axes1_clrbar);
    para.dataprep.data_disptype0=1;   % historical setting
    hld2i=get(hv1,'children');
    for i=1:length(hld2i)
        delete(hld2i(i))
    end
    set(hv1,'xticklabel','','xtick',[])
    set(hv1,'yticklabel','','ytick',[])
    set(hv1,'yticklabel','','xticklabel','');
    set(get(hv1,'xlabel'),'string','');
    set(get(hv1,'ylabel'),'string','');
    if data.in.dim == 3
        set(hv1,'zticklabel','');
        set(get(hv1,'zlabel'),'string','');
    end
    delete(hv1);
    
    hdl.dataprep.axes1 = axes('Parent',hdl.dataprep.h0,	'Color',[1 1 1],'Position',[0.1 0.55 0.5 0.4]);
    hv1=hdl.dataprep.axes1;
end

n=length(data.in.x);
x=para.dataprep.x_norm*data.in.x./para.dataprep.latlonfac+para.dataprep.x_offset;
y=para.dataprep.y_norm*data.in.y+para.dataprep.y_offset;
var=data.in.tv;
sta=1:n;
clr=max(0,floor(clr_max*(var-min(var))/(max(var)-min(var))))+1;
clr=min(clr,clr_max);
cmap=colormap;
inc=max(1,round(n/1600));
inc=1;

%hdl.dataprep.axes1=axes('parent',hdl.dataprep.h0,'position',Pos1,'Tag','dataprepAxes1','NextPlot','replace');
%%%%%%%  Display Type 1

if get(hdl.dataprep.data_type1,'value') == 1  % 2D/3D View of Color-Coded Data
    if ~isempty(data.in.z)
        z=para.dataprep.z_norm*data.in.z+para.dataprep.z_offset;
        for i=1:inc:n
            plot3(x(i),y(i),z(i),'.','color',cmap(clr(i),:),'markersize',para.dispkrig.markersize);
            if i == 1;hold on,end
        end
        if get(hdl.dataprep.zdir,'value')
            set(hdl.dataprep.axes1,'zdir','reverse');
        end
    else
        for i=1:inc:n
            plot(x(i),y(i),'.','color',cmap(clr(i),:),'markersize',para.dispkrig.markersize);
            if i == 1;hold on,end
        end
    end
    
    if get(hdl.dataprep.xdir,'value')
        set(hdl.dataprep.axes1,'xdir','reverse');
    end
    if get(hdl.dataprep.ydir,'value')
        set(hdl.dataprep.axes1,'ydir','reverse');
    end
    hold off
    xstr1=get(hdl.dataprep.xlabel,'string');
    xstr2=get(hdl.dataprep.x_unit,'string');
    indx_str=get(hdl.dataprep.x_unit,'value');
    xlabel_str=[char(xstr1) ' ' char(xstr2(indx_str))];
    ystr1=get(hdl.dataprep.ylabel,'string');
    ystr2=get(hdl.dataprep.y_unit,'string');
    indx_str=get(hdl.dataprep.y_unit,'value');
    ylabel_str=[char(ystr1) ' ' char(ystr2(indx_str))];
    hx1=xlabel(xlabel_str);
    hy1=ylabel(ylabel_str);
    para.dataprep.xlabel=xlabel_str;
    para.dataprep.ylabel=ylabel_str;
    
    hvl1=[hx1 hy1];
    if ~isempty(data.in.z)
        zstr1=get(hdl.dataprep.zlabel,'string');
        zstr2=get(hdl.dataprep.z_unit,'string');
        indx_str=get(hdl.dataprep.z_unit,'value');
        zlabel_str=[char(zstr1) ' ' char(zstr2(indx_str))];
        hz1=zlabel(zlabel_str);
        hvl1=[hvl1 hz1];
        set(hvl1,'fontweight','bold')
        set(get(hdl.dataprep.axes1,'xlabel'),'rotation',15)
        set(get(hdl.dataprep.axes1,'ylabel'),'rotation',-22)
        para.dataprep.zlabel=zlabel_str;
    else
        set(hvl1,'fontweight','bold')
    end
    
    if get(hdl.dataprep.var1,'Value') < 3 | get(hdl.dataprep.var2,'Value') < 3	% either x or y is Lat/Long
        ntick=4;
        [xinc,xdits,yinc,ydits]=get_ninc(hdl.dataprep.axes1,ntick);
        if get(hdl.dataprep.var1,'Value') < 3 & max(data.in.var1) <= 360 & min(data.in.var1) >= -360
            mapax(xinc,xdits,yinc,ydits,hdl.dataprep.axes1,1);
        end
        if get(hdl.dataprep.var2,'Value') < 3  & max(data.in.var2) <= 360 & min(data.in.var2) >= -360
            mapax(xinc,xdits,yinc,ydits,hdl.dataprep.axes1,2);
        end
        tick_len=size(get(hdl.dataprep.axes1,'yticklabel'),2);
        switch tick_len
            case 9
                set(hdl.dataprep.axes1,'fontsize',9)
            case 10
                set(hdl.dataprep.axes1,'fontsize',8)
        end
    end
    %  drawnow
    %% colorbar
    % hdl.dataprep.axes1_clrbar=colorbar;
    pos = get(gca,'Position');
    stripe = 0.075; edge = 0.02;
    [az,el] = view;
    if all([az,el]==[0 90]), space = 0.05; else space = .1; end
    set(hdl.dataprep.axes1,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
    rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
    hdl.dataprep.axes1_clrbar = axes('Position', rect);
    ax=hdl.dataprep.axes1_clrbar;
    image([0 1],[min(var) max(var)],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')');
    set(ax,'Ydir','normal')
    set(ax,'YAxisLocation','right')
    set(ax,'xtick',[])
    
elseif  get(hdl.dataprep.data_type2,'value') == 1  % Sample sequence data
    
    if para.dataprep.data_disptype0 == 1
        hv2=hdl.dataprep.axes1;
    elseif para.dataprep.data_disptype0 == 2
        hv2=hdl.dataprep.axes2;
    end
    para.dataprep.data_disptype0=2;   % historical setting
    Pos2=get(hv2,'position');
    hld2i=get(hv2,'children');
    for i=1:length(hld2i)
        delete(hld2i(i))
    end
    
    set(hv2,'xticklabel','','xtick',[])
    set(hv2,'yticklabel','','ytick',[])
    set(hv2,'yticklabel','','xticklabel','')
    set(get(hv2,'xlabel'),'string','')
    set(get(hv2,'ylabel'),'string','')
    if data.in.dim == 3
        set(hv2,'zticklabel','');
        set(get(hv2,'zlabel'),'string','');
    end
    
    % hdl.dataprep.axes2=axes('parent',hdl.dataprep.h0,'position',Pos2,'Tag','dataprepAxes2','NextPlot','replace');
    hdl.dataprep.axes2 = axes('Parent',hdl.dataprep.h0,	'Color',[1 1 1],'Position',[0.1 0.55 0.5 0.4]);
    plot(sta,var,'-')
    hold on
    for i=1:n
        plot(sta(i),var(i),'.','color',cmap(clr(i),:),'markersize',para.dispkrig.markersize);
    end
    hold off
    
    hx2=xlabel('Sample Number');
    hy2=ylabel('Observation Value');
    hvl2=[hx2 hy2];
    set(hvl2,'fontweight','bold')
    
    delete(hdl.dataprep.axes1_clrbar);
    pos=get(gca,'Position');
    stripe = 0.075; edge = 0.02;
    [az,el] = view;
    if all([az,el]==[0 90]), space = 0.05; else space = .1; end
    set(hdl.dataprep.axes2,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
    rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
    hdl.dataprep.axes1_clrbar = axes('Position', rect);
    ax=hdl.dataprep.axes1_clrbar;
    image([0 1],[min(var) max(var)],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')');
    set(ax,'Ydir','normal')
    set(ax,'YAxisLocation','right')
    set(ax,'xtick',[])
    
end         % Data Display Type

action=check_unitsfig(1);
if get(hdl.dataprep.data_format,'value')  % save data format file
    save_data_format_info(2)
end

if ~isempty(var)
    para.status.dataprep=1;						% data have been loaded
end
variogram3dfig;