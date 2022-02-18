function 	dispkrig3d(dir_index)
%% function dispkrig3d(dir_index) displays the kriging results based on
%% the option of dir_index 
%% dir_index = 0-3   : 3D kriging results
%%             4:  color-scale slider
%              5:  2D
%              6:  trackline
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global data para hdl


if ~isfield(para.dispkrig.trackline,'dispflag')
    para.dispkrig.trackline.dispflag=0;
end
hdls=findobj(hdl.dispkrig3d.h0,'string','Z');
if data.in.dim == 2
   if ~isempty(hdls)	% previous figure settings are for 3D
      dispkrig3dfig(1);
      dispkrig3d(dir_index);
   end
else
   if isempty(hdls)		% previous figure settings are for 2D
      dispkrig3dfig(1);
      dispkrig3d(dir_index);
   end
end

set(0, 'showhidden', 'on')
ch=get(gcf, 'children');

tools_hdl=findobj(ch,'Label','&Tools');	
rotate_hdl=findobj(tools_hdl,'Label','&Rotate 3D');	
if strcmp(get(rotate_hdl,'Checked'),'on') 
    rotation_checked_flag=1;
    [AZ,EL]=view;
    hdl.dispkrig3d.view_AZ=AZ;
    hdl.dispkrig3d.view_EL=EL;
else
    rotation_checked_flag=0;    
end
set(0, 'showhidden', 'off')

%disp(get(rotate_hdl,'Checked'))
clr_max=64;
DEG2RAD=pi/180;
if para.status.kriging == 1
 %  xgt=data.out.krig.xg*para.dataprep.x_norm./para.dataprep.latlonfac+para.dataprep.x_offset;
 %  ygt=data.out.krig.yg*para.dataprep.y_norm+para.dataprep.y_offset;
   xgt=data.out.krig.xg;
   ygt=data.out.krig.yg;
   if data.in.dim == 3
 %     zgt=data.out.krig.zg*para.dataprep.z_norm+para.dataprep.z_offset;
      zgt=data.out.krig.zg;
   end
else
   xmin=-0.5;xmax=0.5;dx=0.05;
   ymin=-0.5;ymax=0.5;dy=0.05;
   zmin=-0.5;zmax=0.5;dz=0.05;
   xgt=xmin:dx:xmax;
   ygt=ymin:dy:ymax;
   zgt=zmin:dz:zmax;
end

max_xindx=length(xgt);
max_yindx=length(ygt);
if data.in.dim == 3
	max_zindx=length(zgt);
   [X,Y,Z]=meshgrid(xgt,ygt,zgt);
   Z=Z*para.dataprep.z_norm+para.dataprep.z_offset;
	data.out.krig.Zg=Z;
else
   [X,Y]=meshgrid(xgt,ygt);
end   

if isfield(para.dataprep,'fileID')
  set(hdl.dispkrig3d.fileID,'string',para.dataprep.fileID);
end

if ~isfield(para.dataprep,'var1_indx')
    krig_message(1,'There are no data to diaplay !!!');
    return
end
%if get(hdl.dataprep.var1,'value') == 1 & get(hdl.dataprep.var2,'value') == 2  %y-axis -> Lat
if para.dataprep.var1_indx == 1 & para.dataprep.var2_indx == 2  %y-axis -> Lat
	Y=Y*para.dataprep.y_norm+para.dataprep.y_offset;			% latitude in degrees
  	X=X./cos(Y*DEG2RAD);
	X=X*para.dataprep.x_norm+para.dataprep.x_offset;		 
    %elseif get(hdl.dataprep.var1,'value') == 2 & get(hdl.dataprep.var2,'value') == 1  %x-axis -> Lat															  
elseif para.dataprep.var1_indx == 2 & para.dataprep.var2_indx == 1  %x-axis -> Lat															  
	X=X*para.dataprep.x_norm+para.dataprep.x_offset;		 	% latitude in degrees
  	Y=Y./cos(X*DEG2RAD);
    Y=Y*para.dataprep.y_norm+para.dataprep.y_offset;			
else
    X=X*para.dataprep.x_norm+para.dataprep.x_offset;		 
    Y=Y*para.dataprep.y_norm+para.dataprep.y_offset;			
end   

data.out.krig.Xg=X;
data.out.krig.Yg=Y;
Xgt=squeeze(mean(X(:,:,1)))';
Ygt=squeeze(Y(:,1,1));
data.out.krig.Xgt=Xgt;
data.out.krig.Ygt=Ygt;

if data.in.dim == 3
  Zgt=squeeze(Z(1,1,:));
  data.out.krig.Zgt=Zgt;
end 

switch hdl.dispkrig3d.var_index
case 1
   DispVar=data.out.krig.Vg;
case 2
   DispVar=data.out.krig.Eg;
case 3
   disp_str={'Show Plot is not valid for Validation.',
          '  ',
        'It is only valid for displaying Krig Map or Variance Map!'};
   krig_message(2,disp_str);
   return
end

if para.Matlab_Version == 7
  if data.in.dim == 3
    hori_size=0.5; 
    hori_x0=0.1;
    hori_clrbar=1.45;
  else      
    hori_size=0.6;
    hori_x0=0.15;
    hori_clrbar=1.15;
  end
else
  if data.in.dim == 3
    hori_size=0.6;
    hori_x0=0.15;
    hori_clrbar=1.4;
  else
    hori_size=0.6;
    hori_x0=0.15;
    hori_clrbar=1.3;
  end
end   

%% dir_index = 0-3 is for 3D display
if dir_index == 0| dir_index == 5
   %%% plot kriging map
%	delete(findobj(hdl.dispkrig3d.h0,'Tag','dispkrig3dAxes1'));
 	figure(hdl.dispkrig3d.h0)
  	delete(findobj(hdl.dispkrig3d.h0,'Tag','colorbar1'));
    delete(hdl.dispkrig3d.axes1);
	hdl.dispkrig3d.axes1 = axes('Parent',hdl.dispkrig3d.h0, ...
	'Color',[1 1 1], ...
    'Tag','dispkrig3dAxes1', ...
	'Position',[hori_x0 0.4 hori_size 0.5]);
%    hdl.dispkrig3d.axes1 = figure(7)
end

switch dir_index		
  case 0								% initial 3D plot
	xmax=max(data.out.krig.Xgt);
	xmin=min(data.out.krig.Xgt);
    slice_index=round(max_xindx/2);
    xdir_val=(xmax-xmin)*slice_index/max_xindx+xmin;
    xdir_val_norm=(xdir_val-xmin)/(xmax-xmin);
    slice(X,Y,Z,DispVar,[Xgt(slice_index)],[Ygt(1)] ,[Zgt(max_zindx)])
    set(hdl.dispkrig3d.xdir_slider,'value',xdir_val_norm);
    set(hdl.dispkrig3d.ydir_slider,'value',0);
    set(hdl.dispkrig3d.zdir_slider,'value',0);
    set(hdl.dispkrig3d.xdir_val,'string',sprintf('%6.5g',xdir_val));
    set(hdl.dispkrig3d.ydir_val,'string',sprintf('%6.5g',min(data.out.krig.Ygt)));
    set(hdl.dispkrig3d.zdir_val,'string',sprintf('%6.5g',max(data.out.krig.Zgt)));
    hdl.dispkrig3d.xindx=slice_index;
    hdl.dispkrig3d.yindx=1;
    hdl.dispkrig3d.zindx=max_zindx;
    set(gca,'zdir','reverse');
    switch hdl.dispkrig3d.shading_index
    case 1
      shading faceted
    case 2
	  shading flat
    case 3
      shading interp
    end
  case 1							% x-dirction
	xmax=max(data.out.krig.Xgt);
	xmin=min(data.out.krig.Xgt);
	slice_index=round(min(max_xindx,max(1,max_xindx*get(hdl.dispkrig3d.xdir_slider,'value'))));
    xdir_val=(xmax-xmin)*slice_index/max_xindx+xmin;
    slice(X,Y,Z,DispVar,[Xgt(slice_index)],[Ygt(hdl.dispkrig3d.yindx)] ,[Zgt(hdl.dispkrig3d.zindx)])
 	set(hdl.dispkrig3d.xdir_val,'String',sprintf('%6.5g',xdir_val));
    hdl.dispkrig3d.xindx=slice_index;
    switch hdl.dispkrig3d.shading_index
    case 1
      shading faceted
    case 2
	  shading flat
    case 3
      shading interp
    end
  case 2							% y-dirction
	ymin=min(data.out.krig.Ygt);
	ymax=max(data.out.krig.Ygt);
	slice_index=round(min(max_yindx,max(1,max_yindx*get(hdl.dispkrig3d.ydir_slider,'value'))));
    ydir_val=(ymax-ymin)*slice_index/max_yindx+ymin;
    slice(X,Y,Z,DispVar,[Xgt(hdl.dispkrig3d.xindx)],[Ygt(slice_index)] ,[Zgt(hdl.dispkrig3d.zindx)])
  	set(hdl.dispkrig3d.ydir_val,'String',sprintf('%6.5g',ydir_val));
    hdl.dispkrig3d.yindx=slice_index;
    switch hdl.dispkrig3d.shading_index
    case 1
      shading faceted
    case 2
	  shading flat
    case 3
      shading interp
    end
  case 3							% z-dirction
	zmax=max(data.out.krig.Zgt);
	zmin=min(data.out.krig.Zgt);
	slice_index=round(min(max_zindx,max(1,max_zindx*(1-get(hdl.dispkrig3d.zdir_slider,'value')))));
    zdir_val=(zmax-zmin)*slice_index/max_zindx+zmin;
    slice(X,Y,Z,DispVar,[Xgt(hdl.dispkrig3d.xindx)],[Ygt(hdl.dispkrig3d.yindx)] ,[Zgt(slice_index)])
    set(hdl.dispkrig3d.zdir_val,'String',sprintf('%6.5g',zdir_val));
    hdl.dispkrig3d.zindx=slice_index;
    switch hdl.dispkrig3d.shading_index
    case 1
      shading faceted
    case 2
	  shading flat
    case 3
      shading interp
    end
  case 4							% colorbar slider
    var1=get(hdl.dispkrig3d.cbar_slider_bot,'value');  
    var2=get(hdl.dispkrig3d.cbar_slider_top,'value');  
    max_var=max(max(max(DispVar)));
    min_var=min(min(min(DispVar)));
    step=(max_var-min_var)/100;
    disp_var1=min_var+(max_var-min_var)*var1;
    disp_var2=max_var-(max_var-min_var)*(1-var2);
    caxis([disp_var1 disp_var2]);
    switch hdl.dispkrig3d.shading_index
    case 1
      shading faceted
    case 2
	  shading flat
    case 3
      shading interp
    end
  case 5   	%%%%%%%% 2D display  %%%%%%%%%%%%%%%%%%%%%
   hd=pcolor(X,Y,DispVar);
   switch hdl.dispkrig3d.shading_index
   case 1
      shading faceted
   case 2
	  shading flat
   case 3
      shading interp
   end
   hold on
   if hdl.dispkrig3d.var_index  == 1   % kriged variable
 	  var=data.in.tv;
   else
 	  var=data.out.krig.err;
   end
   vcontour=auto_contour(var,para.dispkrig.num_of_contour,para.dispkrig.digits_of_contour);
%   cmap=colormap;
   if ~isempty(vcontour)
  		if isnan(vcontour) 
    		[c,H]=contour(X,Y,DispVar,'k');
  		else
    		[c,H]=contour(X,Y,DispVar,vcontour,'k');
  		end
  		hh=clabel(c,H);
   end
  case 6    %%%%%%%%% plot data on tracklines or values on ustomized grids
  % doing nothing
  otherwise
end

%% plot customized grid data
if para.krig.load_griddata_file == 1
   hold on
   cmap=colormap;
   ngd=length(data.out.krig.gx);
   if data.in.dim == 2                              % 2D grids
   %% plot color-coded values on the customized grids
      if hdl.dispkrig3d.var_index  == 1   % kriged variable
         data.out.krig.gv=griddata(X,Y,data.out.krig.Vg,data.out.krig.gx,data.out.krig.gy);
         var_gd=data.out.krig.gv;
      else											% kriging variance
         data.out.krig.ge=griddata(X,Y,data.out.krig.Eg,data.out.krig.gx,data.out.krig.gy);
         var_gd=data.out.krig.ge;
      end    
      if para.dispkrig.trackline.dispflag == 1
        cmap=colormap;
        markersize=para.dispkrig.markersize;
        Vmax=max(max(max(DispVar)));
        Vmin=min(min(min(DispVar)));
        clr_gd=min(max(floor(clr_max*(var_gd-Vmin)/(Vmax-Vmin)),0)+1,clr_max);
	    inc_gd=max(1,round(ngd/1600));
        ni=1;
 	    for i=1:inc_gd:ngd
   		  hdl.dispkrig3d.trackline.hxy(ni)=plot(data.out.krig.gx(i),data.out.krig.gy(i),'.','color',cmap(clr_gd(i),:),'markersize',para.dispkrig.customized_grid_data_markersize);
          ni=ni+1;
	    end   
      end
   else                                                 % 3D grids
   %% plot color-coded values on the customized grids
      if hdl.dispkrig3d.var_index  == 1   % kriged variable
         data.out.krig.gv=griddata(X,Y,Z,data.out.krig.Vg,data.out.krig.gx,data.out.krig.gy,data.out.krig.gy);
         var_gd=data.out.krig.gv;
      else											% kriging variance
         data.out.krig.ge=griddata(X,Y,Z,data.out.krig.Eg,data.out.krig.gx,data.out.krig.gy,data.out.krig.gy);
         var_gd=data.out.krig.ge;
      end    
      if para.dispkrig.trackline.dispflag == 1
        cmap=colormap;
        markersize=para.dispkrig.markersize;
        Vmax=max(max(max(DispVar)));
        Vmin=min(min(min(DispVar)));
        clr_gd=min(max(floor(clr_max*(var_gd-Vmin)/(Vmax-Vmin)),0)+1,clr_max);
	    inc_gd=max(1,round(ngd/1600));
        ni=1;
 	    for i=1:inc_gd:ngd
   		  hdl.dispkrig3d.trackline.hxy(ni)=plot3(data.out.krig.gx(i),data.out.krig.gy(i),data.out.krig.gy(i),'.','color',cmap(clr_gd(i),:),'markersize',para.dispkrig.customized_grid_data_markersize);
          ni=ni+1;
 	    end   
      end
   end
end

% Plot track line
if para.dispkrig.trackline.dispflag == 1
  if isfield(hdl.dispkrig3d.trackline,'hxy') & ~isempty(hdl.dispkrig3d.trackline.hxy)
      if ishandle(hdl.dispkrig3d.trackline.hxy)
         delete(hdl.dispkrig3d.trackline.hxy);
      end
      hdl.dispkrig3d.trackline.hxy=[];
  end
  if isfield(hdl.dispkrig3d.trackline,'htxt')& ~isempty(hdl.dispkrig3d.trackline.htxt)
      if  ~isempty(find(ishandle(hdl.dispkrig3d.trackline.htxt) ==1))
         delete(hdl.dispkrig3d.trackline.htxt);
      end
      hdl.dispkrig3d.trackline.htxt=[];
  end

    hold on
    cmap=colormap;
    markersize=para.dispkrig.markersize;
    Vmax=max(max(max(DispVar)));
    Vmin=min(min(min(DispVar)));
    if hdl.dispkrig3d.var_index  == 1   % kriged variable
 	   var=data.in.tv;
	   clr=min(max(floor(clr_max*(var-Vmin)/(Vmax-Vmin)),0)+1,clr_max);
	   TrackLineType=para.dispkrig.trackline.type_indx;
 	   TrackLineColor=para.dispkrig.trackline.line_color;
    else								   % kriging variance
 	   var=data.out.krig.err;
	   TrackLineType=2;
       TrackLineColor=1;
    end
	n=length(data.in.x0);
	inc=max(1,round(n/1600));
	track_color='ymcrgbwk';
	clr_indx=TrackLineColor+6;
    if data.in.dim == 2                         % 2D track-line
	  switch TrackLineType 
		case 1			% color-coded tracklines
            ni=1;
 			for i=1:inc:n
   			   hdl.dispkrig3d.trackline.hxy(ni)=plot(data.in.x0(i),data.in.y0(i),'.','color',cmap(clr(i),:),'markersize',markersize);
               ni=ni+1;
 			end
   	    case 2		% black/white tracklines
  			if n >= 100000
    			hdl.dispkrig3d.trackline.hxy=plot(data.in.x0,data.in.y0,['-' track_color(clr_indx)],'linewidth',2.5);
  			else
    			hdl.dispkrig3d.trackline.hxy=plot(data.in.x0,data.in.y0,['.' track_color(clr_indx)],'markersize',markersize);
  			end
		case 3		% attribute values
  			hdl.dispkrig3d.trackline.hxy=plot(data.in.x0,data.in.y0,['.' track_color(clr_indx)],'markersize',markersize);
  			hdl.dispkrig3d.trackline.htxt=text(data.in.x0,data.in.y0,num2str(var));
  			set(hdl.dispkrig3d.trackline.htxt,'fontsize',para.dispkrig.trackline.size_indx,'color',track_color(para.dispkrig.trackline.color_indx));
		case 4		% difference between observed and kriged attribute values 
            Is=griddata(X,Y,DispVar,data.in.x0,data.in.y0);
  			hdl.dispkrig3d.trackline.hxy=plot(data.in.x0,data.in.y0,['.' track_color(clr_indx)],'markersize',markersize);
  			hdl.dispkrig3d.trackline.htxt=text(data.in.x0,data.in.y0,num2str(var-Is));
  			set(hdl.dispkrig3d.trackline.htxt,'fontsize',para.dispkrig.trackline.size_indx,'color',track_color(para.dispkrig.trackline.color_indx));
	    end
    else                                         % 3D track-line
	  switch TrackLineType 
		case 1			% color-coded tracklines
            ni=1;
 			for i=1:inc:n
   			   hdl.dispkrig3d.trackline.hxy(ni)=plot3(data.in.x0(i),data.in.y0(i),data.in.z0(i),'.','color',cmap(clr(i),:),'markersize',markersize);
               ni=ni+1;
 			end
   	    case 2		% black/white tracklines
  			if n >= 100000
    			hdl.dispkrig3d.trackline.hxy=plot3(data.in.x0,data.in.y0,data.in.z0,['-' track_color(clr_indx)],'linewidth',2.5);
  			else
    			hdl.dispkrig3d.trackline.hxy=plot3(data.in.x0,data.in.y0,data.in.z0,['.' track_color(clr_indx)],'markersize',markersize);
  			end
		case 3		% attribute values
            krig_message(1,'Not an valid option !!')
		case 4		% difference between observed and kriged attribute values 
            krig_message(1,'Not an valid option !!')
      end   
   end  
else
  if isfield(hdl.dispkrig3d.trackline,'hxy') & ~isempty(hdl.dispkrig3d.trackline.hxy)
      delete(hdl.dispkrig3d.trackline.hxy);
      hdl.dispkrig3d.trackline.hxy=[];
  end
  if isfield(hdl.dispkrig3d.trackline,'htxt') & ~isempty(hdl.dispkrig3d.trackline.htxt)
      delete(hdl.dispkrig3d.trackline.htxt);
      hdl.dispkrig3d.trackline.htxt=[];
  end
end                 % end of displsy_trackline

if dir_index ~= 4 & dir_index ~= 6
	if para.dataprep.var1_indx <= 2 | para.dataprep.var2_indx <= 2	% either x or y is Lat or Long
     ntick=4;
     [xinc,xdits,yinc,ydits]=get_ninc(hdl.dispkrig3d.axes1,ntick);
	  if para.dataprep.var1_indx <= 2
       mapax(xinc,xdits,yinc,ydits,hdl.dispkrig3d.axes1,1);		% x-axis label
	  end
	  if para.dataprep.var2_indx <= 2
       mapax(xinc,xdits,yinc,ydits,hdl.dispkrig3d.axes1,2);		% y-axis label
	  end
	end
end

if get(hdl.dispkrig3d.xdir_reverse,'value') == 1
   set(hdl.dispkrig3d.axes1,'xdir','reverse');
end
if get(hdl.dispkrig3d.ydir_reverse,'value') == 1
   set(hdl.dispkrig3d.axes1,'ydir','reverse');
end
if data.in.dim == 3
 if get(hdl.dispkrig3d.zdir_reverse,'value') == 1
     set(hdl.dispkrig3d.axes1,'zdir','reverse');
 end
end

%%% color bar and color scale
if isfield(hdl.dispkrig3d,'colorbar') & ishandle(hdl.dispkrig3d.colorbar) 
    delete(hdl.dispkrig3d.colorbar);
end
hdl.dispkrig3d.colorbar=colorbar;
set(hdl.dispkrig3d.colorbar,'Position',[0.85 0.4 0.045 0.5],'Tag','colorbar1');

xlabelstr=para.dataprep.xlabel;
ylabelstr=para.dataprep.ylabel;
hx=xlabel(xlabelstr);
hy=ylabel(ylabelstr);


if data.in.dim == 3
	zlabelstr=para.dataprep.zlabel;
	hz=zlabel(zlabelstr);
	set([hx hy hz],'fontsize',12,'fontweight','bold')
	set(hx,'Rotation',16);
    set(hy,'Rotation',-27);
else
	set([hx hy],'fontsize',12,'fontweight','bold')
end

hcs=text(hori_clrbar,1.1,'Color Scale','sc');
set(hcs,'fontsize',10,'fontweight','bold');

if 0
switch hdl.dispkrig3d.shading_index
 case 1
   shading faceted
 case 2
	shading flat
 case 3
   shading interp
end
end

drawnow
hold off
if rotation_checked_flag == 1
    view(hdl.dispkrig3d.view_AZ,hdl.dispkrig3d.view_EL);
end

