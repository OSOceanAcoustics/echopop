function   krig_datachk(opt)
%% function   krig_datachk(opt)
%% datachk.m perform several operations on the original data
%% 1. check raw data to remove overlapping data points	
%% 2. data subsampling using averaged data over support ninc
%% 3. coordinate conversion 
%%	opt = 1			called from data preparation window
%%		  2			called from kriging window using the same data format as in the dataprep window.
%%        3         called from kriging window using a specified data format.
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global para data hdl

DEG2RAD=pi/180;
EPS=1e-10;
    
if  para.dataprep.ext_prog== 1
	x=data.in.var1;
	y=data.in.var2*para.dataprep.ytox;
	z=data.in.var3*para.dataprep.ztox;
	var=data.in.var;
else
   switch para.dataprep.x_axis_indx			% x-axis variable
      case 1
         data.in.var1=data.in.var1_raw;
      case 2
         data.in.var1=data.in.var2_raw;
      case 3
         data.in.var1=data.in.var3_raw;
   end
   switch para.dataprep.y_axis_indx			% y-axis variable
      case 1
         data.in.var2=data.in.var1_raw;
      case 2
         data.in.var2=data.in.var2_raw;
      case 3
         data.in.var2=data.in.var3_raw;
   end
   if data.in.dim == 3   
   	switch para.dataprep.z_axis_indx			% z-axis variable
      	case 1
         	data.in.var3=data.in.var1_raw;
      	case 2
         	data.in.var3=data.in.var2_raw;
      	case 3
         	data.in.var3=data.in.var3_raw;
      end
   else
      data.in.var3=data.in.var3_raw;
   end
   x=data.in.var1;
   y=data.in.var2*para.dataprep.ytox;
   z=data.in.var3*para.dataprep.ztox;
   var=data.in.var_raw;
end


% remove data points that have the same coordinates
% while 1
%   if data.in.dim == 3
%     r=sqrt(diff(x).^2+diff(y).^2+diff(z).^2);
%   else
%     r=sqrt(diff(x).^2+diff(y).^2);
%   end
%   rmin=0.001*mean(r);
%   indxs=find(r < rmin);
%   if ~isempty(indxs)
%      x(indxs)=[];
%      y(indxs)=[];
%      var(indxs)=[];
%      if data.in.dim == 3
%        z(indxs)=[];
%      end
%    else
%      break
%   end
% end

%% remove data points with NaN's
data.in.var1=x;
data.in.var2=y/para.dataprep.ytox;
if data.in.dim == 3
    data.in.var3=z/para.dataprep.ztox;
end
data.in.var=var;

if data.in.dim == 2
	indx1=find(isnan(x)|isnan(y)|isnan(var));
   if ~isempty(indx1)
      x(indx1)=[];
      y(indx1)=[];
      var(indx1)=[];
   end
else
   indx1=find(isnan(x)|isnan(y)|isnan(z)|isnan(var));
   if ~isempty(indx1)
      x(indx1)=[];
      y(indx1)=[];
      z(indx1)=[];
      var(indx1)=[];
   end
end

data.in.x=x;
data.in.y=y;
data.in.z=z;
data.in.v=var;

% data reduction
datareduction(opt);

if para.proc.default_para ~= 1
    %%%%%%  data x-y scale normalization (0-1)%%%%%%
    if para.krig.batch_file_proc == 0 | (para.krig.batch_file_proc == 1 & para.krig.bat_proc_cnt == 1)
        para.dataprep.x_norm=max(x)-min(x);
        para.dataprep.y_norm=max(y)-min(y);
%         para.dataprep.x_offset=mean_nan(x);
%         para.dataprep.y_offset=mean_nan(y);
   %%%% to be consistent with using default kriging parameters, where
   %%%% para.dataprep.x_offset & para.dataprep.y_offset are aet in the
   %%%% default_para file (7/19/2013)
        para.dataprep.x_offset=para.krig.shelf_break_Ref_lon;
        [vario_krig_para,vario_krig_para_str]=xlsread(para.krig.vario_krig_para_filename);
        indx=find(strcmp(vario_krig_para_str,'dataprep.y_offset') == 1);
        para.dataprep.y_offset=vario_krig_para(indx);
        if para.dataprep.x_norm == 0			% x = const.
            para.dataprep.x_norm=max(EPS,x(1));
        end
        if para.dataprep.y_norm == 0			% y = const.
            para.dataprep.y_norm=max(EPS,y(1));
        end
        if data.in.dim == 3
            para.dataprep.z_norm=max(z)-min(z);
            para.dataprep.z_offset=mean_nan(z);
            if para.dataprep.z_norm == 0			% z = const.
                para.dataprep.z_norm=max(EPS,z(1));
            end
        end
    end
end

if para.krig.load_data_format_file == 0 & para.proc.default_para ~= 1
  para.dataprep.var1_indx=get(hdl.dataprep.var1,'value');
  para.dataprep.var2_indx=get(hdl.dataprep.var2,'value');
elseif para.proc.default_para == 1  % (lat, lon) rather than (lon, lat)
    para.dataprep.var1_indx=1;
    para.dataprep.var2_indx=2;
end
%% convert longitude to latitude degrees
if para.dataprep.var1_indx == 1 & para.dataprep.var2_indx == 2   % Long/Lat
   para.dataprep.latlonfac=cos(DEG2RAD*data.in.y);
elseif para.dataprep.var1_indx == 2 & para.dataprep.var2_indx == 1   % Lat/Long
   para.dataprep.latlonfac=cos(DEG2RAD*data.in.x);
else
   para.dataprep.latlonfac=1;
end

%% save the un-normalized coordinates 
data.in.x0=data.in.x;
data.in.y0=data.in.y;
data.in.z0=data.in.z;

if para.dataprep.var1_indx == 1 & para.dataprep.var2_indx == 2  %y-axis -> Lat
  data.in.x=(data.in.x-para.dataprep.x_offset).*para.dataprep.latlonfac./para.dataprep.x_norm;
  data.in.y=(data.in.y-para.dataprep.y_offset)/para.dataprep.y_norm;
elseif para.dataprep.var1_indx == 2 & para.dataprep.var2_indx == 1  		%x-axis -> Lat															  
  data.in.x=(data.in.x-para.dataprep.x_offset)/para.dataprep.x_norm;
  data.in.y=(data.in.y-para.dataprep.y_offset).*para.dataprep.latlonfac/para.dataprep.y_norm;
else														% other variables
  data.in.x=(data.in.x-para.dataprep.x_offset)/para.dataprep.x_norm;
  data.in.y=(data.in.y-para.dataprep.y_offset)/para.dataprep.y_norm;   
end
if data.in.dim == 3
  data.in.z=(data.in.z-para.dataprep.z_offset)/para.dataprep.z_norm;
end
return