function loaddatfile(opt,filename)
%%	opt = 1			called from data preparation window
%%		  2			called from kriging window
%%		  3			called from visualization window
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para data hdl

if opt == 1								% called from data preparation window
 index1=get(hdl.dataprep.x_axis,'value');
 index2=get(hdl.dataprep.y_axis,'value');
 index3=get(hdl.dataprep.z_axis,'value');
 if get(hdl.dataprep.ext_prog,'Value') == 1  & ~isempty(get(hdl.dataprep.dat_conv_fname,'string'))
   eval(['addpath ''' para.file_dir.data_conversion ''''])
   fname_len=length(para.dataprep.dat_conv_fname);
 %  cmd=['[x_axis,y_axis,z_axis,val,data_index]=' para.dataprep.dat_conv_fname '(''' para.dataprep.filename ''');'];
   data_index=[1 2 3 ];
   cmd=['[x_axis,y_axis,z_axis,val]=' para.dataprep.dat_conv_fname '(''' para.dataprep.filename ''');'];
   eval(cmd);
   data.in.var1_raw=x_axis;
   data.in.var2_raw=y_axis;
   set(hdl.dataprep.x_axis,'value',data_index(1));
   set(hdl.dataprep.y_axis,'value',data_index(2));
   if ~isempty(z_axis) & ~isempty(val)
     data.in.var3_raw=z_axis;
     data.in.var_raw=val;
 	 data.in.dim=3;
 %    set(hdl.dataprep.z_axis,'enable','on');  %
     set(hdl.dataprep.var3,'enable','on');  %
     set(hdl.dataprep.zlabel,'enable','on');  %
     set(hdl.dataprep.z_unit,'enable','on');  %
     set(hdl.dataprep.y_axis,'value',data_index(3));
  else
     if isempty(z_axis)
       data.in.var3_raw=[];
       data.in.var_raw=val;
     elseif isempty(val)
       data.in.var3_raw=[];
       data.in.var_raw=z_axis;
     end   
 	 data.in.dim=2;
     set(hdl.dataprep.z_axis,'enable','off');  %
     set(hdl.dataprep.var3,'enable','off');  %
     set(hdl.dataprep.zlabel,'enable','off');  %
     set(hdl.dataprep.z_unit,'enable','off');  %
   end
   data.in.var1=data.in.var1_raw;
   data.in.var2=data.in.var2_raw;
   data.in.var3=data.in.var3_raw;
   data.in.var=data.in.var_raw;
 else									% load data directly
   if para.proc.kriging_auto ~= 1
       eval(['dat=load(' '''' para.dataprep.filename ''');'])
   else
       dat=data.in.US_CAN_dat;
   end
   data.in.var1_raw=dat(:,1);
   data.in.var2_raw=dat(:,2);
   if size(dat,2) >= 4			% 3-D data
     data.in.var3_raw=dat(:,3);
     data.in.var_raw=dat(:,4);
 	 data.in.dim=3;
     set(hdl.dataprep.z_axis,'enable','on');  %
     set(hdl.dataprep.var3,'enable','on');  %
     set(hdl.dataprep.zlabel,'enable','on');  %
     set(hdl.dataprep.z_unit,'enable','on');  %
   else								% 2-D data
     data.in.var3_raw=[];
  	  data.in.dim=2;
     data.in.var_raw=dat(:,3);
     set(hdl.dataprep.z_axis,'enable','off');  %
     set(hdl.dataprep.var3,'enable','off');  %
     set(hdl.dataprep.zlabel,'enable','off');  %
     set(hdl.dataprep.z_unit,'enable','off');  %
   end
   data.in.var1=dat(:,index1);
   data.in.var2=dat(:,index2);
   data.in.var3=dat(:,index3);
   data.in.var=data.in.var_raw;
 end
 para.dataprep.x_axis_indx=get(hdl.dataprep.x_axis,'value');
 para.dataprep.y_axis_indx=get(hdl.dataprep.y_axis,'value');
 if data.in.dim == 3
   para.dataprep.z_axis_indx=get(hdl.dataprep.z_axis,'value');
 end
elseif opt == 2										% called from kriging window
% if ~isfield(para.dataprep,'ext_prog') para.dataprep.ext_prog=0;end
 if para.dataprep.ext_prog == 1							% external file
   eval(['addpath ''' para.file_dir.data_conversion ''''])
   cmd=['[x_axis,y_axis,z_axis,val]=' para.dataprep.dat_conv_fname '(''' para.krig.data_file ''');'];
   eval(cmd);
   data.in.var1_raw=x_axis;
   data.in.var2_raw=y_axis;
   if ~isempty(z_axis) & ~isempty(val)          % 3-D data
     data.in.var3_raw=z_axis;
     data.in.var_raw=val;
 	 data.in.dim=3;
   else                                         % 2-D data
     if isempty(z_axis)
       data.in.var3_raw=[];
       data.in.var_raw=val;
     elseif isempty(val)
       data.in.var3_raw=[];
       data.in.var_raw=z_axis;
     end   
 	  data.in.dim=2;
  	end
   data.in.var1=data.in.var1_raw;
   data.in.var2=data.in.var2_raw;
   data.in.var3=data.in.var3_raw;
   data.in.var=data.in.var_raw;
else																			% no external file
   index1=para.dataprep.x_axis_indx;
   index2=para.dataprep.y_axis_indx;
   eval(['dat=load(' '''' para.krig.data_file ''');'])
   data.in.var1_raw=dat(:,1);
   data.in.var2_raw=dat(:,2);
   if size(dat,2) >= 4			% 3-D data
     data.in.var3_raw=dat(:,3);
     data.in.var_raw=dat(:,4);
  	 index3=para.dataprep.z_axis_indx;
	 data.in.dim=3;
   else								% 2-D data
     data.in.var3_raw=[];
  	 data.in.dim=2;
     data.in.var_raw=dat(:,3);
     index3=[];
   end
   data.in.var1=dat(:,index1);
   data.in.var2=dat(:,index2);
   data.in.var3=dat(:,index3);
   data.in.var=data.in.var_raw; 
 end
elseif opt == 3										% called from visualization window
    mat_file_in=para.file_dir.mat_file_in;
	cmd=['load ''' filename ''''];
    current_matlab_version=para.Matlab_Version;
    HDIR=para.home_dir;
    eval(cmd)
    para.Matlab_Version=current_matlab_version;
    para.home_dir=HDIR;
    para.file_dir.mat_file_in=mat_file_in;
    set(hdl.dispkrig3d.fileID,'string',para.dataprep.fileID);
end
para.dataprep.checkunit_action=0;
para.krig.xmin0=min(data.in.var1);
para.krig.xmax0=max(data.in.var1);
para.krig.dx0=(max(data.in.var1)-min(data.in.var1))/19;
para.krig.ymin0=min(data.in.var2);
para.krig.ymax0=max(data.in.var2);
para.krig.dy0=(max(data.in.var2)-min(data.in.var2))/19;
if data.in.dim == 3   
     para.krig.zmin0=min(data.in.var3);
     para.krig.zmax0=max(data.in.var3);
     para.krig.dz0=(max(data.in.var3)-min(data.in.var3))/19;
end
if hdl.status.krigingfig == 1 & opt == 1
   set(hdl.krig.xmin,'string',para.krig.xmin0);
   set(hdl.krig.xmax,'string',para.krig.xmax0);
   set(hdl.krig.dx,'string',para.krig.dx0);
   set(hdl.krig.ymin,'string',para.krig.ymin0);
   set(hdl.krig.ymax,'string',para.krig.ymax0);
   set(hdl.krig.dy,'string',para.krig.dy0);
   if data.in.dim == 3
     set(hdl.krig.zmin,'string',para.krig.zmin0);
     set(hdl.krig.zmax,'string',para.krig.zmax0);
     set(hdl.krig.dz,'string',para.krig.dz0);
   else
     set(hdl.krig.zmin,'enable','off');
     set(hdl.krig.zmax,'enable','off');
     set(hdl.krig.dz,'enable','off');
  end
end
