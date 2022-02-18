function  get_set_gridfile_para
% get and set grid parameters based on the output- grid data
% from customized grid file
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para data hdl

if para.proc.kriging == 1
    dat=data.in.US_CAN_mesh_cor;
else
    dat=load(para.krig.grid_file);
end

if para.dataprep.x_axis_indx == 1 & para.dataprep.y_axis_indx == 2
    data.out.krig.gx=dat(:,1);
    data.out.krig.gy=dat(:,2);
elseif para.dataprep.x_axis_indx == 2 & para.dataprep.y_axis_indx == 1
    data.out.krig.gx=dat(:,2);
    data.out.krig.gy=dat(:,1);  
else
    disp('Grid Data file Format Inconsistent on line 20 in get_set_gridfile_para.m!');
end

%% only execute if on batch process or if it is the first time with batch process
%if para.krig.batch_file_proc == 0 |(para.krig.batch_file_proc == 1 & bat_proc_cnt == 1)
    
%% get parameters from the default values
xmin=str2num(get(hdl.krig.xmin,'string'));
xmax=str2num(get(hdl.krig.xmax,'string'));
dx=str2num(get(hdl.krig.dx,'string'));
ymin=str2num(get(hdl.krig.ymin,'string'));
ymax=str2num(get(hdl.krig.ymax,'string'));
dy=str2num(get(hdl.krig.dy,'string'));

gxmin=min(min(data.out.krig.gx),xmin);
gxmax=max(max(data.out.krig.gx),xmax);
gdx=min(dx,mean(data.out.krig.gx));
gymin=min(min(data.out.krig.gy,ymin));
gymax=max(max(data.out.krig.gy),ymax);
gdy=min(dy,mean(data.out.krig.gy));

set(hdl.krig.xmin,'string',num2str(gxmin));
set(hdl.krig.xmax,'string',num2str(gxmax));
%set(hdl.krig.dx,'string',num2str(gdx));
set(hdl.krig.ymin,'string',num2str(gymin));
set(hdl.krig.ymax,'string',num2str(gymax));
%set(hdl.krig.dy,'string',num2str(gdy));

if data.in.dim == 3
  data.out.krig.gz=dat(:,3);
  zmin=str2num(get(hdl.krig.zmin,'string'));
  zmax=str2num(get(hdl.krig.zmax,'string'));
  dz=str2num(get(hdl.krig.dz,'string'));
  gzmin=min(min(data.out.krig.gz,ymin));
  gzmax=max(max(data.out.krig.gy),ymax);
  gdz=min(dz,mean(data.out.krig.gz));
  set(hdl.krig.zmin,'string',num2str(gzmin));
  set(hdl.krig.zmax,'string',num2str(gzmax));
  para.krig.zmin0=zmin;
  para.krig.zmax0=zmax;
%  set(hdl.krig.dz,'string',num2str(gdz));
end
para.krig.xmin0=gxmin;
para.krig.xmax0=gxmax;
para.krig.ymin0=gymin;
para.krig.ymax0=gymax;

%end              
