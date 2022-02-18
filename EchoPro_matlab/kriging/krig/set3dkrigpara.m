function set3dkrigpara(opt)
% Set DeFault Kriging Parameters
%% opt = 1   set krig window parameter from krig variable struct
%% opt = 2	 set krig variable struct from krig window settings
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

if 0
if para.krig.batch_file_proc == 1  & para.krig.bat_proc_cnt > 1 
    return
end
end
    
if opt == 1				% set krig window parameter from krig variable struct but data preparation window is active
    %% coordinates
    set(hdl.krig.xmin,'string',para.krig.xmin0);
    set(hdl.krig.xmax,'string',para.krig.xmax0);
    set(hdl.krig.dx,'string',para.krig.dx0);
    set(hdl.krig.ymin,'string',para.krig.ymin0);
    set(hdl.krig.ymax,'string',para.krig.ymax0);
    set(hdl.krig.dy,'string',para.krig.dy0);
    if isfield(para.krig,'zmin0') & ~isempty(para.krig.zmin0) & data.in.dim == 3
        set(hdl.krig.zlabel,'enable','on');
        set(hdl.krig.zmin,'enable','on');
        set(hdl.krig.zmax,'enable','on');
        set(hdl.krig.dz,'enable','on');
        set(hdl.krig.zmin,'string',para.krig.zmin0);
        set(hdl.krig.zmax,'string',para.krig.zmax0);
        set(hdl.krig.dz,'string',para.krig.dz0);
    end
    if data.in.dim == 2
        para.krig.zmin0=[];
        para.krig.zmax0=[];
        para.krig.dz0=[];
        set(hdl.krig.zmin,'string','');
        set(hdl.krig.zmax,'string','');
        set(hdl.krig.dz,'string','');
        set(hdl.krig.zlabel,'enable','off');
        set(hdl.krig.zmin,'enable','off');
        set(hdl.krig.zmax,'enable','off');
        set(hdl.krig.dz,'enable','off');
    end
    
    %% kriging parameters
    set(hdl.krig.model,'value',para.krig.model);
    set(hdl.krig.scheme,'value',para.krig.scheme);
    set(hdl.krig.blk_nx,'string',para.krig.blk_nx);
    set(hdl.krig.blk_ny,'string',para.krig.blk_ny);
    set(hdl.krig.blk_nz,'string',para.krig.blk_nz);
    if para.status.kriging == 1
        %		para.krig.srad=0.3*sqrt(para.dataprep.x_norm^2+para.dataprep.y_norm^2);
        para.krig.srad=0.3;
    end
    set(hdl.krig.srad,'string',para.krig.srad);
    set(hdl.krig.kmin,'string',para.krig.kmin);
    set(hdl.krig.kmax,'string',para.krig.kmax);
    set(hdl.krig.elim,'string',para.krig.elim);
    
    %% loading/saving parameters file
    set(hdl.krig.load_para,'value',para.krig.load_para);
    set(hdl.krig.save_para,'value',para.krig.save_para);
    if para.krig.load_para == 0
        set(hdl.krig.vario_para,'value',para.krig.vario_para,'enable','off');
        set(hdl.krig.krig_para,'value',para.krig.krig_para,'enable','off');
        set(hdl.krig.both_para,'value',para.krig.both_para,'enable','off');
        if para.krig.save_para == 0
            set(hdl.krig.para_file_browser,'enable','off');
        else
            set(hdl.krig.para_file_browser,'enable','on');
        end
    end
    
    %% load data file
    set(hdl.krig.load_data_file,'value',para.krig.load_data_file);
    set(hdl.krig.data_file_browser,'enable','off');
    set(hdl.krig.data_file,'visible','off');
elseif opt == 2                                     % set krig parameters based on the defination on the krig panel
    xmin=str2num(get(hdl.krig.xmin,'string'));
    xmax=str2num(get(hdl.krig.xmax,'string'));
    dx=str2num(get(hdl.krig.dx,'string'));
    ymin=str2num(get(hdl.krig.ymin,'string'));
    ymax=str2num(get(hdl.krig.ymax,'string'));
    dy=str2num(get(hdl.krig.dy,'string'));
    para.krig.xmin0=xmin;
    para.krig.xmax0=xmax;
    para.krig.dx0=dx;
    para.krig.ymin0=ymin;
    para.krig.ymax0=ymax;
    para.krig.dy0=dy;
    if para.dataprep.var1_indx == 1 & para.dataprep.var2_indx == 2  %y-axis -> Lat
        para.krig.xmin=min((xmin-para.dataprep.x_offset)*para.dataprep.latlonfac/para.dataprep.x_norm);
        para.krig.xmax=max((xmax-para.dataprep.x_offset)*para.dataprep.latlonfac/para.dataprep.x_norm);
        para.krig.ymin=(ymin-para.dataprep.y_offset)/para.dataprep.y_norm;
        para.krig.ymax=(ymax-para.dataprep.y_offset)/para.dataprep.y_norm;
    elseif para.dataprep.var1_indx == 2 & para.dataprep.var2_indx == 1  %x-axis -> Lat
        para.krig.xmin=(xmin-para.dataprep.x_offset)/para.dataprep.x_norm;
        para.krig.xmax=(xmax-para.dataprep.x_offset)/para.dataprep.x_norm;
        para.krig.ymin=min((ymin-para.dataprep.y_offset)*para.dataprep.latlonfac/para.dataprep.y_norm);
        para.krig.ymax=max((ymax-para.dataprep.y_offset)*para.dataprep.latlonfac/para.dataprep.y_norm);
    else
        para.krig.xmin=(xmin-para.dataprep.x_offset)/para.dataprep.x_norm;
        para.krig.xmax=(xmax-para.dataprep.x_offset)/para.dataprep.x_norm;
        para.krig.ymin=(ymin-para.dataprep.y_offset)/para.dataprep.y_norm;
        para.krig.ymax=(ymax-para.dataprep.y_offset)/para.dataprep.y_norm;
    end
    para.krig.dx=dx/para.dataprep.x_norm;
    para.krig.dy=dy/para.dataprep.y_norm;
    if data.in.dim == 3
        zmin=str2num(get(hdl.krig.zmin,'string'));
        zmax=str2num(get(hdl.krig.zmax,'string'));
        dz=str2num(get(hdl.krig.dz,'string'));
        para.krig.zmin0=zmin;
        para.krig.zmax0=zmax;
        para.krig.dz0=dz;
        para.krig.zmin=(zmin-para.dataprep.z_offset)/para.dataprep.z_norm;
        para.krig.zmax=(zmax-para.dataprep.z_offset)/para.dataprep.z_norm;
        para.krig.dz=dz/para.dataprep.z_norm;
    else
        
    end
    para.krig.model=get(hdl.krig.model,'value');
    para.krig.scheme=get(hdl.krig.scheme,'value');
    if para.krig.scheme == 2
        para.krig.blk_nx=str2num(get(hdl.krig.blk_nx,'string'));
        para.krig.blk_ny=str2num(get(hdl.krig.blk_ny,'string'));
        if data.in.dim == 3
            para.krig.blk_nz=str2num(get(hdl.krig.blk_nz,'string'));
        end
    end
    para.krig.srad=str2num(get(hdl.krig.srad,'string'));
    para.krig.kmin=str2num(get(hdl.krig.kmin,'string'));
    para.krig.kmax=str2num(get(hdl.krig.kmax,'string'));
    para.krig.elim=str2num(get(hdl.krig.elim,'string'));
end
