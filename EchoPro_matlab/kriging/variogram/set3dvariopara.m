function set3dvariopara(opt)
%% function set3dvariopara(opt) sets variogram/correlogram parameters 
%% based on para.vario struct array
%% opt = 0          set default anisotropy display
%% opt = 1			set slider position only
%% opt = 2			set slider position, edit field and all other settings
%% opt = 3			set anisotropy parameters only (from panel)
%% opt = 4			set anisotropy parameters only (from parameter struct)
%% opt = 5          set all parameters from the parameter struct
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl

if opt == 0
   set(hdl.vario.ang_beg_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_end_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_res_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_wd_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_rot_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_beg_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_end_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_res_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_wd_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ang_rot_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ytox_ratio,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   set(hdl.vario.ztox_ratio,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
   return
end
%% change slider positions, load parameter files, default setings
set(hdl.vario.range_slider,'Value',para.vario.range/para.vario.max_range);
set(hdl.vario.nugt_slider,'Value',para.vario.nugt/para.vario.max_nugt);
set(hdl.vario.sill_slider,'Value',para.vario.sill/para.vario.max_sill);
set(hdl.vario.lscl_slider,'Value',para.vario.lscl/para.vario.max_lscl);
set(hdl.vario.powr_slider,'Value',para.vario.powr/para.vario.max_powr);
set(hdl.vario.hole_slider,'Value',para.vario.hole/para.vario.max_hole);
set(hdl.vario.disp_range_slider,'Value',para.vario.disp_range/para.vario.max_range);
set(hdl.vario.disp_value_slider,'Value',para.vario.disp_value/para.vario.max_value);
if opt == 1
   return
end
if opt >= 2				% load parameter files, default setings
	set(hdl.vario.variogram,'Value',para.vario.vario);
    set(hdl.vario.correlogram,'Value',para.vario.corr);
    set(hdl.vario.model,'Value',para.vario.model);
	set(hdl.vario.nugt_edit,'String',num2str(para.vario.nugt));
	set(hdl.vario.range_edit,'String',num2str(para.vario.range));
	set(hdl.vario.sill_edit,'String',num2str(para.vario.sill));
	set(hdl.vario.lscl_edit,'String',num2str(para.vario.lscl));
	set(hdl.vario.powr_edit,'String',num2str(para.vario.powr));
	set(hdl.vario.hole_edit,'String',num2str(para.vario.hole));
	set(hdl.vario.res,'String',num2str(para.vario.res));
	set(hdl.vario.load_para_file,'value',0);
    set(hdl.vario.load_file_browser,'enable','off');
	set(hdl.vario.ang_beg_azm,'string',para.vario.azm_beg);
end
if opt < 3
   return
else	%  anisotropy parameters only
    if para.vario.anisotropy ~= 1 return;end
  if 0
	if ~isempty(findobj(hdl.vario.h0,'Tag','variogramcolorbar'))
			set(hdl.vario.plot2d3d,'visible','on');
			set(hdl.vario.colorbar2,'visible','on');
			set(hdl.vario.htxt2d,'visible','on');
         set(get(hdl.vario.colorbar2,'children'),'visible','on');
			subplot(hdl.vario.axes2);
			axis off;
    end
  end
	if opt == 3				% get parameters from the panel
   %% Azimuth
	  para.vario.azm_beg=str2num(get(hdl.vario.ang_beg_azm,'string'));
	  para.vario.azm_end=str2num(get(hdl.vario.ang_end_azm,'string'));
	  para.vario.azm_res=str2num(get(hdl.vario.ang_res_azm,'string'));
      para.vario.ang_wd_azm=str2num(get(hdl.vario.ang_wd_azm,'string'));
      para.vario.ang_rot_azm=str2num(get(hdl.vario.ang_rot_azm,'string'));
   % Elevation
	  para.vario.dip_beg=str2num(get(hdl.vario.ang_beg_dip,'string'));
	  para.vario.dip_end=str2num(get(hdl.vario.ang_end_dip,'string'));
	  para.vario.dip_res=str2num(get(hdl.vario.ang_res_dip,'string'));
      para.vario.ang_wd_dip=str2num(get(hdl.vario.ang_wd_dip,'string'));
      para.vario.ang_rot_dip=str2num(get(hdl.vario.ang_rot_dip,'string'));
   %% Ratios
	  para.vario.ytox_ratio=str2num(get(hdl.vario.ytox_ratio,'string'));
	  para.vario.ztox_ratio=str2num(get(hdl.vario.ztox_ratio,'string'));
	  para.vario.res=str2num(get(hdl.vario.res,'string'));
    else		% set anisotropy parameters from the parameter structure
   %% Azimuth
	  set(hdl.vario.ang_beg_azm,'string',num2str(para.vario.azm_beg));
	  set(hdl.vario.ang_end_azm,'string',num2str(para.vario.azm_end));
	  set(hdl.vario.ang_res_azm,'string',num2str(para.vario.azm_res));
      set(hdl.vario.ang_wd_azm,'string',num2str(para.vario.ang_wd_azm));
      set(hdl.vario.ang_rot_azm,'string',num2str(para.vario.ang_rot_azm));
   % Elevation
      set(hdl.vario.ang_beg_dip,'string',num2str(para.vario.dip_beg));
	  set(hdl.vario.ang_end_dip,'string',num2str(para.vario.dip_end));
	  set(hdl.vario.ang_res_dip,'string',num2str(para.vario.dip_res));
      set(hdl.vario.ang_wd_dip,'string',num2str(para.vario.ang_wd_dip));
      set(hdl.vario.ang_rot_dip,'string',num2str(para.vario.ang_rot_dip));
   %% Ratios
      set(hdl.vario.ytox_ratio,'string',num2str(para.vario.ytox_ratio));
	  set(hdl.vario.ztox_ratio,'string',num2str(para.vario.ztox_ratio));
	  set(hdl.vario.res,'string',num2str(para.vario.res));    
	end
	return
end




   