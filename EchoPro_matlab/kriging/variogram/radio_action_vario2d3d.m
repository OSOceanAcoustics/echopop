function radio_action_vario2d3d(index)
% function radio_action_vario2d3d(index) process the display options
% for the 2D/3D semi-variogram/correlogram vidualization
%   index = 1  azimuth angle
%         = 2  dip angle
%         = 3  shading faceted
%         = 4  shading flat
%         = 5  shading interp
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

switch index
 case 1				% adjust diaplay azimuth angle
   set(hdl.dispvario2d3d.dip_radio,'value',0);
   set(hdl.dispvario2d3d.dip_val,'string','-90 -- 90 deg');
   set(hdl.dispvario2d3d.azm_val,'enable','on')
   set(hdl.dispvario2d3d.azm_slider,'enable','on')
   mean_azm_ang=(para.vario.azm_beg+para.vario.azm_end)/2;
   set(hdl.dispvario2d3d.azm_slider,'value',0.5)
   set(hdl.dispvario2d3d.azm_val,'string',[num2str(mean_azm_ang)])
   set(hdl.dispvario2d3d.dip_val,'enable','off')
   set(hdl.dispvario2d3d.dip_slider,'enable','off')
   plotvariogram2d3d(2)
 case 2				% adjust diaplay dip angle
   set(hdl.dispvario2d3d.azm_radio,'value',0);
   set(hdl.dispvario2d3d.azm_val,'string','0 -- 360 deg');
   set(hdl.dispvario2d3d.dip_val,'enable','on')
   set(hdl.dispvario2d3d.dip_slider,'enable','on')
   mean_dip_ang=(para.vario.dip_beg+para.vario.dip_end)/2;
   set(hdl.dispvario2d3d.dip_slider,'value',0.5)
   set(hdl.dispvario2d3d.dip_val,'string',[num2str(mean_dip_ang)])
   set(hdl.dispvario2d3d.azm_val,'enable','off')
   set(hdl.dispvario2d3d.azm_slider,'enable','off')
   plotvariogram2d3d(3)
 case 3
   shading faceted
   set(hdl.dispvario2d3d.shading_radio2,'value',0);
   set(hdl.dispvario2d3d.shading_radio3,'value',0);
 case 4
	shading flat
   set(hdl.dispvario2d3d.shading_radio1,'value',0);
   set(hdl.dispvario2d3d.shading_radio3,'value',0);
 case 5
   shading interp
   set(hdl.dispvario2d3d.shading_radio1,'value',0);
   set(hdl.dispvario2d3d.shading_radio2,'value',0);
end

return