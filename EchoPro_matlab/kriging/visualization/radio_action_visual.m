function radio_action_visual(index)
%%  function radio_action_visual(index) processes visualization options
%%  
%%   index =  1   Kriging map
%%            2   Kriging variance map
%%            3   Cross validation
%%            4   Shading Faceted
%%            5   Shading Flat
%%            6   Shading interp
%%            7   X-axis direction - normal/reverse
%%            8   Y-axis direction - normal/reverse
%%            9   Z-axis direction - normal/reverse
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl data

if index <= 3
   hdl.dispkrig3d.var_index=index;
elseif index <= 6
   hdl.dispkrig3d.shading_index=index-3;
end
if data.in.dim == 3
   var_index=0;
else
   var_index=5;
end

switch index
 case 1
   set(hdl.dispkrig3d.disp_var2,'value',0);
   set(hdl.dispkrig3d.disp_var3,'value',0);
   dispkrig3d(var_index);
 case 2
   set(hdl.dispkrig3d.disp_var1,'value',0);
   set(hdl.dispkrig3d.disp_var3,'value',0);
   dispkrig3d(var_index);
 case 3
   set(hdl.dispkrig3d.disp_var1,'value',0);
   set(hdl.dispkrig3d.disp_var2,'value',0);
   validationfig;
 case 4
   shading faceted
   set(hdl.dispkrig3d.shading_radio2,'value',0);
   set(hdl.dispkrig3d.shading_radio3,'value',0);
 case 5
	shading flat
   set(hdl.dispkrig3d.shading_radio1,'value',0);
   set(hdl.dispkrig3d.shading_radio3,'value',0);
 case 6
   shading interp
   set(hdl.dispkrig3d.shading_radio1,'value',0);
   set(hdl.dispkrig3d.shading_radio2,'value',0);
case 7
   if get(hdl.dispkrig3d.xdir_reverse,'value') == 1
      set(hdl.dispkrig3d.axes1,'xdir','reverse')
   else
      set(hdl.dispkrig3d.axes1,'xdir','normal')
	end      
case 8
   if get(hdl.dispkrig3d.ydir_reverse,'value') == 1
      set(hdl.dispkrig3d.axes1,'ydir','reverse')
   else
      set(hdl.dispkrig3d.axes1,'ydir','normal')
	end      
case 9
   if get(hdl.dispkrig3d.zdir_reverse,'value') == 1
      set(hdl.dispkrig3d.axes1,'zdir','reverse')
   else
      set(hdl.dispkrig3d.axes1,'zdir','normal')
	end      
end

return