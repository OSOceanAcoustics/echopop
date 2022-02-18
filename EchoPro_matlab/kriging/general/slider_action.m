function   slider_action(window_index,option)
%%  change corresponding edit field based on the slider position and update theoretical
%%  		modeling curve.
% window_index = index of process task 
%            1 - Data Preparation
%            2 - Variogram
%            3 - Krig
%            2 - Visualization
% option = index for specific options
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl 

switch window_index
	case 2				% variogram / correlogram window
		switch option
			case 1			% nugget
				para.vario.nugt=para.vario.max_nugt*get(hdl.vario.nugt_slider,'value');
				set(hdl.vario.nugt_edit,'String',num2str(para.vario.nugt));			
			case 2			% sill
				para.vario.sill=para.vario.max_sill*get(hdl.vario.sill_slider,'value');
				set(hdl.vario.sill_edit,'String',num2str(para.vario.sill));
			case 3			% length scale
				para.vario.lscl=para.vario.max_lscl*get(hdl.vario.lscl_slider,'value');
				set(hdl.vario.lscl_edit,'String',num2str(para.vario.lscl));
			case 4			% power
				para.vario.powr=para.vario.max_powr*get(hdl.vario.powr_slider,'value');
				set(hdl.vario.powr_edit,'String',num2str(para.vario.powr));
			case 5			% hole effect
				para.vario.hole=para.vario.max_hole*get(hdl.vario.hole_slider,'value');
				set(hdl.vario.hole_edit,'String',num2str(para.vario.hole));
			case 6			% range
				para.vario.range=para.vario.max_range*get(hdl.vario.range_slider,'value');
				set(hdl.vario.range_edit,'String',num2str(para.vario.range));
			case 7			% display range
				para.vario.disp_range=get(hdl.vario.disp_range_slider,'Value')*para.vario.max_range;
				set(hdl.vario.axes1,'xlim',[0 para.vario.disp_range]);
				return
			case 8			% display value
				para.vario.disp_value=get(hdl.vario.disp_value_slider,'Value')*para.vario.max_value;
				set(hdl.vario.axes1,'ylim',[min(get(hdl.vario.axes1,'ylim')) para.vario.disp_value]);
				return
		end
end     
if para.status.variogram == 2	% data-based semi-variogram/correlogram has been computed
	variogram_theo(1);				% get other parameters and update the theoretical curve  
end
