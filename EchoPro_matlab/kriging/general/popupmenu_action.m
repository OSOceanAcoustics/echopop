function	popupmenu_action(window_index,popupmenu_index)
% popup menu action
% window_index = index of process task 
%            1 - Data Preparation
%            2 - Variogram
%            3 - Krig
%            2 - Visualization
% popupmenu_index = index for specific options
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

switch	window_index
	case 1					% data preperation window
      data_col={'Data Col. 1', 'Data Col. 2', 'Data Col. 3'};
      label_str={'LONGITUDE','LATITUDE','DEPTH','X','Y','TIME','Other'};
		unit_str={'(deg)','(km)','(m)','(cm)','(mm)','(Day)','(hour)','(min)','(sec)','(other)'};
		switch	popupmenu_index
  			case 1		% data column for the x-axis variable  
  			case 2		% data column for the y-axis variable  
  			case 3		% data column for the z-axis variable  
			case 4		% x-axis variable --> xlabel
				xlabel_indx=get(hdl.dataprep.var1,'Value');
			   set(hdl.dataprep.xlabel,'string',label_str(xlabel_indx));
				if xlabel_indx <= 6
					set(hdl.dataprep.x_unit,'style','popupmenu','string',unit_str);
					if xlabel_indx <= 2
						set(hdl.dataprep.x_unit,'value',1);
					elseif xlabel_indx == 3
						set(hdl.dataprep.x_unit,'value',3);
					elseif xlabel_indx <= 5
						set(hdl.dataprep.x_unit,'value',2);
					elseif xlabel_indx == 6
						set(hdl.dataprep.x_unit,'value',6);
					end
				else
					set(hdl.dataprep.x_unit,'style','edit','string','');
			   end
		   case 5		% ylabel
				ylabel_indx=get(hdl.dataprep.var2,'Value');
			   set(hdl.dataprep.ylabel,'string',label_str(ylabel_indx));
				if ylabel_indx <= 6
					set(hdl.dataprep.y_unit,'style','popupmenu','string',unit_str);
					if ylabel_indx <= 2
						set(hdl.dataprep.y_unit,'value',1);
					elseif ylabel_indx == 3
						set(hdl.dataprep.y_unit,'value',3);
					elseif ylabel_indx <= 5
						set(hdl.dataprep.y_unit,'value',2);
					elseif ylabel_indx == 6
						set(hdl.dataprep.y_unit,'value',6);
					end
				else
					set(hdl.dataprep.y_unit,'style','edit','string','');
			   end
		   case 6		% zlabel
				zlabel_indx=get(hdl.dataprep.var3,'Value');
			   set(hdl.dataprep.zlabel,'string',label_str(zlabel_indx));
				if zlabel_indx <= 6
					set(hdl.dataprep.z_unit,'style','popupmenu','string',unit_str);
					if zlabel_indx <= 2
						set(hdl.dataprep.z_unit,'value',1);
					elseif zlabel_indx == 3
						set(hdl.dataprep.z_unit,'value',3);
					elseif zlabel_indx <= 5
						set(hdl.dataprep.z_unit,'value',2);
					elseif zlabel_indx == 6
						set(hdl.dataprep.z_unit,'value',6);
					end
				else
					set(hdl.dataprep.z_unit,'style','edit','string','');
			   end
        	case 7		% data filter 
  				if para.status.dataprep == 1  dataprep3d(2);	end
         case 8		% data transformation 
				if para.status.dataprep == 1	dataprep3d(3);	end
		end
   case 2					% variogram/covariance window 
      vario_opt_button;
   case 3					% kriging wondow
      scheme=get(hdl.krig.scheme,'value');
      if scheme == 1				% point to point
         set(hdl.krig.blksize,'enable','off');
         set(hdl.krig.blk_X,'enable','off');
         set(hdl.krig.blk_nx,'enable','off');
         set(hdl.krig.blk_Y,'enable','off');
         set(hdl.krig.blk_ny,'enable','off');
         set(hdl.krig.blk_Z,'enable','off');
         set(hdl.krig.blk_nz,'enable','off');
      else
         set(hdl.krig.blksize,'enable','on');
         set(hdl.krig.blk_X,'enable','on');
         set(hdl.krig.blk_nx,'enable','on');
         set(hdl.krig.blk_Y,'enable','on');
         set(hdl.krig.blk_ny,'enable','on');
         set(hdl.krig.blk_Z,'enable','on');
         set(hdl.krig.blk_nz,'enable','on');
      end
	case 4					% visualization window
		if popupmenu_index == 1
			cross_validation(get(hdl.dispkrig.validationmodel,'value'),0);
		end
end