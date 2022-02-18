function		advanced_dispkrig3d(option,index1,index2,index3)
% function		advanced_dispkrig3d(option,index1,index2,index3) sets track lines and color map 
% option  = 1	 trackline
%		index1 = 1	color-coded trakline
%				 2	black/white trackline
%				 3	attribute values
%				 4	difference between observed and predicted attribute values
%                5  none
%		index2 = 	color index
%		index3 =	font size
% option = 2	 colormap
% advanced display operation
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

para.dispkrig.trackline.dispflag=1;

switch	option	
	case 1			% track line
	  	set(eval(sprintf('hdl.dispkrig3d.trackline.c%g',para.dispkrig.trackline.type_indx)),'checked','off');
		set(eval(sprintf('hdl.dispkrig3d.trackline.c%g',index1)),'checked','on');
		para.dispkrig.trackline.type_indx=index1;
        if nargin == 2
            switch index1
                case 5
                   para.dispkrig.trackline.dispflag=0;
                otherwise
                   para.dispkrig.trackline.dispflag=1;
            end
        elseif nargin == 3
	   	  set(eval(sprintf('hdl.dispkrig3d.trackline.c2color%g',para.dispkrig.trackline.line_color)),'checked','off');
	      set(eval(sprintf('hdl.dispkrig3d.trackline.c2color%g',index2)),'checked','on');
		  para.dispkrig.trackline.line_color=index2;
	    elseif nargin == 4
		  if index2 == 3			% color
				% value
	   		set(eval(sprintf('hdl.dispkrig3d.trackline.c3color%g',para.dispkrig.trackline.color_indx)),'checked','off');
	      	set(eval(sprintf('hdl.dispkrig3d.trackline.c3color%g',index3)),'checked','on');
	  			% difference
		 		set(eval(sprintf('hdl.dispkrig3d.trackline.c4color%g',para.dispkrig.trackline.color_indx)),'checked','off');
	      	set(eval(sprintf('hdl.dispkrig3d.trackline.c4color%g',index3)),'checked','on');
				para.dispkrig.trackline.color_indx=index3;
		  elseif index2 == 4		% size
				% value
	   		set(eval(sprintf('hdl.dispkrig3d.trackline.c3size%g',para.dispkrig.trackline.size_indx)),'checked','off');
	      	set(eval(sprintf('hdl.dispkrig3d.trackline.c3size%g',index3)),'checked','on');
				% difference
	   		set(eval(sprintf('hdl.dispkrig3d.trackline.c4size%g',para.dispkrig.trackline.size_indx)),'checked','off');
	      	set(eval(sprintf('hdl.dispkrig3d.trackline.c4size%g',index3)),'checked','on');
			para.dispkrig.trackline.size_indx=index3;
		  end
		end
		dispkrig3d(6);
	case 2
		if index1 < 18
	     clrmap_str=get(eval(sprintf('hdl.dispkrig3d.colormap.c%g',index1)),'label');
		  set(eval(sprintf('hdl.dispkrig3d.colormap.c%g',para.dispkrig.colormap_indx)),'checked','off');
		  set(eval(sprintf('hdl.dispkrig3d.colormap.c%g',index1)),'checked','on');
		  para.dispkrig.colormap_indx=index1;
		  eval(['colormap(' clrmap_str ');']);
		  dispkrig3d(6);
		else							% customized colormap
        file_browser3d(4,2);
		end
end
		