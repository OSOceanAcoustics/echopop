function		radio_action(window_index,option_index)
% action corresponding to the radio selection
% window_index = index of process task 
%            1 - Data Preparation
%            2 - Variogram
%            3 - Krig
%            2 - Visualization
% option_index = index for specific options
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

switch window_index
	case 1					% data preparation window
		switch(option_index)  
 		  case 1				% external data processing program
			if get(hdl.dataprep.ext_prog,'value');
              set(hdl.dataprep.dat_conv,'enable','on')
              set(hdl.dataprep.dat_conv_fname,'visible','on');
     		  para.dataprep.ext_prog=1;
              set(hdl.dataprep.x_axis,'enable','off');
              set(hdl.dataprep.y_axis,'enable','off');
              set(hdl.dataprep.z_axis,'enable','off');
       	    else
              set(hdl.dataprep.dat_conv,'enable','off')
              set(hdl.dataprep.dat_conv_fname,'visible','off');
			  para.dataprep.ext_prog=0;
              set(hdl.dataprep.x_axis,'enable','on');
              set(hdl.dataprep.y_axis,'enable','on');
              set(hdl.dataprep.z_axis,'enable','on');
            end
		  case 2          % x-axis direction
			 if get(hdl.dataprep.xdir,'value');
			   set(hdl.dataprep.axes1,'xdir','reverse');
             else
               set(hdl.dataprep.axes1,'xdir','normal');
			 end
		  case 3          % y-axis direction
			 if get(hdl.dataprep.ydir,'value');
			   set(hdl.dataprep.axes1,'ydir','reverse');
             else
               set(hdl.dataprep.axes1,'ydir','normal');
			 end
		  case 4          % z-axis direction
			 if get(hdl.dataprep.zdir,'value');
			   set(hdl.dataprep.axes1,'zdir','reverse');
             else
               set(hdl.dataprep.axes1,'zdir','normal');
			 end
          case 5            % save data format to a file
              if get(hdl.dataprep.data_format,'value')
                  set(hdl.dataprep.data_format_browser,'enable','on');
              else
                  set(hdl.dataprep.data_format_browser,'enable','off');
              end
          case 6
              if get(hdl.dataprep.data_type1,'value') == 1
                  set(hdl.dataprep.data_type2,'value',0);
                  para.dataprep.data_disptype=1;   % current setting
                  dataprep3d(2);
              end
          case 7
              if get(hdl.dataprep.data_type2,'value') == 1
                  set(hdl.dataprep.data_type1,'value',0);
                  para.dataprep.data_disptype=2;   % current setting
                  dataprep3d(2);
            end                 
	   end
	case 2		% semi-variogram/correlogram window
		switch option_index 
			case 1								% select variogram 
				if get(hdl.vario.variogram,'value')
					set(hdl.vario.correlogram,'value',0);
				else
					set(hdl.vario.correlogram,'value',1);
				end	
            if para.status.dataprep == 1
               plotvariogram1d(1);
            elseif para.status.varigram == 1
               plotvariogram2d3d(1);
               variogram_theo(1); 
            end
			case 2								% select correlogram
				if get(hdl.vario.correlogram,'value')
					set(hdl.vario.variogram,'value',0);
				else
					set(hdl.vario.variogram,'value',1);
				end
            if para.status.dataprep == 1
               plotvariogram1d(1);
            elseif para.status.variogram == 1
				  corr_y_min=1-max(get(hdl.vario.axes1,'ylim'));
              plotvariogram1d(1);
              variogram_theo(1); 
				  corr_y_max=max(get(hdl.vario.axes1,'ylim'));
				  set(hdl.vario.axes1,'ylim',[corr_y_min corr_y_max] );
            end
			case 3							% load parameter file
				if get(hdl.vario.load_para_file,'value');
					para.vario.load_para=1;
					set(hdl.vario.load_file_browser,'enable','on');
                    set(hdl.vario.save_para_file,'value',0);
                    set(hdl.vario.save_file_browser,'enable','off');
				else
					set(hdl.vario.load_file_browser,'enable','off');
					para.vario.load_para=0;
				end	
           case 4                          % save parameters to a file
				if get(hdl.vario.save_para_file,'value');
					para.vario.save_para=1;
					set(hdl.vario.save_file_browser,'enable','on');
                    set(hdl.vario.load_para_file,'value',0);
                    set(hdl.vario.load_file_browser,'enable','off');
				else
					set(hdl.vario.save_file_browser,'enable','off');
					para.vario.save_para=0;
				end	
           case 5							% Enable 1-D variogram - isotropy
               para.vario.anisotropy=0;
 			   set(hdl.vario.enable2d3d,'value',0);
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
               set(hdl.vario.ang_azm_2d3d,'enable','off');
               set(hdl.vario.ang_dip_2d3d,'enable','off');
           case 6							% Enable 2-D variogram  - anisotropy
                para.vario.anisotropy=1;
 				set(hdl.vario.enabled,'value',0);
                set(hdl.vario.ang_dip_2d3d,'enable','on');
                set(hdl.vario.ang_azm_2d3d,'enable','on');
                if get(hdl.vario.ang_azm_2d3d,'value')
  				   set(hdl.vario.ang_beg_azm,'enable','on','BackgroundColor',[1 1 1]);
 				   set(hdl.vario.ang_end_azm,'enable','on','BackgroundColor',[1 1 1]);
 				   set(hdl.vario.ang_res_azm,'enable','on','BackgroundColor',[1 1 1]);
 				   set(hdl.vario.ang_wd_azm,'enable','on','BackgroundColor',[1 1 1]);
 				   set(hdl.vario.ang_rot_azm,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ytox_ratio,'enable','on','BackgroundColor',[1 1 1]);
                   if data.in.dim == 3
                     set(hdl.vario.ang_dip_2d3d,'enable','on');
                     set(hdl.vario.ang_beg_dip,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ang_end_dip,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ang_res_dip,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ang_wd_dip,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ang_rot_dip,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ztox_ratio,'enable','on','BackgroundColor',[1 1 1]);
                   else
                     set(hdl.vario.ang_dip_2d3d,'enable','off');
                   end
                else 
       if 0                 % not effective at this time
                   set(hdl.vario.ang_beg_dip,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ang_end_dip,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ang_res_dip,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ang_wd_dip,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ang_rot_dip,'enable','on','BackgroundColor',[1 1 1]);
                   set(hdl.vario.ztox_ratio,'enable','on','BackgroundColor',[1 1 1]);
                   if data.in.dim == 3
 				     set(hdl.vario.ang_beg_azm,'enable','on','BackgroundColor',[1 1 1]);
 				     set(hdl.vario.ang_end_azm,'enable','on','BackgroundColor',[1 1 1]);
 				     set(hdl.vario.ang_res_azm,'enable','on','BackgroundColor',[1 1 1]);
 				     set(hdl.vario.ang_wd_azm,'enable','on','BackgroundColor',[1 1 1]);
 				     set(hdl.vario.ang_rot_azm,'enable','on','BackgroundColor',[1 1 1]);
                     set(hdl.vario.ytox_ratio,'enable','on','BackgroundColor',[1 1 1]);
                 end
        end
               end 
			%	set3dvariopara(3);
           case 7							%  anisotropic case - select azimuth plane as the primary plane
               set(hdl.vario.ang_dip_2d3d,'value',0);
  			   set(hdl.vario.ang_beg_azm,'enable','on','BackgroundColor',[1 1 1]);
 			   set(hdl.vario.ang_end_azm,'enable','on','BackgroundColor',[1 1 1]);
 			   set(hdl.vario.ang_res_azm,'enable','on','BackgroundColor',[1 1 1]);
 			   set(hdl.vario.ang_wd_azm,'enable','on','BackgroundColor',[1 1 1]);
 			   set(hdl.vario.ang_rot_azm,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ytox_ratio,'enable','on','BackgroundColor',[1 1 1]);
      %         set(hdl.vario.ang_beg_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
      %         set(hdl.vario.ang_end_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
      %         set(hdl.vario.ang_res_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
      %         set(hdl.vario.ang_wd_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
      %         set(hdl.vario.ang_rot_dip,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
      %         set(hdl.vario.ztox_ratio,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
               if data.in.dim == 3
                  set(hdl.vario.ang_beg_dip,'enable','on','BackgroundColor',[1 1 1]);
                  set(hdl.vario.ang_wd_dip,'enable','on','BackgroundColor',[1 1 1]);
                  set(hdl.vario.ztox_ratio,'enable','on','BackgroundColor',[1 1 1]);
                  set(hdl.vario.ang_end_dip,'enable','on','BackgroundColor',[1 1 1]);
                  set(hdl.vario.ang_res_dip,'enable','on','BackgroundColor',[1 1 1]);
                 set(hdl.vario.ang_rot_dip,'enable','on','BackgroundColor',[1 1 1]);
              end
           case 8
               set(hdl.vario.ang_azm_2d3d,'value',0);
               set(hdl.vario.ang_beg_dip,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ang_end_dip,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ang_res_dip,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ang_wd_dip,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ang_rot_dip,'enable','on','BackgroundColor',[1 1 1]);
               set(hdl.vario.ztox_ratio,'enable','on','BackgroundColor',[1 1 1]);
 %  			   set(hdl.vario.ang_beg_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
 %			   set(hdl.vario.ang_end_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
 %			   set(hdl.vario.ang_res_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
 %			   set(hdl.vario.ang_wd_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
 %			   set(hdl.vario.ang_rot_azm,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
 %              set(hdl.vario.ytox_ratio,'enable','off','BackgroundColor',[0.65 0.65 0.65]);
               if data.in.dim == 3
 				  set(hdl.vario.ang_beg_azm,'enable','on','BackgroundColor',[1 1 1]);
  			      set(hdl.vario.ang_end_azm,'enable','on','BackgroundColor',[1 1 1]);
 			      set(hdl.vario.ang_res_azm,'enable','on','BackgroundColor',[1 1 1]);
				  set(hdl.vario.ang_wd_azm,'enable','on','BackgroundColor',[1 1 1]);
 			      set(hdl.vario.ang_rot_azm,'enable','on','BackgroundColor',[1 1 1]);
                  set(hdl.vario.ytox_ratio,'enable','on','BackgroundColor',[1 1 1]);
               end 
		   end
	case 3				% kriging window
		switch option_index 
			case 1			% load parameter file
                set(hdl.krig.save_para,'value',0);
				if get(hdl.krig.load_para,'value') 
					para.krig.load_para=1;
					para.krig.save_para=0;
				    set(hdl.krig.save_para,'value',0);
					set(hdl.krig.para_file_browser,'enable','on');
					set(hdl.krig.vario_para,'value',para.krig.vario_para,'enable','on');
					set(hdl.krig.krig_para,'value',para.krig.krig_para,'enable','on');
					set(hdl.krig.both_para,'value',para.krig.both_para,'enable','on');
				else
					para.krig.load_para=0;
					set(hdl.krig.para_file_browser,'enable','off');
					set(hdl.krig.vario_para,'value',para.krig.vario_para,'enable','off');
					set(hdl.krig.krig_para,'value',para.krig.krig_para,'enable','off');
					set(hdl.krig.both_para,'value',para.krig.both_para,'enable','off');
				end	
			case 2		% save parameter file
                set(hdl.krig.load_para,'value',0);
				if  get(hdl.krig.save_para,'value')
					para.krig.save_para=1;
					para.krig.load_para=0;
					set(hdl.krig.load_para,'value',0);
					set(hdl.krig.para_file_browser,'enable','on');
					set(hdl.krig.vario_para,'value',para.krig.vario_para,'enable','off');
					set(hdl.krig.krig_para,'value',para.krig.krig_para,'enable','off');
					set(hdl.krig.both_para,'value',para.krig.both_para,'enable','off');
				else
					set(hdl.krig.para_file_browser,'enable','off');
					para.krig.save_para=0;
				end	
 			case 3		% variogram parameter file only
				if  get(hdl.krig.vario_para,'value')
					para.krig.vario_para=1;
					para.krig.krig_para=0;
					para.krig.both_para=0;
					set(hdl.krig.krig_para,'value',0);
					set(hdl.krig.both_para,'value',0);
				end
			case 4		% kriging parameter file only
				if get(hdl.krig.krig_para,'value')
					para.krig.vario_para=0;
					para.krig.krig_para=1;
					para.krig.both_para=0;
					set(hdl.krig.vario_para,'value',0);
					set(hdl.krig.both_para,'value',0);
				end
			case 5		% both variogram and kriging parameter files
				if get(hdl.krig.both_para,'value')
					para.krig.vario_para=0;
					para.krig.krig_para=0;
					para.krig.both_para=1;
					set(hdl.krig.vario_para,'value',0);
					set(hdl.krig.krig_para,'value',0);
				end
			case 6		% load data file
                if hdl.status.dataprepfig ~= 1 & isempty(para.krig.data_format_file)
                    krig_message(1,'Need a Data Format File first or Load a data file in Data Preparation Window !!');
                    set(hdl.krig.load_data_file,'value',0)
                    return
                end
				if get(hdl.krig.load_data_file,'value') 
					set(hdl.krig.data_file_browser,'enable','on');
					set(hdl.krig.data_file,'visible','on');
				else
					para.krig.load_data_file=0;
					set(hdl.krig.data_file_browser,'enable','off');
					set(hdl.krig.data_file,'visible','off');
				end	
 			case 7  		% batch file processing
 				if get(hdl.krig.batch_file_proc,'value') 
					para.krig.batch_file_proc=0;
					set(hdl.krig.krig_button,'enable','off');
					set(hdl.krig.batch_data_file_browser,'enable','on');
					set(hdl.krig.batch_data_file,'enable','on');
					set(hdl.krig.batch_log_file_browser,'enable','on');
					set(hdl.krig.batch_log_file,'enable','on');
					set(hdl.krig.batch_krig,'enable','on');
				else
					para.krig.batch_file_proc=0;
					set(hdl.krig.krig_button,'enable','on');
					set(hdl.krig.batch_data_file_browser,'enable','off');
					set(hdl.krig.batch_data_file,'enable','off');
					set(hdl.krig.batch_log_file_browser,'enable','off');
					set(hdl.krig.batch_log_file,'enable','off');
					set(hdl.krig.batch_krig,'enable','off');
				end	
            case 8      % custormized grid points
 				   if get(hdl.krig.customized_gridfile,'value') 
					   set(hdl.krig.gridfile_browser,'enable','on');
				   else
					   para.krig.load_griddata_file=0;
                       para.krig.bat_proc_cnt=0;
                       para.krig.grid_file='';
					   set(hdl.krig.gridfile_browser,'enable','off');
				   end	  
            case 9      % load data format file
				if get(hdl.krig.load_data_format_file,'value') 
					set(hdl.krig.data_format_file_browser,'enable','on');
					set(hdl.krig.data_file,'visible','on');
				else
					para.krig.load_data_format_file=0;
					set(hdl.krig.data_format_file_browser,'enable','off');
					set(hdl.krig.data_file,'visible','off');
				end	
			end
	case 4			%visualization window			
		if option_index == 1			% data forward/inverse transformation
			para.dispkrig.Qcheck=0;
			para.dispkrig.JKcheck=0;
			dispkrig(1)
		end
end