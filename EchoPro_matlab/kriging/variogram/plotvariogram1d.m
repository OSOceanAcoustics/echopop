function		plotVariogram1d(opt)
%% plot 1-D semi-variogram or correlogram 
%		opt = 1		plot data only
%             2     plot theoretical model only
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

subplot(hdl.vario.axes1);
switch opt
	case 1				% data only
      if para.status.variogram >= 1  
			hdl_vario1d=findobj(hdl.vario.axes1,'type','line');
			for i=1:length(hdl_vario1d)
				if strcmp(get(hdl_vario1d(i),'linestyle'),'none')
					delete(hdl.vario.data_plot); 
				elseif strcmp(get(hdl_vario1d(i),'linestyle'),'-')
					delete(hdl.vario.theo_plot); 
				end
			end
			hold on
		end
      if get(hdl.vario.variogram,'Value') == 1 
         hdl.vario.data_plot=plot(data.out.vario.lag,data.out.vario.gammah,'o');
         hy=ylabel('Semi-Variogram');
         if isfield(data.out.vario,'lag_theo')
            hold on
				hdl.vario.theo_plot=plot(data.out.vario.lag_theo,data.out.vario.gammah_theo,'r', ...
					'linewidth',1.5);
            hold off
            para.status.variogram=2;
        end
  		else
      	hdl.vario.data_plot=plot(data.out.vario.lag,1-data.out.vario.gammah,'o');
         if isfield(data.out.vario,'lag_theo')
            hold on
				hdl.vario.theo_plot=plot(data.out.vario.lag_theo,1-data.out.vario.gammah_theo,'r', ...
					'linewidth',1.5);
            hold off
            para.status.variogram=2;
         end
      	hy=ylabel('Correlogram');   
      end
      hold off
   case 2				% theoretical curve
      if para.status.variogram >= 1  
			hdl_vario1d=findobj(hdl.vario.axes1,'type','line');
			for i=1:length(hdl_vario1d)
				if strcmp(get(hdl_vario1d(i),'linestyle'),'-')
					delete(hdl.vario.theo_plot); 
				end
			end
		end
      hold on
		if get(hdl.vario.variogram,'Value') == 1
			hdl.vario.theo_plot=plot(data.out.vario.lag_theo,data.out.vario.gammah_theo,'r', ...
					'linewidth',1.5);
      	hy=ylabel('Semi-Variogram');
  		else
			hdl.vario.theo_plot=plot(data.out.vario.lag_theo,1-data.out.vario.gammah_theo,'r', ...
					'linewidth',1.5);
      	hy=ylabel('Correlogram');   
  		end
		hold off
   end
%para.vario.max_value=max(get(hdl.vario.data_plot,'Ydata'))*1.2;
para.vario.max_value=max(data.out.vario.gammah)*1.2;
hx=xlabel(para.vario.xlabel);
set([hx hy],'fontsize',12,'fontweight','bold');
para.vario.disp_range=get(hdl.vario.disp_range_slider,'Value')*para.vario.max_range;
set(hdl.vario.axes1,'xlim',[0 para.vario.disp_range]);
para.vario.disp_value=get(hdl.vario.disp_value_slider,'Value')*para.vario.max_value;
ymin=min(get(hdl.vario.data_plot,'Ydata'));
set(hdl.vario.axes1,'ylim',[ymin para.vario.max_value]);

set(hdl.vario.axes1,'xlim',[0 para.vario.disp_range0]);
disp_slider_value=para.vario.disp_range0/para.vario.max_range;
set(hdl.vario.disp_range_slider,'Value',disp_slider_value);

