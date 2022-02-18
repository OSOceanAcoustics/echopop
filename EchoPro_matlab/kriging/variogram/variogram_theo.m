function 			variogram_theo(opt)
%% compute theoretical semi-variogram/correlogram using specified parameters on the panel
%%	opt = 1		get model parameters from variogram window panel
%%		  2		get model parameters from variable struct
%%
%%  Kriging Software Package  version 3.0,   December 29, 2001
%%  Copyright (c) 1998, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.
global hdl para data

if hdl.status.variogramfig == 1			% window is open
   %% Data-based Semi-Variogram or Correlogram has been computed
	if opt == 1
		[model,nugt,sill,lscl,powr,hole]=get3dvariopara_theo;
   else
		model=para.vario.model;
		nugt=para.vario.nugt;
		sill=para.vario.sill;
		lscl=para.vario.lscl;
		powr=para.vario.powr;
		hole=para.vario.hole;
		set(hdl.vario.nugt_edit,'String',num2str(para.vario.nugt));
		set(hdl.vario.sill_edit,'String',num2str(para.vario.sill));
		set(hdl.vario.lscl_edit,'String',num2str(para.vario.lscl));
		set(hdl.vario.powr_edit,'String',num2str(para.vario.powr));
		set(hdl.vario.hole_edit,'String',num2str(para.vario.hole));
      get3dvariopara_theo;
	end
	model_para=[nugt sill lscl powr hole];
	set3dvariopara(1);							% set slider positions only
	lag=data.out.vario.lag;
	range=str2num(get(hdl.vario.range_edit,'String'));
	indx=find(lag <= range);
	data.out.vario.lag_theo=data.out.vario.lag(indx);
	data.out.vario.gammah_theo=variogrammodel3d(model,data.out.vario.lag_theo,model_para);
	plotvariogram1d(2)
	para.status.variogram=2;
else
   krig_message(1,['Select ' '''' 'Compute' '''' ' Button to compute Data-based Semi-Variogram/Correlogram from the data']);  
end

