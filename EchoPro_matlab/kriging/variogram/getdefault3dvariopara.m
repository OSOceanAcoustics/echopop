function getdefault3dvariopara(opt)
%% function getdefault3dvariopara(opt) sets default variogram parameters
%% opt = 1			start variogram/correlogram window
%%     = 2			refresh the window parameter
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

if para.status.variogram >= 1
   para.vario.max_nugt=max(data.out.vario.gammah);
	para.vario.max_sill=min(1.5,1.5*max(data.out.vario.gammah));
else
	para.vario.max_nugt=1;
   para.vario.max_sill=1.5;
end

para.vario.max_powr=4.0;
para.vario.max_range=sqrt(2);			% normalized range
para.vario.max_hole=4*pi/para.vario.max_range;
para.vario.max_lscl=para.vario.max_range;
para.vario.max_value=3;				% maximum value of semi-variogram

para.vario.range=0.95;
para.vario.res=0.025;								% relative to the total range sqrt(2)
para.vario.vario=1;
para.vario.corr=0;
para.vario.nugt=0;
para.vario.sill=min(1.0,0.7*para.vario.max_sill);
para.vario.lscl=0.1;
para.vario.powr=1.5;
para.vario.hole=0;
para.vario.model=13;
para.vario.disp_range=1.0;
para.vario.disp_value=2.0;

para.vario.azm_beg=0;				% begin azimuth angle in degree
para.vario.azm_end=0;				% end azimuth angle in degree
para.vario.azm_res=200;			% azimuth angle resolution in degree
para.vario.dip_beg=0;				% begin dip angle in degree
para.vario.dip_end=0;				% end dip angle in degree
para.vario.dip_res=100;			% dip angle resolution in degree
para.vario.ytox_ratio=1;				% aspect ratio of y to x
para.vario.ztox_ratio=1;				% angle ratio of z to x

%% other parameters for multiple attributes and additional constraint
para.vario.atol=90;
para.vario.bandwh=1;
para.vario.dtol=90;
para.vario.bandwd=1;
para.vario.nvarg=1;
para.vario.nvar=1;
para.vario.ivtail=1;
para.vario.ivhead=1;
para.vario.ivtype=4;
para.vario.isill=1;
	
set3dvariopara(2);						% set all default parameters including edit field and slider positions
%set3dvariopara(4);						% set anisotropy parameters from para. struct
if para.status.dataprep == 1  & hdl.status.variogramfig == 1	& opt ~= 1 % semi-variogram/correlogram  & 
   														% data-based semi-variogram/correlogram is plotted
   plotvariogram1d(1);							% re-plot data-based semi-variogram/correlogram to
   														% restore semi-variogram
else
   para.vario.corr=0;
end

if opt == 1												% first time to computer semi-variogram/correlogram
   if ~isfield(data,'in') 
      krig_message(1,'Data have not been loaded yet !!!');
   end
else
   if para.status.variogram == 2						% delete previously plotted model-based semi-variogram/correlogram
     delete(hdl.vario.theo_plot);
     variogram_theo(1);
   end
end




