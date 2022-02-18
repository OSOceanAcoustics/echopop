function non_GUI_load_data

global para	 		% all setting and processing parameters
global data  		% input and output data

non_GUI_loaddatfile;				% load data into memory
krig_datachk(2);		    % delete nans, filtering and data reduction

%% Data Transformation
data.in.tvar=datatransform(1,data.in.var,1);		% forward transformation on the original data
data.in.tv=datatransform(1,data.in.v,1);			% forward transformation on the normalized data
para.vario.c0=std(data.in.tv)^2;
data.out.vario.c0=para.vario.c0;

return