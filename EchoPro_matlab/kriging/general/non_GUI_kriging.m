function non_GUI_kriging

global para	 		% all setting and processing parameters
global data  		% input and output data

%% get grid cell
dat=data.in.US_CAN_mesh_cor;
data.out.krig.gx=dat(:,2);
data.out.krig.gy=dat(:,1);  

%model_para=[para.vario.nugt para.vario.sill para.vario.lscl para.vario.powr para.vario.hole];
%krig3dmanager;
non_GUI_krig3dmanager;
return