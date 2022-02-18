function  para=auto_load_vario_krig_para(filename)

global para	data 		% all setting and processing parameters

[dat, name]=xlsread(filename);
% get variogram & kriging parameters
for i=1:length(dat)
   cmd=char({['para.' name{i} '=' sprintf('%g',dat(i)) ';']});
%    fprintf('%d\t %s\n',i,cmd)
   eval(cmd);
end

%    para.dataprep.y_offset=45;
%    para.dataprep.x_offset= -124.1;
%    para.dataprep.y_offset=46.4;
%     para.dataprep.x_offset= -124;

%     para.krig.xmin0=xmin;
%     para.krig.xmax0=xmax;
%     para.krig.dx0=dx;
%     para.krig.ymin0=ymin;
%     para.krig.ymax0=ymax;
%     para.krig.dy0=dy;
% dat=data.in.US_CAN_dat;
% data.in.var1_raw=dat(:,1);
% data.in.var2_raw=dat(:,2);
% data.in.var3_raw=[];
% data.in.dim=2;
% data.in.var_raw=dat(:,3);

para.dataprep.x_offset=para.krig.shelf_break_Ref_lon;
y=data.in.US_CAN_dat(:,1);
x=data.in.US_CAN_dat(:,2);
para.dataprep.latlonfac_min=cos(para.krig.ymax0*pi/180);     % corresponding to minimum longitude
para.dataprep.latlonfac_max=cos(para.krig.ymin0*pi/180);     % corresponding to maximum longitude
para.dataprep.x_norm=max(x)-min(x);
para.dataprep.y_norm=max(y)-min(y);
para.krig.xmin=min((para.krig.xmin0-para.dataprep.x_offset)*para.dataprep.latlonfac_min/para.dataprep.x_norm);
para.krig.xmax=max((para.krig.xmax0-para.dataprep.x_offset)*para.dataprep.latlonfac_max/para.dataprep.x_norm);
para.krig.ymin=(para.krig.ymin0-para.dataprep.y_offset)/para.dataprep.y_norm;
para.krig.ymax=(para.krig.ymax0-para.dataprep.y_offset)/para.dataprep.y_norm;
para.krig.dx=para.krig.dx0/para.dataprep.x_norm;
para.krig.dy=para.krig.dy0/para.dataprep.y_norm;

return