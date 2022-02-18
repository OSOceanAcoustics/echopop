function non_GUI_loaddatfile

global para data

dat=data.in.US_CAN_dat;
data.in.var1_raw=dat(:,1);
data.in.var2_raw=dat(:,2);
data.in.var3_raw=[];
data.in.dim=2;
data.in.var_raw=dat(:,3);

data.in.var1=dat(:,2);
data.in.var2=dat(:,1);
data.in.var3=dat(:,3);
data.in.var=data.in.var_raw;
para.dataprep.x_axis_indx=2;
para.dataprep.y_axis_indx=1;
para.dataprep.checkunit_action=0;



