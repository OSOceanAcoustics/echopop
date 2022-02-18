function		[x,y,z,v]=remove_nan(filename)
%%% remove NaN's
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl
dat=load(filename);
if (get(hdl.dataprep.var1,'Value') == 2 & get(hdl.dataprep.var2,'Value') == 1) | ...
	(get(hdl.dataprep.var1,'Value') == 5 & get(hdl.dataprep.var2,'Value') == 4)
	para.dataprep.xy_switch=1;
	x=dat(:,2);
	y=dat(:,1);
else
	para.dataprep.xy_switch=0;
	x=dat(:,1);
	y=dat(:,2);
end
if size(dat,2) >= 4
   z=dat(:,3);
   v=dat(:,4);
	indx=find( isnan(x) | isnan(y) | isnan(z) | isnan(v));
	x(indx)=[];
	y(indx)=[];
	z(indx)=[];
   v(indx)=[];
else
   v=dat(:,3);
	indx=find( isnan(x) | isnan(y) | isnan(v));
	x(indx)=[];
	y(indx)=[];
   v(indx)=[];
   z=[];
end