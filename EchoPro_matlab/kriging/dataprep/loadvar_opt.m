function	[x,y,z,v]=loadvar_opt(filename);
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


dat=load(filename);
[n,m]=size(dat);
x=dat(:,1);
y=dat(:,2);
if m >= 4
  z=dat(:,3);
  v=dat(:,4);
else
  v=zeros(n,1);
  z=dat(:,3);
end

indx=find(isnan(x)|isnan(y)|isnan(z)|isnan(v));
x(indx)=[];
y(indx)=[];
z(indx)=[];
v(indx)=[];

if m <= 4
   v=[];
end
