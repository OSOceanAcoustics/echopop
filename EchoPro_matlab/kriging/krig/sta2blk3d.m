function Csb=sta2blk(x1,y1,z1,x2,y2,z2,nx,ny,nz,Lx,Ly,Lz,model,model_para)
%% stations to a square block 
%% (x1,y1)  coordinates of stations
%% (x2,y2)  coordinates of the center of block 
%% nx x ny   no. of elements of the block
%% Lx, Ly     block size  
%%  model:  semi-variogram model index
%%  model_para model parameters
%
% Author: Jim Ledwell, 10/97 Woods Hole Oceanographic Institution
% Revised by Dezhang Chu,   10-29-98
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

ns=length(x1);
blk_x_res=Lx/nx;				% horizontal resolution within the block
blk_y_res=Ly/ny;				% horizontal resolution within the block
blk_z_res=Lz/nz;				% vertical resolution within the block
%% create block matrices 
xl2=linspace(x2-Lx/2+blk_x_res/2,x2+Lx/2-blk_x_res/2,nx);
yl2=linspace(y2-Ly/2+blk_y_res/2,y2+Ly/2-blk_y_res/2,ny);
zl2=linspace(z2-Lz/2+blk_z_res/2,z2+Lz/2-blk_z_res/2,nz);

[X2,Y2,Z2]=meshgrid(xl2,yl2,zl2);
CVV=0;
for i=1:ns
   dX=x1(i)-X2;
   dY=y1(i)-Y2;
   dZ=z1(i)-Z2;
   R=sqrt(dX.*dX + dY.*dY + dZ.*dZ);
   A=variogrammodel3D(model,R,model_para);
   Csb(i)=sum(sum(sum(A)))/(nx*ny*nz);						% simple summation
end
