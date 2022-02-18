function CVV=blk2blk(x1,y1,z1,x2,y2,z2,nx,ny,nz,Lx,Ly,Lz,model,model_para)
%% a square block to a square block 
%% (x1,y1,z1)  coordinates of the center of block 1
%% (x2,y2,z2)  coordinates of the center of block 2
%% nx x ny x nz   no. of elements of the block
%%  Lx,Ly,Lz     block size  
%%  model:  semi-variogram model index
%%  model_para model parameters
% Author: Jim Ledwell, 10/97 Woods Hole Oceanographic Institution
% Revised by Dezhang Chu,   10-29-99
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

blk_x_res=Lx/nx;				% horizontal resolution within the block
blk_y_res=Ly/ny;				% horizontal resolution within the block
blk_z_res=Lz/nz;				% horizontal resolution within the block
%% create block matrices 
xl1=linspace(x1-Lx/2+blk_x_res/2,x1+Lx/2-blk_x_res/2,nx);
yl1=linspace(y1-Ly/2+blk_y_res/2,y1+Ly/2-blk_y_res/2,ny);
zl1=linspace(z1-Lz/2+blk_z_res/2,z1+Lz/2-blk_z_res/2,nz);
xl2=linspace(x2-Lx/2+blk_x_res/2,x2+Lx/2-blk_x_res/2,nx);
yl2=linspace(y2-Ly/2+blk_y_res/2,y2+Ly/2-blk_y_res/2,ny);
zl2=linspace(z2-Lz/2+blk_z_res/2,z2+Lz/2-blk_z_res/2,nz);

[X2,Y2,Z2]=meshgrid(xl2,yl2,zl2);
CVV=0;
for i=1:nx
   dX=xl1(i)-X2;
   for j=1:ny
     dY=yl1(j)-Y2;
     for k=1:nz
       dZ=zl1(j)-Z2;
       R=sqrt(dX.*dX + dY.*dY +dZ.*dZ);
       A=variogrammodel3D(model,R,model_para);
       CVV=CVV+sum(sum(sum(A)));						% simple summation
     end
   end
end
CVV=CVV/(nx*nx*ny*ny*nz*nz);
