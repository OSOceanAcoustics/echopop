function Csb=sta2blk2d(x1,y1,x2,y2,nx,ny,Lx,Ly,model,model_para)
%% function Csb=sta2blk(x1,y1,x2,y2,nx,ny,Lx,Ly,model,model_para)
%% computes covariance of stations to a square block 
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
%% create block matrices 
xl2=linspace(x2-Lx/2+blk_x_res/2,x2+Lx/2-blk_x_res/2,nx);
yl2=linspace(y2-Ly/2+blk_y_res/2,y2+Ly/2-blk_y_res/2,ny);

[X2,Y2]=meshgrid(xl2,yl2);
CVV=0;
for i=1:ns
   dX=x1(i)-X2;
   dY=y1(i)-Y2;
   R=sqrt(dX.*dX + dY.*dY);
   A=variogrammodel(model,R,model_para);
   Csb(i)=sum(sum(A))/(nx*ny);						% simple summation
end
