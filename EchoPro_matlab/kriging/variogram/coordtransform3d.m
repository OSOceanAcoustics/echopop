function [cx_out]=coordtransform3d(cx_in)
%  function [cx_out]=coordtransform3d(cx_in) takes cx_in as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in transform_para (defined by the structure para.vario)
% Rotations are all performed anticlockwise with the observer located on the positive side of 
% the axis and looking toward the origin. In 3D, rotations are performed first along z,
% then along rotated y and then along twice rotated x.
% Author: D. Marcotte (Version 2.1  97/aug/1)  
%%
%% Modified by Dezhang Chu,  April 3, 2003. some constants are defined
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para data

[n,d]=size(cx_in);


% perform rotation counterclockwise

   if d==2
      ang=para.vario.azm_rot; cang=cos(ang/180*pi); sang=sin(ang/180*pi);
      rot=[cang,-sang;sang,cang];
      t=[1 para.vario.ytox_ratio];
   else
      % rotation matrix in 3-D is computed around z, y and x in that order
      angz=para.vario.azm_rot; cangz=cos(angz/180*pi); sangz=sin(angz/180*pi);
      angy=para.vario.dip_rot; cangy=cos(angy/180*pi); sangy=sin(angy/180*pi);
      angx=0; cangx=cos(angx/180*pi); sangx=sin(angx/180*pi);
      rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
      roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
      rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
      rot=rotz*roty*rotx;
  	   t=[1 para.vario.ytox_ratio para.vario.ztox_ratio];
   end
% rotation is performed around z, y and x in that order, the other
% coordinates are left unchanged.
   cx_out=cx_in*rot;
   t=diag(t);

% perform contractions or dilatations (reduced h)
cx_out=cx_out/t;
para.vario.rotmat=rot;
para.vario.ampmat=t;
