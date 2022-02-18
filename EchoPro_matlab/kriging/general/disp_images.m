function		disp_images
% display globec (Georges bank), whoi logos and a kriging example
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


axes('position',[0.625 0.5 0.25 0.3])
[imgdat,cc1]=imread('whoi_logo1.tif','tif');
himg1=image(imgdat);
colormap(cc1)
cc1=colormap;
axis off

%% EXAMPLE IMAGES
axes('position',[0.375 0.2 0.25 0.3])
[imgdat,cc2]=imread('sample.tif','tif');
himg=image(imgdat);
axis([180 1080 90 800])
colormap(cc2);
axis off


axes('position',[0.125 0.5 0.25 0.3])
[imgdat,cc3]=imread('globec1.tif','tif');
himg2=image(imgdat);
colormap(cc3)
axis off
