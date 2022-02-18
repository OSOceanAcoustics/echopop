function check_unit_action(h0)
% convert strings to numbers
%
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global data hdl para

para.dataprep.ytox=str2num(get(hdl.unit_conv.ytox,'string'));
para.dataprep.ztox=str2num(get(hdl.unit_conv.ztox,'string'));
close;
dataprep3d(2)
