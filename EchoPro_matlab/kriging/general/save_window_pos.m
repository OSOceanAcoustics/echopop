function	save_window_pos(target_window_h0)
% save window position including size and position information
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

%win_pos=get(hdl.navigator.h0,'position');
win_pos=get(target_window_h0,'position');

if para.Matlab_Version == 7
   save window_position win_pos  -nounicode
else
   save window_position win_pos  
end
hdl.window_position=win_pos;
