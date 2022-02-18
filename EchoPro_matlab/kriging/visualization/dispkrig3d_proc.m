function		dispkrig3d_proc
%%  determine whether is a 2-D kriging or a 3-D kriging
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

if data.in.dim == 3
   dispkrig3d(0);
else
   dispkrig3d(5);
end
