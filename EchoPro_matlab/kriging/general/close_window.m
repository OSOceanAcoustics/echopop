function			close_window(window_index0)
% funstion			close_window(window_index0)
% close window based on the window index
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl 

for i=1:length(window_index0)
    window_index=window_index0(i);
    switch window_index
        case 0
            if isfield(hdl,'navigator') & ishandle(hdl.navigator.h0)
                close(hdl.navigator.h0);
            end
        case 1			% data preparation window
            if ishandle(hdl.dataprep.h0)
                close(hdl.dataprep.h0);
                hdl.status.dataprepfig=0;
            end
        case 2			% semi-variogram/correlogram window
            if ishandle(hdl.vario.h0)
                close(hdl.vario.h0);
                hdl.status.variogramfig=0;
            end
        case 3			% kriging	window
            if ishandle(hdl.krig.h0)
                close(hdl.krig.h0);
                hdl.status.krigingfig=0;
            end
        case 4		% visulization window
            if isfield(hdl,'dispkrig3d') & ishandle(hdl.dispkrig3d.h0)
                close(hdl.dispkrig3d.h0);
                hdl.status.dispkrigfig=0;
            end
            if hdl.status.krigingfig == 1
                p=findobj(hdl.krig.h0,'type','axes');
                if ~isempty(p) delete(p);end
                return
            end
        case 5		% 2D-3D variogram/correlogram visulization window
            if isfield(hdl,'dispvario2d3d') & ishandle(hdl.dispvario2d3d.h0)
                close(hdl.dispvario2d3d.h0);
                hdl.status.dispvario2d3d=0;
            end
    end
end

return