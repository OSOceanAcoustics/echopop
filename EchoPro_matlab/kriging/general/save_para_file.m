function		save_para_file(filename,type)
% Save parameters from a file
%% type = 1   		parameters
%%		= 2			parameters and data
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data



switch type
	case 1
		cmd=['save ''' filename '''  para'];
	case 2
		cmd=['save ''' filename '''  para   data'];
end
if para.Matlab_Version == 7
%   cmd=[cmd ' -nounicode'];
end
eval(cmd)

return