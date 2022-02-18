function      navigator_help(index)
%% function      navigator_help(index)
%% displays help file for navigator window
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

window_pos=[];

switch index
  case 1
     varname='Task - I';
     content={
        '   Load Data  --  Load the 2-D/3-D data file. ',
        '   ',
        '   Variogram  --  Compute semi-variogram/correlogram.',
        '   ',
        '            Krig  --  Perform kriging or batch kriging. ',
        '    ',
        'Visualization --  Display and save the kriged results. ', 
              };
  case 2
     varname='Task - II';
     content={
        'Save Window Position  --  Save the current window settings as the default settings,',
        '                             including window size and position. This settings will affect all',
        '                             of other task windows.',
        '       ',
        '    Quit    -- Quit Easykrig.'};
     
end
%p=general_message('no_action',varname,content);

if para.Matlab_Version == 5
     help_message(content,varname,[0 0 0]);
elseif para.Matlab_Version >= 6
     p=general_message('no_action',varname,content);
end




