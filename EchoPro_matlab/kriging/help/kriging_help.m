function      krig_help(index,matlab_ver)
%% help file for kriging panel

%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para

varname=[];
content=[];

switch index
  case 1
     varname='Coordinates of Kriging Grids';
     content={
              '  Minimum    --  Mininum value (lower boundary) of the kriging grids',
              '   ',
              '  Maximum    --  Maximum value (upper boundary) of the kriging grids',
              '   ',
              '  Resolution --  Resolution of the kriging grids',
              '   ',
              '    (X,Y,Z)  --  Settings corresponding to coordinates in x, y, and z directions',
              '                 The units are those specified in the ''Data Preparation'' window'};
 case 2
     varname='Kriging Parameters - I';
     content={'Kriging Model:',
              '   ',
			  'Simple Kriging (SK):',
              '        Kriging with zero mean.',
              '     ',
 			  'Ordinary Kriging (OK):',
     		  '        Kriging with a non-zero mean with an additional constraint that the',
       		  '        summation of the weighting coefficents is unity. OK can handle the ',
              '        situation with an unknown constant mean.',
              '    ',
			  'Universal Kriging (Linear drift):',
			  '        Similar to the Ordinary Kriging, but with a linear',
              '        drift mean within the kriging region.'};
  case 3
     varname='Kriging Parameters - II';
     content={'Kriging Scheme:',
         '      ',
        	'   Point to Point Kriging:'
			'         To estimate the value at a specified location.'
            '    ',
			'   Point to Block Kriging:'
         '         To estimate the value of a block centered at a specified location.'
         '         When it is selected, block sizes,  or number of sub-divisions in x,'
         '         y, and z directions need to be provided.'}
 case 5
     varname='Kriging Parameters - III';
     content={'Other Kriging Parameters, including',
        '       '
        '   Search Radius: Search radius for kriging. This parameter specifies a circular',
        '         region within which the observed values are used to estimate',
        '         the value at a given location. However, if the number of observed data',
        '         points is less than the specified minimum kriging points or greater than the',
        '         specified maximum kriging points, it will be modified automatically to',
        '         satisfy the conditions involving the minimum and maximum kriging points.',
        '         This is a normalized dimensionless value, the maximum value is 1.414, or sqrt(2).',
        ' ',
        '   Minimum Kriging Points: The minimum number of points used in kriging to estimate the'
        '         value at a specified location. If within a previously specified search'
        '         radius, the number of data points is less than the specified minimum '
        '         kriging points, the program automatically enlarges the search radius until',
        '         the condition of the minimum kriging point is met.',
        ' ',
 		'   Maximum Kriging Points: The maximum number of points used in kriging to estimate the'
        '           value at a specified location. If within a previously specified search'
        '           radius, the number of data points is greater than the specified maximum '
        '           kriging points, the program automatically reduces the search radius until',
        '           the condition of maximum kriging point is met.',
        ' ',
        '  Relative Variance: Kriging variance normalized by the data variance at the'
        '           location where the kriging estimate is made.'};
  case 6
     varname='Parameter File';
     content={'Load parameter file using file browser -- There are three options:',
           '    ',
        	'   (1) semi-variogram parameters only',
        	'   (2) kriging parameters only',
		 	'   (3) both variogram and kriging parameters.',
           '    ',
           'Save parameter file -- Save both semi-variogram and kriging parameters.'};
case 7
     varname='Load Data Format File';
     content={'Load data format file without loading the data in the ''Data Preparation'' ', ...
         'window. The data format file should be saved previously in the ''Data Preparation'' ', ...
         'window.'};        
  case 8
     varname='Load Data File';
     content={'Load data file without changing the variogram/correlogram parameters.'};
  case 9
     varname='Batch File Processing';
     content={'This option allows user(s) to process data file using the same ',
        'variogram/correlogram, and kriging parameters.',
        '      ',
        '    Load Filename-list File: load the text file in which a list of data files is ',
        '          provided. Each line in the file specifies a file with a full path being included.',
        '   ',
        '    Save Batch Process Log File: specify a filename and location of the log file.'
        '   ',
        '    Batch kriging pushbutton: klick on this button to start the batch kriging.'};
  case 10
     varname='Push Buttons';
     content={
     '  Krig - Start Kriging.',
     '   ',
     '  Refresh - Update data file and variogram parameters, especially when a new data file is loaded.',
     '   ',
     '  Navigator - return to the base window (Navigator).',
     '   ',
     '  Quit - close this window.'};
case 20
     varname='Coordinates of Kriging Grids: Customized (user-defined) Grid File';
     content={'Specify a file that contains the customized grids. There are 2 columns (X,Y) '...
             'in the file for a 2D problem and 3 coluumns (X,Y,Z) for a 3D problem. Note that '...
             'the order of the column vectors (variables) in the file should be the same as the'...
             'order of X, Y, or Z defined in the ''Data Preparation'' window. For example, ', ...
              'for a 2D problem, the first and second columns in the original data file are Latitude '...
              'and Longitude, respectively, but in the ''Data Preparation'' window, the conventional '...
              'way is that the X axis is longitude and the Y axis is latitude. This can be done by '...
              'selecting ''Data Col.2'' for the column of the X-axis and selecting ''Data Col.1'' for '...
              'the column of the Y-axis in the ''Data Preparation'' window. ' ...
              '  '...
              'In the customized grid file, the first and second columns ' ... 
              'should have ''Longitude'' on the first column and ''Latitude'' on the second column. ' ...
              'This way the generic grid file will be independent of the original data file and can '...
              'be applied to situations when the data file formats are different.'...
              '  '...
              'The values on the customized grids are obtained using interpolation of the regular grids'...
              'defined by the ''Coordinates'' in the Kriging window.'};
end
if ~isempty(varname) & ~isempty(content)
  if para.Matlab_Version == 5
     help_message(content,varname,[0 0 0]);
  elseif para.Matlab_Version >= 6
     p=general_message('no_action',varname,content);
  end
end
