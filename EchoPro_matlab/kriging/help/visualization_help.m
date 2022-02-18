function      displaykrighelp(index)
%% function      displaykrighelp(index)
%% displays help file for display-kriging-result window
%%  index = index of options
%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

varname=[];
content=[];
switch index
  case 0						% not currently available
     varname='Interp to Mesh Node and Save To S2R File';
     content={'Interpolate the values for kriging image (regular gridded data ' ...
        'to the mesh nodes specified in a chosen file and save the interpolated ' ...
        'values to a file with S2R file format.'};
  case 1					
     varname='(X,Y,Z) Coordinate Control';
     content={'Reverse  --  reverse the display direction of corresponding axis. ', 
              '    ',
              ' Slider  --  choose the slice position in the corresponding coordinate ', 
              '             (for 3D case only).',
              '    ',
		      '  Value  --  value of the corresponding coordinate of the slice.',
              '             (for 3D case only).'};
  case 2
	  varname='Colorbar Control';
	  content={'Change the sliders will change the upper and lower limits of the colormap. '};
  case 3
	  varname='Displayed Variable';
	  content={'Select what variable to display: ',
               '    ',
			   '         Krig Map  --  Display the kriged values.' ,
               '    ',
               '         Variance Map   --  Display the krig variance normalized by the variance,',
               '                       C(0), from the covariance function.',
               '    ',
               '   Validation  --  Display the cross-validation results.'};
  case 5
	  varname='Cross-Validation: Q-1 & Q-2';
	  content={'Q-1 and Q-2 cross validations are used to check the statistical distribution ',
               'of the residuals between the observed data and kriged values at the original ',
               'observation locations by using the same kriging parameters and variogram model ',
               'parameters. To perform Q-1 and Q-2 cross validations, a normalized residual ',
               'array (Er) needs to be constructed (see ''Introduction to geostatistics, ',
               'Application in Hydrogeology '' by Kitanidis, Cambridge Univ. Press, 1997, ',
              'pp 86-95).',
              'Q-1  --   checks the statistics of the mean of the residual Er and approximately ',
              '          follows the normal distribution.',
              'Q-2  --   checks the statistics of the variance of Er. (Q-2)*(n-1) approximately',
              '          follows the chi-square distribution with parameter n-1. ',
              '   ',
              'The acceptable region defined in the program (two black vertical lines) is the ',
              '0.025 and 0.975 percentiles.'};
  case 7
	  varname='Cross-Validation: Double Kriging';
	  content={'This cross-validation scheme is to evaluate the level of agreement between the'...
              'kriged or predicted values and the original observations at all observation '...
              'locations. The predicted data at grids obtained from the kirging (first pass) '...
              'are served as ''input data''. The mean value at the original observation locations are '...
              'estimated by kriging (second pass) with the same kriging parameters and variogram model '...
              'parameters. The results from the second kriging are then compared with the original '...
              'observed data in a separate plot.'};
  case 8
      varname='Cross-Validation: Leave One Out';
	  content={'Same as '' Double Kriging'' except that, for each location, all observed data '...
              'but the one at that location are served as ''input data'' in performing second '...
              'pass kriging. The results from this kriging are compared with the original '...
              'observed data in a separate plot.'};    
  case 9
      varname='Cross-Validation: Compute/Re-Compute';
	  content={'  Compute  --  Perform the cross-validation computation and plot the results.',
              '    ',
               'Re-Compute --  Using the new kriging/variogram parameters to perform the new ',
               '                computation.'};
  case 10
      varname='Shading';
	  content={'This option allows the user to select how the kriging map is shaded.'...
               'Detailed explanation of this option can be found by typing ''help shading''  in '...
              'the Matlab command window.'};
  case 11
      varname='Push Button';
	  content={'  Show Plot  --  Show kriged/variance map ',
               '    ',
               '       Load  --  Load kriging results from previously saved file (can also',
               '                         be loaded from the File option from the Manu-bar).',
               '    ',
               '   Navigator --  Return to the base window (Navigator).',
               '    ',
               '       Quit  --  Close the current window.'};
  case 20
      varname='Selection ''Data'' on the Menu-Bar';
	  content={'This option will give you flexibility of showing the raw data and',
               'display the kriged results on the customized grids.',
               '  ',
               '  Colored-Coded --  The colors of the data are coded based on the data values ',
               '                          defined by the colormap.',
               '    ',
               '   Black/White  --  The colors of the data are either black or white.'};
  case 21
      varname='Selection ''File'' on the Menu-Bar';
	  content={'This option will allow user to load and save data/kriging results.',
               '  ',
               '  Load --  Load kriging results from previously saved file (the data can ',
               '             also be loaded using the ''Load'' pushbutton.',
               '    ',
               '   Save Figure As  --  Save the entire window to a Matlab fig file.',
               '    ',
               '   Save File As  --  Save the raw data and the kriging results to a *.mat file, along',
               '                      with all necessary parameters including variogram and kriging ',
               '                      parameters. The file is a binary file in Matlab format.',
               '    ',
               '   Export ...  --  Save the entire window, or figure in a different format format, such as',
               '                    eps, tiff, jpeg, etc.'};
  case 22
      varname='Selection ''Data'' on the Tool-Bar';
	  content={'This option will give you flexibility of showing the raw data and',
               'display the kriged results on the customized grids.',
               '  ',
               '   Colored-Coded --  The colors of the data are coded based on the data values ',
               '                          defined by the colormap.',
               '    ',
               '    Black/White  --  The colors of the data are either black or white.',
               '  ',
               '          Value  --  Display the numerical numbers of the data.',
               '    ',
               '          Error  --  Display the numerical numbers of the difference between',
               '                          the data and the kriged results at the data positions.',
               '    ',
               '          None   --  Raw data are not displayed on the map',
               '    ',
               '          Text Color and the Text Fontsize are the face color of the text and the fontsize of the',
               '          selected variables, either ''Value'' or ''Error''. '}; 
  case 23
      varname='Selection ''Data'' on the Tool-Bar';
	  content={'None -- Do not display the data and/or the kriged results ',
               '        on the customized grids.'};
  case 24
      varname='Selection ''Color Map'' on the Tool-Bar';
	  content={'Color Map option will allow user to select the appropriate color map.',
               '        '};
  end
  if ~isempty(varname) & ~isempty(content)
   if para.Matlab_Version == 5
     help_message(content,varname,[0 0 0]);
   elseif para.Matlab_Version >= 6
    p=general_message('no_action',varname,content);
   end
  end


