function      variogram_help(index)
%% function      variogram_help(index)
%% displays help file for variogram window
%%
%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

varname=[];
content=[];

switch index
  case 1
     varname='Semi-Variogram/Correlogram';
     content={
        'Variogram  --  The semi-variogram normalized by the autocorrelation ',
        '                function at lag = 0 (variance)',
        'Correlogram  -- The covariance function normalized by its ',
        '                value at lag = 0 (variance).', 
        '  ',
        'The relation between the normalized semi-variogram and the correlogram is:',
        '   ',
        '      Vg(h) = 1 - Cr(0)'};
  case 15
     varname='Nugget';
     content={'A quantity due to nugget effect. ',
        'For a semi-variogram, it is the value ' ,
        'at lag = 0, while for a correlogram it is one minus the ' ,
        'value at lag = 0'};
  case 16
     varname='Sill';
     content={'The asymptotic value as lag approaches infinity for ' ,
        'a semi-variogram. In most cases, this value is close to 1.'};
  case 17 
     varname='Length Scale';
     content={'It is a normalized length scale that reflects the ascending rate for a semi-variogram ' ,
        'or the descending rate for a correlogram. The maximum normalized length scale is sqrt(2), i.e. ' ,
        'both horizontal and vertical coordinates are normalized to unity'};
  case 18
     varname='Power';
     content={'This is a parameter used for general exponential-bessel model ' ,
        'only, and its value controls the ascending rate for a semi-variogram or ' ,
        'descending rate for a correlogram. Particularly, a power of 1 gives ' ,
        'the exponential-bessel model while a power of 2 gives the Gaussian-bessel model'};
  case 19
     varname='Scale of Hole Effect';
     content={'It is a normalized length scale resulting from a hole structure ' ,
        'in the data, i.e., a region whose overall values are very small ' ,
        'compared with those of the adjacent regions. Such a feature will ' ,
        'generate an oscillatory pattern in either a semi-variogram ' ,
        'or a correlogram.'};
  case 20
     varname='Range';
     content={'The maximum range for computing the semi-variogram/correlogram. ' ,
	     'It is normalized by the maximum (diagonal) distance of the rectangular data ',
         'observation region.'};
  case 21
     varname='Resolution:';
     content={'The lag resolution used in computing the semi-variogram ' ,
        'using either measured data or the theoretical variogram model. ' ,
        'It is a normalized quantity.'};
  case 22
     varname='Load/Save Parameter file';
     content={'Load and Save a variogram parameter file using file browser.' }
  case 25
	  varname='Isotropy/Anisotropy:';
     content={'Enable 2-D/3-D semi-variogram/correlogram computations. ',
			'The 2-D semi-variogram/correlogram map provides anisotropy information ' ,
			'which helps user(s) to choose (angle, aspect ratio) of the anisotropic ' ,
  			'feature for 2-D variogram, and (azimuth, dip, two aspect ratios) of the anisotropic' ,
            'feature for 3-D variogram.',
            '    ',
            'Note that once ''Isotropy'' is enabled, the program will compute an EQUIVALENT',
            '1-D variogram/correlogram using the parameters: ''Rotation (deg)'' and',
            '''Ratio'' including "Y/x" and "Z/x" (3-D) set once ''Anisotropy'' is enabled.',
            'The program transfers the original anisotropic data to isotropic data using these',
            'parameters.'};
  case 26
	 varname='Anisotropic Semi-variogram/Correlogram - I';
     content={'Principal Data Plane: ',
               '      ',
			   'Azimuth  --  Data are processed in the azimuthal plane',
               '      ',
               '    Dip  --  Data are processed in the vertical plane',
               '      ',
               'For the 2-D case, one can only choose ''Azimuth'' as the computational domain,',
               'while for the 3-D case, parameters for both will need to be defined.'};
  case 27
     varname='Anisotropic Semi-variogram/Correlogram - II';
     content={'Begin Angle  --  0 to 180 for azimuthal angle, -90 to 90 for dip angle',
               '      ',
              '  End Angle  --  0 to 180 for azimuthal angle, -90 to 90 for dip angle',
               '      ',
              'Angle Increment  --  Angle increment in computing semi-variogram/Correlogram'};
  case 28
     varname='Anisotropic Semi-variogram/Correlogram - III';
     content={'Tolerance  --  This parameter specifies the range where the maximum width of the ',
             'semi-variogram/Correlogram computation region in a particular angle is defined as ',
             'the projected arc-length at this range (see Figure III.2 on p49 in GSLIB geostatistical',
             'Software Library and User''s Guide by Deutsch and Journel, 1998, Oxford' ,
             ' Univ. Press, New York). Note that the range is the normalized range (see Help ''Range'').'};
  case 29
     varname='Anisotropic Semi-variogram/Correlogram - IV';
	  content={'Rotation Angle (degree):',
               '      ',
               'Azimuth  --  Angle rotated in the azimuthal (horizontal) plane',
               '      '
               '    Dip  --  Angle rotated in the vertical plane',
               '      ',
               'These two angles should correspond to the pricipal axes of the',
               'Semi-variogram/Correlogram. Semi-variogram/correlogram and ',
		       'kriging are performed using transformed data, but the final kriging ',
		       'results are still displayed in the original (unrotated) coordinates, ',
               'i.e., user will not see the intermediate results in their isotropic form. ',
               'The angle (in degree) is measured counterclockwise in azimuth but clockwise',
               'in the vertical plane (dip angle)'};
  case 30
     varname='Anisotropic Semi-variogram/Correlogram -V';
     content={'Aspect Ratio of the Anisotropic Semi-variogram/Correlogram: ',
              '   ',
              'This parameter specifies a ratio for transforming the raw anisotropic',
              'data to isotropic data. Semi-variogram/Correlogram and kriging are ',
		      'performed using transformed data, but the final kriging results are',
		      'still displayed in the original (untransformed) coordinates, ',
              '        '
		      'Y/X  --  This parameter is defined as the ratio of the latitudinal scale',
              '             to the longitudinal scale in the azimuthal (horizontal)',
              '             plane in the rotated coordinates.',
              '        ',
		      'Z/X  --  This parameter is defined as the ratio of the vertical scale',
              '             to the horizontal scale in the vertical plane in the',
              '             rotated coordinates.',
               };
  case 40
     varname='Display Range and Value';
     content={'Lag: displayed normalized range (see help on ''Range'');',
              'Value: displayed vertical range. ' };  
  case {50,51,52,53}
     varname='Push Buttons ';
     content={ 'Compute - Compute either the semi-variogram or correlogram using ' ,
               '          observed (measured) data.' ,
               '     ',
               'LSQ Fit - Perform a Least-Square Fit to find a set of parameters ' ,
               '          that characterize either the semi-variogram or the correlogram.' ,
               '     ',
               'Default - Restore the default parameters for computing the ' ,
               '          semi-variogram or correlogram.',
               '     ',
               '  Apply - Compute the theoretical semi-variogram/correlogram using a ' ,
               '          chosen model.'};
  case 55
     varname='Other Push Buttons';
     content={' Navigator -- return to the base window (Navigator). ',
              '   ',
              '  Quit  -- close this window.'};
 case 100
     varname='Advanced';
     content={'Advanced parameter setting for least-square-fit.'};
end
if index < 2 | index > 14
  if ~isempty(varname) & ~isempty(content)
    if para.Matlab_Version == 5
       help_message(content,varname,[0 0 0]);
    elseif para.Matlab_Version >= 6
       p=general_message('no_action',varname,content);
    end
  end
else
   model_index=index-1;
   modelstring(model_index);
end




