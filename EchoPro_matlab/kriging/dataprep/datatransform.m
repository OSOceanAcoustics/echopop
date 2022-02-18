function    yout=datatransform(trans_flag,yin,trans_type)
%%% Data Transformation/Inverse Data Transformation
%%			trans_flag = 1  Forward Transform
%%                     = 2  Inverse Transform
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global data para hdl

if nargin  == 2
  trans_type=get(hdl.dataprep.transform,'value');
end
para.dataprep.transform_index=trans_type;
if trans_flag == 1  % Data Transformation with specified type
   switch trans_type
	  case 1
		  yout=yin;
     case 2
        yout=log10(yin+1);
     case 3
        yout=log(yin+1);
     case 4
        yout=10*log10(yin+1);
     case 5
        yout=10*log(yin+1);
     case 6
        yout=log10(abs(yin)+eps);
     case 7
        yout=log(abs(yin)+eps);
    end
  else					%% Inverse Data Transformation
    switch trans_type
 	  case 1
		  yout=yin;
     case 2
        yout=10.^yin-1;
     case 3
        yout=exp(yin)-1;
     case 4
        yout=10.^(0.1*yin)-1;
     case 5
        yout=exp(0.1*yin)-1;
     case 6
        yout=10.^yin;
     case 7
        yout=exp(yin);
    end
  end