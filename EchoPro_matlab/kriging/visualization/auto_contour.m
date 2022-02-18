function			contour_lines=auto_contour(var,n,digits)
%% function	contour_lines=auto_contour(var,n,digits)
%% automatically determine the values of contour lines
%% INPUT:
%%   var  =  variables to be contoured
%%     n  =  number of contour lines with values  
%%        at ratios of i/(n+1) of the data range, where i = 1, ...,n
%%        i.e.,  vi= i/(n+1)*(vmax-vmin)+vmin
%% digits =  	non-zero digits 
%% OUTPUT:
%%  contour_line - values of the automatically determined contour lines
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

i=1:n;
tmp=(i/(n+1)*(max(var)-min(var))+min(var));
fac=10.^(-(floor(log10(tmp))-(digits-1)));
contour_lines=round(tmp.*fac)./fac;
