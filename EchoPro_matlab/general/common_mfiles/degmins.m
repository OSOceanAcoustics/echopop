function degstrings=degmins(degrees,ndigit);
% function degstrings=degmins(degrees,ndigit);
% creates a degrees and minutes label for use in MAPAX routine.
%
% Usage:  degstring=degmins(degrees,ndigit);
%
%    Inputs:  degrees = decimal degrees
%             ndigit  = number of decimal places for minutes
%
%    Outputs: degstring = string containing label
% Version 1.0 Rocky Geyer (rgeyer@whoi.edu)
% Version 1.1 J.List (jlist@usgs.gov)
%             fixed bug
%%
%% Modified by Dezhang Chu, Sept. 31, 1999

degrees=degrees(:);

for i=1:length(degrees);

	if degrees(i)<= 0 
  		minus_sign=1;
  		degrees(i)=-degrees(i);
   else
      minus_sign=0;
	end

	deg(i,1)=floor(degrees(i));
	deg(i,2)=round((degrees(i)-deg(i,1))*600)/10;
   
	if ndigit==0;
 	   if deg(i,1) >= 100 
    		tmp_str=sprintf('%3d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'''');
      elseif deg(i,1) >= 10 
    		tmp_str=sprintf('%2d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'''');
      else
    		tmp_str=sprintf('%1d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'''');
	   end
	elseif ndigit==1;
  		deg(i,3)=round(10*abs(floor(deg(i,2))-deg(i,2)));
  	   if deg(i,1) >= 100 
 			tmp_str=sprintf('%3d%s%2.2d%s%1.1d%s',deg(i,1),...
                 setstr(176),floor(deg(i,2)),'.',deg(i,3),'''');
      elseif deg(i,1) >= 10 
  			tmp_str=sprintf('%2d%s%2.2d%s%1.1d%s',deg(i,1),...
                 setstr(176),floor(deg(i,2)),'.',deg(i,3),'''');
      else 
  			tmp_str=sprintf('%1d%s%2.2d%s%1.1d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'.',deg(i,3),'''');
      end
	elseif ndigit==2;
  		deg(i,3)=round(100*abs(round(deg(i,2))-deg(i,2)));
  	   if deg(i,1) >= 100 
  			tmp_str=sprintf('%3d%s%2.2d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'.',deg(i,3),'''');
  	   elseif deg(i,1) >= 10 
  			tmp_str=sprintf('%2d%s%2.2d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'.',deg(i,3),'''');
  	   else
  			tmp_str=sprintf('%1d%s%2.2d%s%2.2d%s',deg(i,1),...
                 setstr(176),round(deg(i,2)),'.',deg(i,3),'''');
      end
	end;
   degstring(i,1:length(tmp_str))=tmp_str;
	if minus_sign == 0
 		tmp=['-' degstring(i,:)];
 	 	degstrings(i,1:length(tmp))=tmp;
   else
 	 	degstrings(i,1:length(degstring(i,:)))=degstring(i,:);
	end
end
