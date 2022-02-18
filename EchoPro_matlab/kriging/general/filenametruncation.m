function  stringname=filenametruncation(filename)
% function  stringname=filenametruncation(filename)
%% remove path from the full path-filename string and the file extensions 
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


pp=find(filename == '/' | filename == '\' | filename == ':');
if length(pp) > 0
   indx1=pp(length(pp))+1;
else
   indx1=1;
end
qq=find(filename == '.');
if length(qq) > 0 
  indx2=qq(length(qq))-1;
else
  indx2=length(filename);
end
stringname=filename(indx1:indx2);
