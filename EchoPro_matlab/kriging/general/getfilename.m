% function fname_out=getfilename(str_in,para) converts window-based
% filename to Linux/Unix-based filename
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

function fname_out=getfilename(str_in,para)

fname_out=str_in;
if para.platform == 1
    return
else
    slash_str='/';
    indx=find(str_in == '\');
    fname_out(indx)=slash_str(ones(size(indx)));
end