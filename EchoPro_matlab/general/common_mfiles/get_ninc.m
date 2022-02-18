function	[xinc, xdigits, yinc, ydigits]=get_ninc(h, ntick)
% get tick increment for x, y and z axes automatically
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

xlim=get(h,'xlim');
if xlim(2)-xlim(1) > 0.6
   xinc=ceil(ceil((xlim(2)-xlim(1))*15)/10)*10;
   xdigits=0;
elseif xlim(2)-xlim(1) > 0.06
   xinc=ceil((xlim(2)-xlim(1))*15);
   xdigits=0;
elseif xlim(2)-xlim(1) > 0.006
   xinc=ceil((xlim(2)-xlim(1))*150)/10;
   xdigits=1;
end

ylim=get(h,'ylim');
if ylim(2)-ylim(1) > 0.6
  yinc=ceil(ceil((ylim(2)-ylim(1))*15)/10)*10;
  ydigits=0;
elseif ylim(2)-ylim(1) > 0.06
   yinc=ceil((ylim(2)-ylim(1))*15);
   ydigits=0;
elseif ylim(2)-ylim(1) > 0.006
   yinc=ceil((ylim(2)-ylim(1))*150)/10;
   ydigits=1;
end
