function []=mapax(nminlong,ndiglong,nminlat,ndiglat,hdl,axis_opt);
% function []=mapax(nminlong,ndiglong,nminlat,ndiglat,hdl,axis_opt);
%   puts degrees and minutes on map axes instead of decimal degrees
%
%   Usage: mapax(nnminlong,ndiglong,nminlat,ndiglat);
%
%   Inputs:
%          nminlon  = minutes of spacing between longitude labels
%          ndiglong = number of decimal places for longitude minute label
%
%          nminlon  = minutes of spacing between latitude labels
%          ndiglong = number of decimal places for latitude minute label
%
% Example:  mapax(15,1,20,0);   
%               labels lon every 15 minutes with 1 decimal place (eg 70 40.1')
%           and labels lat every 20 minutes with no decimal place (eg 42 20')
%
% Version 1.0 Rocky Geyer  (rgeyer@whoi.edu)
% Version 1.1 J. List (6/5/95) had apparent bug with
%  ndigit being set to 0: routine degmins blows up.
%  Fixed by adding arguments specifying number of decimal
%  digits (can vary from 0 to 2)
%  Modified by D. Chu 7-16-99, add optional graphic handle 
%  Modified by D. Chu 9-3-99,  add axis label option
%		axis_opt = 1		x-axis label only
%					  2		y-axis label only
%					  3		both axes labels

if nargin < 5
  hdl=gca;
end
EPS=1.e-4;
nfaclong=60/nminlong;
nfaclat=60/nminlat;

if axis_opt ~= 2
    
    if nminlong>0;
        xlim=get(hdl,'xlim');
        xlim(1)=floor(xlim(1)*nfaclong)/nfaclong;
        xtick=xlim(1):1/nfaclong:xlim(2);
        set(hdl,'xtick',xtick);
        xtick0=str2num(char(get(hdl,'xticklabel')));	% added to avoid inconsistency between
        % xticks and xticklabels, D. Chu, 9-3-99
        %% added by D. Chu on 5-11-2001 to correct the inconsistency between
        %% the number of xtick and the number of xticklab
        for i=1:length(xtick0)
            indx=find(abs(xtick-xtick0(i)) < EPS);
            xtick_indx(i)=indx;
        end
        set(hdl,'xtick',xtick(xtick_indx))
        
        % modified 6/5/95 J.List:
        xticklab=degmins(-xtick0,ndiglong);
        set(hdl,'xticklabel',xticklab);
    end;
    
end

if axis_opt ~= 1
    if nminlat>0;
        ylim=get(hdl,'ylim');
        ylim(1)=floor(ylim(1)*nfaclat)/nfaclat;
        ytick=ylim(1):1/nfaclat:ylim(2);
        set(hdl,'ytick',ytick);
        ytick0=str2num(char(get(hdl,'yticklabel')));			% added to avoid inconsistency between
        % yticks and yticklabels, D. Chu, 9-3-99
        %% added by D. Chu on 5-11-2001 to correct the inconsistency between
        %% the number of xtick and the number of xticklab
        for i=1:length(ytick0)
            indx=find( abs(ytick - ytick0(i)) < EPS);
            ytick_indx(i)=indx;
        end
        set(hdl,'ytick',ytick(ytick_indx))
        
        % modified 6/5/95 J.List:
        yticklab=degmins(-ytick0,ndiglat);
        set(hdl,'yticklabel',yticklab,'fontsize',10);
    end
end