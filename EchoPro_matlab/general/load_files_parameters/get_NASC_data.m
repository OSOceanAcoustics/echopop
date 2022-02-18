function out=get_NASC_data(filename,survey_year,transect_offset)
% load VL interval-based NASC table
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Modification:       4/4/2013
%% Modification:       11/17/2016 new 1995 Export data that only have 8 columns


fprintf('Load Acoustic NASC Data ...\n')

fid=fopen(filename(1,:));
if fid < 0
    out=-1;
    fprintf('file %s not found !!',filename);
    return
end
fclose(fid);

%header=textscan(fid,'%s',1);    % header line
[dat, str_dat,raw]=xlsread(filename);

[n,m]=size(dat);
% if strcmp(survey_year,'1995')  
%  %%%% merge juvenile hake
%  %% merge juvenile only to combined stratum
%    indx=find(dat(:,5) == -1 & dat(:,9) ~= -1);
%    dat(indx,5)=dat(indx,9);
%    dat(indx,9)=-1;
%  %% add juvenile to the bottom of the table
%    indx1=find(dat(:,9) ~= -1);
%    dat1=dat(indx1,:);
%    dat=[dat;dat1];
%    dat(:,9:m)=[];
% end
if str2num(survey_year) < 2003 
    [n,m]=size(dat);
    out(1:n,1)=dat(:,1);                    % transect number   
    ind=find(out(1:n,1) > 1000);            % 1000 is a fixed number added to the first transect of the CAN survey transect
    if ~isempty(ind)
        out(ind,1)= out(ind,1)-1000+transect_offset;             % modify the first transect line number to Tnum+offset
        out(ind+1:end,1)=out(ind+1:end,1)+transect_offset;       % modify the rest CAN transect line numbers to Tnum+offset
    end
    out(1:n,2)=999*ones(n,1);               % region number - not useful
    out(1:n,3)=dat(:,2);                    % VL start 
    out(1:n,4)=dat(:,2)+0.5;                % VL stop
    out(1:n,5:m+2)=dat(:,3:m);              % 
    ind = find(out(1:n,6) > 0  );           % convert longitude to negative
    out(ind,6) = -out(ind,6);
    if str2num(survey_year) == 1995
       ind_nan = find(isnan(out(:,7)) == 1);
       ind_good = find(isnan(out(:,7)) == 0);
       out(ind_nan, 7) = floor(interp1(ind_good, out(ind_good,7), ind_nan));
       ind7 = find(out(:,7) == 7);
       ind7_lt5000 = ind7(ind7 < 5000);
       ind7_gt5000 = ind7(ind7 >= 5000);
       out(ind7_lt5000,7) = 6;
       out(ind7_gt5000,7) = 8;
    end
else
    out=dat; 
end
if size(filename,1)> 1
end
    
return