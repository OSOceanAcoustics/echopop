function   out=load_biocatch_ORACLE_database(filename)
% load biocatch data in ORECLE data format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013


fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Catch file: %s !!',filename));
    return
end

[d,s]=xlsread(filename);

haul=d(:,3);
ind=find(diff(haul) > 0);

for i=1:length(ind)+1
    if  i== 1
        indx=1:ind(i);
    elseif i == length(ind)+1
        indx=ind(i-1)+1:length(haul);
    else
        indx=ind(i-1)+1:ind(i);
    end
    out(i).trawl_no=d(indx(1),3);
    for j=1:length(indx)
        out(i).species(j).ID=d(indx(j),4);
        out(i).species(j).name=d(indx(j),5);
        if isnan(out(i).species(j).name)
            out(i).species(j).name=s(indx(j)+1,5);      % consider the header line
        end
        out(i).species(j).exp_cnt=d(indx(j),6);        
        out(i).species(j).exp_wgt=d(indx(j),7);        
    end
end

return
%out.filename=filename;
% line=fgetl(fid);            % information line
% 
% i=0;
% j=0;
% line=fgetl(fid);
% while line > 0
% % process the first line of any trawl
%    ind=find(line == ',');
%    trawl_num=str2num(line(ind(2)+1:ind(3)-1));
%    if i ~= 0 & out(i).trawl_no == trawl_num
%            j=j+1;
%    else
%       j=1;
%       i=i+1;
%       out(i).trawl_no=trawl_num;
%    end
%    out(i).species(j).ID=line(ind(3)+1:ind(4)-1);
%    out(i).species(j).name=line(ind(4)+1:ind(5)-1);
%    cnt=line(ind(5)+1:ind(6)-1);
%    if ~isempty(cnt)
%      out(i).species(j).exp_cnt=cnt;
%    else
%      out(i).species(j).exp_cnt=nan;
%    end
%    wgt=str2num(line(ind(6)+1:ind(7)-1));
%    if ~isempty(wgt)
%       out(i).species(j).exp_wgt=wgt;   
%    else
%       out(i).species(j).exp_wgt=nan;   
%    end
%    line=fgetl(fid);
% end
% 
% fclose(fid);