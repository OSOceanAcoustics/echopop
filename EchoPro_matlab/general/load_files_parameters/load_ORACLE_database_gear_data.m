function out=load_ORACLE_database_gear_data(filename)
% load trawl gear data saved in Oracle database format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para

fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Gear file: %s !!',filename));
    return
end
fclose(fid);

dat=xlsread(filename);

if para.proc. exclude_age1 == 1
    %% exclude age-1 hauls
    [intersect_hauls,IA,IB]= intersect(dat(:,3),para.proc.age1_haul);
    if ~isempty(intersect_hauls)
        [selected_hauls,IA,IB]= setxor(dat(:,3),intersect_hauls);
        ind=[];
        for i=1:length(IA)
            ind0=find(dat(:,3) == selected_hauls(i));
            ind=[ind; ind0];
        end
        dat=dat(ind,:);
    end
end

out.transect=dat(:,8);
%% commented out the next two lines for 2011 processed data
if str2num(char(para.survey_year))  ~= 2011 & str2num(char(para.survey_year))  ~= 2015 & str2num(char(para.survey_year))  ~= 1995
    ind=find(diff(out.transect)< -10);
    out.transect(ind+1:end)=out.transect(ind+1:end)+out.transect(ind);
end
%% for specified transects, all trawls not on the specified transects will still be included, modified on 12-2-2012
% ind=find(out.transect >= para.proc.start_transect &  out.transect <= para.proc.end_transect);

ind=1:length(out.transect);
out.transect=out.transect(ind);
dat=dat(ind,:);
out.trawl_no=dat(:,3);
ind=find(diff(out.trawl_no)< -10);
out.trawl_no(ind+1:end)=out.trawl_no(ind+1:end)+para.bio.haul_no_offset;
out.footropeD=dat(:,4);
out.Tsurf=dat(:,5);
out.Tdep=dat(:,6);
out.ave_wireout=dat(:,7);
out.ave_netopen=dat(:,10);

return