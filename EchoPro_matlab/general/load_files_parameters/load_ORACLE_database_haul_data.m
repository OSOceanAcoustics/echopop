function out=load_ORACLE_database_haul_data(filename)
% load haul data saved in Oracle database format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para data

fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Haul file: %s !!',filename));
    return
end
fclose(fid);

[dat,dat_str]=xlsread(filename);
ind=find(isnan(dat(:,18)) == 1);
if  length(ind) == size(dat,1) 
    dat(:,18)=str2num(char(dat_str(2:end,18)));  % for some reasons, the column 'AVERAGE_BOTTOM_DEPTH' is read as 'string' instead of number
end

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

out.trawl_no=dat(:,3);
ind=find(diff(out.trawl_no)< -10);
out.trawl_no(ind+1:end)=out.trawl_no(ind+1:end)+para.bio.haul_no_offset;
out.haultype=dat(:,7);
out.performance=dat(:,8);
out.duration=dat(:,11);
out.distance=dat(:,12);
out.stratum=dat(:,13);
if para.proc.hemisphere(1) == 'N'
   out.EQlat=abs(dat(:,14));
elseif para.proc.hemisphere(1) == 'S'
   out.EQlat=-abs(dat(:,14));
else
   fprintf('Wrong N/S Hemisphere definition!!\n');
end
if para.proc.hemisphere(2) == 'W'
   out.EQlon=-abs(dat(:,15));
elseif para.proc.hemisphere(2) == 'E'
   out.EQlon=abs(dat(:,15));
else
   fprintf('Wrong E/W Hemisphere definition!!\n');
end
out.ave_dep=dat(:,18);
if size(dat,2) >= 24
    out.VLstart=dat(:,23);
    out.VLstop=dat(:,24);
else
    out.VLstart=nan*ones(size(dat,1),1);
    out.VLstop=out.VLstart;
end
return