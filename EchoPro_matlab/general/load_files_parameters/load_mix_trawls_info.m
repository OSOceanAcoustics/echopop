function   out=load_mix_trawls_info(filename)
% load mix trawls information data
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Mix Tralws file: %s !!',filename));
    return
end

out.filename=filename;
line=fgetl(fid);            % header - first line
i=1;
while ~feof(fid)
    line=fgetl(fid);
    ind=find(line == ',');
    out.stratum(i)=str2num(line(1:ind(1)-1));
    out.trawl_no(i)=str2num(line(ind(1)+1:ind(2)-1));
    out.VL_start(i)=str2num(line(ind(2)+1:ind(3)-1));
    out.VL_stop(i)=str2num(line(ind(3)+1:ind(4)-1));
    out.TX_no(i)=str2num(line(ind(4)+1:ind(5)-1));
    out.region_ID(i)=str2num(line(ind(5)+1:ind(6)-1));
    out.species1_sa_proportion(i)=str2num(line(ind(6)+1:ind(7)-1))/100;
    if length(ind) == 7
       out.species2_sa_proportion(i)=str2num(line(ind(7)+1:end))/100;
    else
       out.species2_sa_proportion(i)=str2num(line(ind(7)+1:ind(8)-1))/100;
    end
    i=i+1;
end
fclose(fid);
return