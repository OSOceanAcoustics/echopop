function   out=load_xtr_info_rev(filename)
% load transect-trawl-region file that contain hake region information data
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013


fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Transect-Trawl-Region file: %s !!',filename));
    return
end

out.filename=filename;
i=1;
while ~feof(fid)
    out(i).Xnum=fgetl(fid);
    line=fgetl(fid);
    ind=find(line == ',');
    j=1;
    while line > 0
        out(i).region(j).num=str2num(line(1:ind(1)-1));
        out(i).region(j).n_trawls=str2num(line(ind(1)+1:ind(2)-1));
        for k=1:out(i).region(j).n_trawls
            out(i).region(j).trawl_no(k)=str2num(line(ind(k+1)+1:ind(k+2)-1));
        end
        out(i).region(j).class=line(ind(k+1)+1:end);
        j=j+1;
        line=fgetl(fid);
        ind=find(line == ',');
    end
    i=i+1;
end
fclose(fid);