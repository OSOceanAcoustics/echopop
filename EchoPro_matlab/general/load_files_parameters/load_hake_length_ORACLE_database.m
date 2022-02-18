function     out=read_hake_length_ORACLE_database(filename,species_code_id,haul_num_offset)
% read hake length data in ORECLE data format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013
global para

if nargin < 3
    haul_num_offset=0;
end

if isempty(filename)
    out=[];
    return
end
d=xlsread(filename);

ind=find(d(:,4) == species_code_id);
%% extract target species
d=d(ind,:);
d(:,3)=d(:,3)+haul_num_offset;

if para.proc. exclude_age1 == 1
    %% exclude age-1 hauls
    [intersect_hauls,IA,IB]= intersect(d(:,3),para.proc.age1_haul);
    if ~isempty(intersect_hauls)
        [selected_hauls,IA,IB]= setxor(d(:,3),intersect_hauls);
        ind=[];
        for i=1:length(IA)
            ind0=find(d(:,3) == selected_hauls(i));
            ind=[ind; ind0];
        end
        d=d(ind,:);
    end
end


[haul, sort_haul_indx]=sort(d(:,3));

d=d(sort_haul_indx,:);
ind=find(diff(haul) > 0);

for i=1:length(ind)+1
    if  i== 1
        if isempty(ind)
            indx=1:length(haul);
        else
            indx=1:ind(i);
        end
    elseif i == length(ind)+1
        indx=ind(i-1)+1:length(haul);
    else
        indx=ind(i-1)+1:ind(i);
    end
    out(i).trawl_no=d(indx(1),3);
    tmp(i).Gender=d(indx,5);
    tmp(i).length=d(indx,6);
    tmp(i).freq=d(indx,7);
end
for i=1:length(out)
    out(i).length=[];
    out(i).Gender=[];
    for j=1:length(tmp(i).length)
        out(i).length=[out(i).length tmp(i).length(j)*ones(1,tmp(i).freq(j))];
        out(i).Gender=[out(i).Gender tmp(i).Gender(j)*ones(1,tmp(i).freq(j))];
    end
    len=out(i).length;
    sex=out(i).Gender;
    out(i).Male_ind=find(sex== 1);
    out(i).Female_ind=find(sex== 2);
    out(i).nM=length(out(i).Male_ind);
    out(i).nF=length(out(i).Female_ind);
    out(i).n=length(len);
    out(i).meanlen=nanmean(len);
    out(i).stdlen=nanstd(len);
    TS0=20*log10(len)-68;
    TS=10*log10(nanmean(10.^(TS0/10)));
    out(i).TS_lin=TS;
    out(i).TS_sd=std(TS0);
    out(i).TS_log=mean(TS0);
end

return