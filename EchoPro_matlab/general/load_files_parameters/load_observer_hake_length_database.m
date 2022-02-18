function     out=read_observer_hake_length_database(filename,species_code_id,haul_num_offset)
% read observer hake length data 
%% A-Shop data from Vanessa Tutter on June 1, 2018 (N:\Survey.Acoustics\Other Data\A-SHOP data for Chu_060718
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/23/2018

global para

if nargin < 3
    haul_num_offset=0;
end

if isempty(filename)
    out=[];
    return
end
[d0, d_str]=xlsread(filename);

%% convert gender letters to numbers
if str2num(para.survey_year) > 2007
    gender_col = 4;
else
    gender_col = 5;
end    
for i = 1:length(d0)
    switch char(d_str(i+1,5))
        case 'M'
            d0(i,gender_col)=1;
        case 'F'
            d0(i,gender_col)=2;
        otherwise
            d0(i,gender_col)=3;
    end
end

if str2num(para.survey_year) > 2007
    d = d0;
    %% generate ship_ID according to the ascending numbers
    ship_id = cell2mat(d_str(2:end,1));
    haul_id = num2str2(d0(:,1), 4);
    % ship_haul_id0 = [ship_id haul_id];
    % ship_haul_id_num0 = str2num(ship_haul_id0);
    ship_haul_id0 = ship_id(:, [1:11 17:20]);
    ship_haul_id_num0 = int64(str2num(ship_haul_id0));    
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(ship_haul_id_num0);
    ship_haul_id = ship_haul_id0(sort_ind,:);
    d(:, 1) = ship_haul_id_num;
    d(:, 2:end) =  d0(sort_ind, 2:end);
    ship_haul_id_num = int64(str2num(ship_haul_id));
else
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(d0(:,1));
    d = d0(sort_ind,[1 3:end]);
    ship_haul_id = num2str(ship_haul_id_num);   
end

if para.proc.exclude_age1 == 1
    %% exclude age-1 hauls
    [intersect_hauls,IA,IB]= intersect(d(:,1),para.proc.age1_haul);
    if ~isempty(intersect_hauls)
        [selected_hauls,IA,IB]= setxor(d(:,1),intersect_hauls);
        ind=[];
        for i=1:length(IA)
            ind0=find(d(:,1) == selected_hauls(i));
            ind=[ind; ind0];
        end
        d=d(ind,:);
    end
end


[haul, sort_haul_indx]=sort(d(:,1));

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
%     out(i).trawl_no=d(indx(1),1);
    out(i).trawl_no=ship_haul_id_num(indx(1));   % modified on 11-10-2018 to handle the larger integers
    tmp(i).Gender=d(indx,4);
    tmp(i).length=d(indx,2);
    tmp(i).freq=d(indx,3);
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