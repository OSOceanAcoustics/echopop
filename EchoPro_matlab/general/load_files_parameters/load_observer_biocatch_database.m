function   out=load_observer_biocatch_database(filename)
% load observer biocatch data 
%% A-Shop data from Vanessa Tutter on June 1, 2018 (N:\Survey.Acoustics\Other Data\A-SHOP data for Chu_060718
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/23/2013

global para data
fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Catch file: %s !!',filename));
    return
end

[d0,d_str0]=xlsread(filename);
if str2num(para.survey_year) > 2007
    d = d0;
    %% generate ship_ID according to the ascending numbers
    ship_id = cell2mat(d_str0(2:end,1));
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
    d_str = d_str0(sort_ind+1, :);              % header line has been removed
else
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(d0(:,1));
    d = d0(sort_ind,[1 3:end]);
    ship_haul_id = num2str(ship_haul_id_num);
    d_str = d_str0(sort_ind+1, :);              % header line has been removed
end
hauls = ship_haul_id_num;

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

ind=find(diff(hauls) > 0);
for i=1:length(ind)+1
    if  i== 1
        indx=1:ind(i);
    elseif i == length(ind)+1
        indx=ind(i-1)+1:length(hauls);
    else
        indx=ind(i-1)+1:ind(i);
    end
    out(i).trawl_no=d(indx(1),1);
    for j=1:length(indx)
        out(i).species(j).ID=d(indx(j),2);
        out(i).species(j).name=d(indx(j),3);
        if isnan(out(i).species(j).name)
            out(i).species(j).name=d_str(indx(j),4);      % consider the header line
        end
        out(i).species(j).exp_cnt=d(indx(j),4);        
        out(i).species(j).exp_wgt=d(indx(j),5);        
    end
end

return
