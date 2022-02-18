function     [out, out_all]=load_observer_hake_length_weight_sex_age_data(filename,species_code_id,database_type,excluding_age1_flag,age_data_flag,haul_num_offset)
% read observer hake length, weight, sex, and age data from CSV file recorded with
%% A-Shop data from Vanessa Tutter on June 1, 2018 (N:\Survey.Acoustics\Other Data\A-SHOP data for Chu_060718
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para

if nargin < 6
    haul_num_offset=0;
    age_data_flag=1;
end

[d0, d_str]=xlsread(filename);

if age_data_flag == 1 & size(d0, 2) < 6
   age_data_flag = 2;
   msgbox('Missing age data! Calculate un-aged biomass!', 'NO AGE DATA', 'warn')
end


%% convert gender letters to numbers
if str2num(para.survey_year) > 2007
    gender_col = 5;
else
    gender_col = 6;
end    
for i = 1:length(d0)
    switch char(d_str(i+1,6))
        case 'M'
            d0(i,gender_col)=1;
        case 'F'
            d0(i,gender_col)=2;
        otherwise
            d0(i,gender_col)=3;
    end
end
% for i = 1:length(d0)
%     switch char(d_str(i+1,6))
%         case 'M'
%             d0(i,5)=1;
%         case 'F'
%             d0(i,5)=2;
%         otherwise
%             d0(i,5)=3;
%     end
% end
d = d0;
if str2num(para.survey_year) > 2007
    %% generate ship_ID according to the ascending numbers
    ship_id = cell2mat(d_str(2:end,1));
    ship_haul_id0 = ship_id(:, [1:11 17:20]);
    ship_haul_id_num0 = int64(str2num(ship_haul_id0));
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(ship_haul_id_num0);
    ship_haul_id = ship_haul_id0(sort_ind,:);
    d(:, 1) = ship_haul_id_num;
    d(:, 2:end) =  d0(sort_ind, 2:end);
else
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(d0(:,1));
    d = d0(sort_ind,[1 3:end]);
    ship_haul_id = num2str(ship_haul_id_num);
end

haul = ship_haul_id_num;

ind=find(diff(haul) > 0);
out_all.len=[];
out_all.wgt=[];
out_all.sex=[];
out_all.age=[];
out_all.trawl_no=[];
for i=1:length(ind)+1
    if  i== 1
        indx=1:ind(i);
    elseif i == length(ind)+1
        indx=ind(i-1)+1:length(haul);
    else
        indx=ind(i-1)+1:ind(i);
    end
%     out(i).trawl_no=d(indx(1),1);
    out(i).trawl_no=int64(ship_haul_id_num(indx(1)));  % modified on 11-10-2018 to handle the larger integers
    out(i).sex=d(indx,5)';
    out(i).length=d(indx,3)';
    out(i).wgt=d(indx,4)';
    out(i).age_barcode=d(indx,2)';    
    switch age_data_flag
        case 0
        % faked age data
            out(i).age=max(1,round(3+randn(1,length(indx))));   % mean age = 3 and std = 2
        case 1
        % actual aged data
            out(i).age=d(indx,6)';
        case 2
        % from historical age vs length
        %% 2003-2013 trawl data: C:\Projects\EchoPro\EchoProGUI_Current\other_programs\age_vs_length\age_vs_len_raw.m
        %%   male: age = 4.9553 - 0.3775*len + 0.0089 * len^2,  where len in cm
        %% female: age = 0.0430 - 0.0749*len + 0.0044 * len^2,  where len in cm
        %%  unsex: age = 1.1475 - 0.1427*len + 0.0054 * len^2,  where len in cm (all)
        male_ind=find(out(i).sex == 1);
        female_ind=find(out(i).sex == 2);
        unsex_ind=find(out(i).sex == 3);
        out(i).age(male_ind)  =round(4.9553 - 0.3775*out(i).length(male_ind)   + 0.0089 * out(i).length(male_ind).^2);
        out(i).age(female_ind)=round(0.0431 - 0.0749*out(i).length(female_ind) + 0.0044 * out(i).length(female_ind).^2);
        out(i).age(unsex_ind) =round(1.1475 - 0.1427*out(i).length(unsex_ind) + 0.0054 * out(i).length(unsex_ind).^2);
    end
    %% overall arrays regardless of trawls
    out_all.len=[out_all.len out(i).length];
    out_all.wgt=[out_all.wgt out(i).wgt];
    out_all.sex=[out_all.sex out(i).sex];
    out_all.age=[out_all.age out(i).age];
%     out_all.trawl_no=[out_all.trawl_no out(i).trawl_no*ones(1,length(out(i).length))];
    % modified on 11-10-2018 to handle the larger integers
    out_all.trawl_no=[out_all.trawl_no repmat(out(i).trawl_no, length(out(i).trawl_no), length(out(i).length))];
end

if length(find(~isnan(out_all.age) ==1)) < 0.1*length(out_all.age) 
    fprintf('Warning ------ aged data are less than 10%%!! ------------------\n')
%     out=[];
end  

return