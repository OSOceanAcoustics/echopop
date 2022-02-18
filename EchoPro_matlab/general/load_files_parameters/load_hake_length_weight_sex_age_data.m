function     [out, out_all]=load_hake_length_weight_sex_age_data(filename,species_code_id,database_type,excluding_age1_flag,age_data_flag,haul_num_offset)
% read hake length, weight, sex, and age data from CSV file recorded with
% ORACLE datbase format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para

if nargin < 6
    haul_num_offset=0;
    age_data_flag=1;
end
d=xlsread(filename);

ind=find(d(:,4) == species_code_id);
%% extract target species
d=d(ind,:);
d(:,3)=d(:,3)+haul_num_offset;

if para.proc.exclude_age1 == 1
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
out_all.len=[];
out_all.wgt=[];
out_all.sex=[];
out_all.age=[];
out_all.trawl_no=[];
for i=1:length(ind)+1
%     if i == 30
%         disp(i)
%     end
    if  i== 1
        indx=1:ind(i);
    elseif i == length(ind)+1
        indx=ind(i-1)+1:length(haul);
    else
        indx=ind(i-1)+1:ind(i);
    end
    out(i).trawl_no=d(indx(1),3);
    out(i).sex=d(indx,5)';
    out(i).length=d(indx,6)';
    out(i).wgt=d(indx,7)';
    out(i).age_barcode=d(indx,8)';    
    switch age_data_flag
        case 0  % not used anymore
        % faked age data
            out(i).age=max(1,round(3+randn(1,length(indx))));   % mean age = 3 and std = 2
        case 1
        % actual aged data
            out(i).age=d(indx,9)';
        case 2   % selected method
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
    out_all.trawl_no=[out_all.trawl_no out(i).trawl_no*ones(1,length(out(i).length))];
end

if length(find(~isnan(out_all.age) ==1)) < 0.1*length(out_all.age) 
    fprintf('Warning ------ aged data are less than 10%%!! ------------------\n')
%     out=[];
end  

return