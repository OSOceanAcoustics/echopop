function update_data_bio_struct(TX_selected, data0)
%% update data.bio struct for every bootstrap realization
%% created   3/9/2020

global data para

CAN_trawl_ind = find(TX_selected > 200);
TX_selected(CAN_trawl_ind) = TX_selected(CAN_trawl_ind) - para.bio.haul_no_offset;

j0_ls = 1;   % start index of len_sex
j0_lwsa = 1; % start index of len_wgt_sex_age
j0_c = 1;    % start index of catch
j0_t = 1;    % start index of trawl
j0_g = 1;    % start index of gear
cnt_ls = 1;   % conut index of len_sex
cnt_lwsa = 1; % conut index of len_wgt_sex_age
cnt_c = 1;    % conut index of catch
cnt_t = 1;    % conut index of trawl
cnt_g = 1;    % conut index of gear
%% data.bio struct initialization
data.bio.hake_length_sex = [];
data.bio.hake_length_weight_sex_age = [];
data.bio.catch = [];
data.bio.trawl = [];
data.bio.gear = [];

%% updated gear and possible haul numbers
for i = 1:length(TX_selected)
    for j = j0_g:length(data0.bio.gear.trawl_no)
        if data0.bio.gear.transect(j) == TX_selected(i)
            data_bio.gear.trawl_no(cnt_g) = data0.bio.gear.trawl_no(j);
            data_bio.gear.transect(cnt_g) = data0.bio.gear.transect(j);
            data_bio.gear.footropeD(cnt_g) = data0.bio.gear.footropeD(j);
            data_bio.gear.Tsurf(cnt_g) = data0.bio.gear.Tsurf(j);
            data_bio.gear.Tdep(cnt_g) = data0.bio.gear.Tdep(j);
            data_bio.gear.ave_wireout(cnt_g) = data0.bio.gear.ave_wireout(j);
            data_bio.gear.ave_netopen(cnt_g) = data0.bio.gear.ave_netopen(j);
            %% get the trawl_no
            TX_selected_trawl_no_indx(cnt_g) = data_bio.gear.trawl_no(cnt_g);
            j0_g = j + 1;
            cnt_g = cnt_g + 1;
        end
    end
end

data.tmp.selected_TX = TX_selected;
data.tmp.TX_selected_trawl_no_indx = TX_selected_trawl_no_indx;

%% find hake trawls
for i = 1:length(TX_selected_trawl_no_indx)
%% updated hake_length_sex
    found_trawl = 0;
    for j = j0_ls:length(data0.bio.hake_length_sex)        
        if data0.bio.hake_length_sex(j).trawl_no == TX_selected_trawl_no_indx(i)
            data_bio.hake_length_sex = data0.bio.hake_length_sex(j);
            j0_ls = j + 1;
            found_trawl = 1;
            break
        end
    end
    if j0_ls == j + 1 & found_trawl == 1
        if  cnt_ls == 1
            data.bio.hake_length_sex = data_bio.hake_length_sex;
        else
            data.bio.hake_length_sex(cnt_ls) = data_bio.hake_length_sex;
        end
%         fprintf('cnt = %d\t len_sex -> trawl is found on transect: %d \n', cnt_ls, TX_selected_trawl_no_indx(i));
        cnt_ls = cnt_ls + 1;
    else
%         fprintf('cnt = %d\t len_sex -> no trawl is found on transect: %d ... \n', cnt_ls, TX_selected_trawl_no_indx(i));
    end
    
%% updated hake_length_weight_sex_age
    found_trawl = 0;
    for j = j0_lwsa:length(data0.bio.hake_length_weight_sex_age)
        if data0.bio.hake_length_weight_sex_age(j).trawl_no == TX_selected_trawl_no_indx(i)
            data_bio.hake_length_weight_sex_age = data0.bio.hake_length_weight_sex_age(j);
            j0_lwsa = j + 1;
            found_trawl = 1;
            break
        end
    end
    if j0_lwsa == j + 1 & found_trawl == 1
        if  cnt_lwsa == 1
            data.bio.hake_length_weight_sex_age = data_bio.hake_length_weight_sex_age;
        else
            data.bio.hake_length_weight_sex_age(cnt_lwsa) = data_bio.hake_length_weight_sex_age;
        end
%         fprintf('cnt = %d\t len_wgt_sex_age -> trawl is found on transect: %d \n', cnt_lwsa, TX_selected_trawl_no_indx(i));
        cnt_lwsa = cnt_lwsa + 1;
    else
%         fprintf('cnt = %d\t len_wgt_sex_age -> no trawl is found on transect: %d ... \n', cnt_lwsa, TX_selected_trawl_no_indx(i));
    end
    
%% updated catch
    found_trawl = 0;
    for j = j0_c:length(data0.bio.catch)
        if data0.bio.catch(j).trawl_no == TX_selected_trawl_no_indx(i)
            data_bio.catch = data0.bio.catch(j);
            j0_c = j + 1;
            found_trawl = 1;
            break
        end
    end
    if j0_c == j + 1 & found_trawl == 1
        if  cnt_c == 1 
            data.bio.catch = data_bio.catch;
        else
            data.bio.catch(cnt_c) = data_bio.catch;
        end
%         fprintf('cnt = %d\t catch -> no trawl is found on transect: %d \n', cnt_c, TX_selected_trawl_no_indx(i));
        cnt_c = cnt_c + 1;
    else
%         fprintf('cnt = %d\t catch -> no trawl is found on transect: %d ... \n', cnt_c, TX_selected_trawl_no_indx(i));
    end
    
%% updated trawl
    for j = j0_t:length(data0.bio.trawl.trawl_no)
        if data0.bio.trawl.trawl_no(j) == TX_selected_trawl_no_indx(i)
            data_bio.trawl.trawl_no(cnt_t) = TX_selected_trawl_no_indx(i);
            data_bio.trawl.haultype(cnt_t) = data0.bio.trawl.haultype(j);
            data_bio.trawl.performance(cnt_t) = data0.bio.trawl.performance(j);
            data_bio.trawl.duration(cnt_t) = data0.bio.trawl.duration(j);
            data_bio.trawl.distance(cnt_t) = data0.bio.trawl.distance(j);
            data_bio.trawl.stratum(cnt_t) = data0.bio.trawl.stratum(j);
            data_bio.trawl.EQlat(cnt_t) = data0.bio.trawl.EQlat(j);
            data_bio.trawl.EQlon(cnt_t) = data0.bio.trawl.EQlon(j);
            data_bio.trawl.ave_dep(cnt_t) = data0.bio.trawl.ave_dep(j);
            data_bio.trawl.VLstart(cnt_t) = data0.bio.trawl.VLstart(j);
            data_bio.trawl.VLstop(cnt_t) = data0.bio.trawl.VLstop(j);
            cnt_t = cnt_t + 1;
            j0_t = j + 1;
            break
        end
    end
end
%% update the trawl and gear structs
data.bio.trawl = data_bio.trawl;
data.bio.gear = data_bio.gear;

return