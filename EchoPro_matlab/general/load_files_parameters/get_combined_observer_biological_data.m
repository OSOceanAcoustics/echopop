function   get_combined_observer_biological_data
% get US observer combined biological trawl data
%% A-Shop data from Vanessa Tutter on June 1, 2018 (N:\Survey.Acoustics\Other Data\A-SHOP data for Chu_060718)
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/23/2018

global  para data

%%% if US & CAN, load US biological data
if  para.proc.source <= 3 %  both US&CAN data but read US portion data first
    hake_length_sex_filename=para.bio.filename.length;
    hake_length_weight_sex_age_filename=para.bio.filename.specimen;
    trawl_filename=para.bio.filename.trawl;
    gear_filename=para.bio.filename.gear;
    catch_filename=para.bio.filename.catch;
end

if ~isempty(trawl_filename)
    data.bio.trawl=load_observer_haul_database(trawl_filename);
else
    data.bio.trawl=[];
end

data.bio.hake_length_sex=load_observer_hake_length_database(hake_length_sex_filename,para.bio.species_code_ID);
data.bio.catch=load_observer_biocatch_database(catch_filename);

if ~isempty(gear_filename)
    data.bio.gear=load_observer_database_gear_data(gear_filename);
else
    data.bio.gear=[];
end

[out1,out2]=load_observer_hake_length_weight_sex_age_data(hake_length_weight_sex_age_filename,para.bio.species_code_ID,para.bio.database_type,para.proc.exclude_age1,para.proc.age_data_status,0);
data.bio.hake_length_weight_sex_age=out1;

%% combine trawl data collected at both stations 1 and 2
data.bio.len_wgt_all=out2;

