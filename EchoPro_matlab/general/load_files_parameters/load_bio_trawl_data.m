function load_bio_trawl_data
% load trawl biological data in ORECLE dataformat
%
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013


global para data 
% read catch data
if strcmp(para.bio.database_type,'ORACLE') | strcmp(para.bio.database_type,'Oracle')
  data.bio.catch=load_biocatch_ORACLE_database(para.bio.filename.catch);
else
  data.bio.catch=load_biocatch_FSCS_database(para.bio.filename.catch);
end

%% real length-sex
if strcmp(para.bio.database_type,'ORACLE') | strcmp(para.bio.database_type,'Oracle')
    data.bio.hake_length_sex=load_hake_length_ORACLE_database(para.bio.filename.length,para.bio.species_code_ID);
else
    data.bio.hake_length_sex=load_hake_length_FSCS_database(para.bio.filename.length,para.bio.species_code_ID);
end
% read length-weight-sex-age
[out1,out2]=load_hake_length_weight_sex_age_data(para.bio.filename.specimen,para.bio.species_code_ID,para.bio.database_type,para.proc.exclude_age1,para.proc.age_data_status);
data.bio.hake_length_weight_sex_age=out1;
data.bio.len_wgt_all=out2;

%% load transect-haul file
out=load_transect_trawl_info(para.bio.filename.transect_haul);
date.bio.TX=out.TX;
proc_len_wgt_data;
get_strata_data;
