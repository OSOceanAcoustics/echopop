function load_bio_catch_data
% load biological catch data in ORECLE/FSCS data format
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

