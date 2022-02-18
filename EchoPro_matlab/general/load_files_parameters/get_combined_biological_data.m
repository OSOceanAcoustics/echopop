function   get_combined_biological_data
% get US & CAN combined biological trawl data
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global  para data

%%% if US & CAN, load US biological data
if  para.proc.source == 3 %  both US&CAN data but read US portion data first
    hake_length_sex_filename=para.bio.filename.length_US;
    hake_length_weight_sex_age_filename=para.bio.filename.specimen_US;
    trawl_filename=para.bio.filename.trawl_US;
    gear_filename=para.bio.filename.gear_US;
    catch_filename=para.bio.filename.catch_US;
else % US or CAN data only
    hake_length_sex_filename=para.bio.filename.length;
    hake_length_weight_sex_age_filename=para.bio.filename.specimen;
    trawl_filename=para.bio.filename.trawl;
    gear_filename=para.bio.filename.gear;
    catch_filename=para.bio.filename.catch;
end
%% Oracle or FSCS data format 
if strcmp(para.bio.database_type,'ORACLE') | strcmp(para.bio.database_type,'Oracle')
    data.bio.hake_length_sex=load_hake_length_ORACLE_database(hake_length_sex_filename,para.bio.species_code_ID);
    data.bio.catch=load_biocatch_ORACLE_database(catch_filename);
else % FSCS format
    data.bio.hake_length_sex=load_hake_length_FSCS_database(hake_length_sex_filename,para.bio.species_code_ID);
    data.bio.catch=load_biocatch_FSCS_database(catch_filename);
end
if ~isempty(trawl_filename)
    data.bio.trawl=load_ORACLE_database_haul_data(trawl_filename);
else
    data.bio.trawl=[];
end
if ~isempty(gear_filename)
    data.bio.gear=load_ORACLE_database_gear_data(gear_filename);
else
    data.bio.gear=[];
end

[out1,out2]=load_hake_length_weight_sex_age_data(hake_length_weight_sex_age_filename,para.bio.species_code_ID,para.bio.database_type,para.proc.exclude_age1,para.proc.age_data_status,0);
data.bio.hake_length_weight_sex_age=out1;

%% combine trawl data collected at both stations 1 and 2
data.bio.len_wgt_all=out2;

%%% if US & CAN, load Canadian biological data
if para.proc.source == 3 % both US&CAN data but read CAN portion data now since US portion data have been loaded already
    hake_length_sex_filename=para.bio.filename.length_CAN;
    hake_length_weight_sex_age_filename=para.bio.filename.specimen_CAN;
    catch_filename=para.bio.filename.catch_CAN;
    %% Oracle or FSCS data format 
    if strcmp(para.bio.database_type,'ORACLE') | strcmp(para.bio.database_type,'Oracle')
        out=load_hake_length_ORACLE_database(hake_length_sex_filename,para.bio.species_code_ID,para.bio.haul_no_offset);
        out_catch=load_biocatch_ORACLE_database(catch_filename);
    else   % FSCS format
        out=load_hake_length_FSCS_database(hake_length_sex_filename,para.bio.species_code_ID,para.bio.haul_no_offset);
        out_catch=load_biocatch_FSCS_database(catch_filename);
    end
    [out1,out2]=load_hake_length_weight_sex_age_data(hake_length_weight_sex_age_filename,para.bio.species_code_ID,para.bio.database_type,para.proc.exclude_age1,para.proc.age_data_status,para.bio.haul_no_offset);
%% all length data from station #1 and 2 of US and CAN trawl data
    combine_two_biological_data_output(out,out1,out2);
    trawl_filename=para.bio.filename.trawl_CAN;
    gear_filename=para.bio.filename.gear_CAN;
    out_gear=load_ORACLE_database_gear_data(gear_filename);
    %% deal with trawl number offset
    max_US_haul_no=max(data.bio.trawl.trawl_no);
    if max_US_haul_no > para.bio.haul_no_offset
        CAN_haul_no_offset = 100*ceil(max_US_haul_no/100);
    else
        CAN_haul_no_offset = para.bio.haul_no_offset;
    end
    max_US_Transect=max(data.bio.gear.transect);
%     if max_US_Transect > para.proc.transect_offset & min(out_gear.transect) < max_US_Transect  % transect overlap
%         CAN_Transect_offset = 100*ceil(max_US_Transect/100);
%     else
%         CAN_Transect_offset = para.proc.transect_offset;
%     end
%%  transect offset is set to zero since we will not have overlap transects
    CAN_Transect_offset=0;
% combine US & CAN trawl files
    out_trawl=load_ORACLE_database_haul_data(trawl_filename);
    if isstruct(out_trawl)
        out_trawl.trawl_no=out_trawl.trawl_no+CAN_haul_no_offset;
        ind=find(out_trawl.EQlon> 0);
        out_trawl.EQlon(ind)=-out_trawl.EQlon(ind);   % change lonW to -lon
        data.bio.trawl=merge_two_struct(data.bio.trawl,out_trawl);
    end
% combine US & CAN gear files    
    if isstruct(out_gear)
        out_gear.trawl_no=out_gear.trawl_no+CAN_haul_no_offset;
        out_gear.transect=out_gear.transect+CAN_Transect_offset;
        data.bio.gear=merge_two_struct(data.bio.gear,out_gear);
    end
% combine US & CAN catch files
    if isstruct(out_catch)
        n=length(data.bio.catch);
        n1=length(out_catch);
        for i=n+1:n+n1
            data.bio.catch(i).trawl_no=out_catch(i-n).trawl_no+CAN_haul_no_offset;
            data.bio.catch(i).species=out_catch(i-n).species;
        end
    end
end



