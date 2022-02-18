function update_filename_window_settings
%% update filenames based on the settings in 'Reload Parameter Files' windows
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013
global hdl para data

h=hdl. redefine_filenames;

survey_yr_cell=get(h.popup_survey_year,'string');
str_sy=char(survey_yr_cell{get(h.popup_survey_year,'value')});
para.survey_year=str_sy(1:4);

if get(h.radio_US,'value') == 1
    para.proc.source=1;        % US data
elseif get(h.radio_CAN,'value') == 1
    para.proc.source=2;        % Canadian data
else
    para.proc.source=3;        % US & Canadian data    
end
para_setting_filename=['proc_parameters_' para.survey_year];
%disp(para_setting_filename)
eval(para_setting_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Biological data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(get(h.pb_bio_gear,'UserData'))
    para.bio.filename.gear=get(h.pb_bio_gear,'UserData');
end
if ~isempty(get(h.pb_bio_trawl,'UserData'))
    para.bio.filename.trawl=get(h.pb_bio_trawl,'UserData');
end
if ~isempty(get(h.pb_bio_catch,'UserData'))
    para.bio.filename.catch=get(h.pb_bio_catch,'UserData');
end
if ~isempty(get(h.pb_bio_length,'UserData'))
    para.bio.filename.length=get(h.pb_bio_length,'UserData');
end
if ~isempty(get(h.pb_bio_specimen,'UserData'))
    para.bio.filename.specimen=get(h.pb_bio_specimen,'UserData');
end
if ~isempty(get(h.pb_bio_transect_haul,'UserData'))
    para.bio.filename.transect_haul=get(h.pb_bio_transect_haul,'UserData');
end
if ~isempty(get(h.pb_bio_strata,'UserData'))
    para.acoust.filename.strata=get(h.pb_bio_strata,'UserData');
end
% data.bio.trawl=load_ORACLE_database_haul_data(para.bio.filename.trawl);
% data.bio.gear=load_ORACLE_database_gear_data(para.bio.filename.gear);
% get_combined_biological_data  %replacing load_bio_trawl_data
% construct_catch_trawl_output_matrices



if get(h.radio_Oracle,'value') == 1
    para.bio.database_type='Oracle';
else 
    para.bio.database_type='FSCS';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Biological data parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(h.radio_acoust_EK60_dir,'value') == 1       % EK60
    para.acoust.file_sys=1;
    if ~isempty(get(h.pb_acoust_datatype,'Userdata'))
        para.acoust.dir.EK60_root=get(h.pb_acoust_datatype,'UserData');
    end
elseif get(h.radio_acoust_EK500_dir,'value') == 1       % EK500
    para.acoust.file_sys=2;
    if ~isempty(get(h.pb_acoust_datatype,'Userdata'))
        para.acoust.dir.EK500_root=get(h.pb_acoust_datatype,'UserData');
    end
else
   para.acoust.file_sys=3;
   if ~isempty(get(h.pb_acoust_datatype,'Userdata'))
       para.acoust.filename.processed_data=get(h.pb_acoust_datatype,'UserData');
   end
end

if ~isempty(get(h.pb_acoust_cal,'Userdata'))
   para.acoust.filename.cal=get(h.pb_acoust_cal,'Userdata');
end

if ~isempty(get(h.pb_acoust_evr,'Userdata'))
   para.acoust.dir.evr_root=get(h.pb_acoust_evr,'Userdata');
end

if get(h.radio_acoust_species_mix,'value') ==  0
    para.acoust.mix_region = 0;
    if ~isempty(get(h.pb_acoust_transect_region,'Userdata'))
       para.acoust.filename.transect_region=get(h.pb_acoust_transect_region,'Userdata');
    end
else
    para.acoust.mix_region = 1;
    if ~isempty(get(h.pb_acoust_transect_mix_region,'Userdata'))
       para.acoust.filename.transect_mix_region=get(h.pb_acoust_transect_mix_region,'Userdata');
    end
end

if ~isempty(get(h.pb_acoust_VL_spacing,'Userdata'))
   para.acoust.filename.VL_spacing=get(h.pb_acoust_VL_spacing,'Userdata');
end

if ~isempty(get(h.pb_acoust_VL_lat_lon,'Userdata'))
   para.acoust.dir.VL_LatLon=get(h.pb_acoust_pb_VL_lat_lon,'Userdata');
end

if ~isempty(get(h.pb_acoust_Ev_bot,'Userdata'))
   para.acoust.dir.transect_bot=get(h.pb_acoust_Ev_bot,'Userdata');
end

if str2num(para.survey_year) > 2001 & get(hdl.redefine_filenames.radio_acoust_processed,'value') == 0
% read transect-trawl-region data
    if para.acoust.mix_region == 0
       xtr=load_xtr_info(para.acoust.filename.transect_region);
    else
       xtr=load_xtr_info_rev(para.acoust.filename.transect_region);
       data.acoust.mix_trawls=load_mix_trawls_info(para.acoust.filename.transect_mix_region);
    end
    data.acoust.transect_trawl_region=xtr;
    out=load_cal_table(para.acoust.filename.cal,para.acoust.freq_ind);
    data.acoust.cal=out;
end


return