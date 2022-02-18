function  proc_parameters_2021
%% this file is to provide all required input filenames & some process parameter settings
%%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       8/02/2021

global para data

% para.proc.source=3;         % 1 = US; 2 CAN; 3 = all
if para.proc.source == 1     %% US only
    %% Biological trawl data
    para.bio.filename.trawl=['' para.data_root_dir 'Biological\2021\US\2021_biodata_haul.xlsx'];
    para.bio.filename.gear=['' para.data_root_dir 'Biological\2021\US\2021_biodata_gear.xlsx'];
    para.bio.filename.catch=['' para.data_root_dir 'Biological\2021\US\2021_biodata_catch.xlsx'];
    para.bio.filename.length=['' para.data_root_dir 'Biological\2021\US\2021_biodata_length.xlsx'];
    para.bio.filename.specimen=['' para.data_root_dir 'Biological\2021\US\2021_biodata_specimen.xlsx'];
%     para.bio.filename.trawl=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\202106_biodata_haul.xlsx'];
%     para.bio.filename.gear=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\202106_biodata_gear.xlsx'];
%     para.bio.filename.catch=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\202106_biodata_catch.xlsx'];
%     para.bio.filename.length=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\202106_biodata_length.xlsx'];
%     para.bio.filename.specimen=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\202106_biodata_specimen.xlsx'];
    %% Stratification files
    para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2021\US strata 2021_11-Nov-2021.xlsx'];
    para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\2021\US20210927_age2+_exports.xlsx'];
    para.acoust.filename.processed_data_age2=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\US&CAN_detailsa_2021_table2y+_ALL_11-9-2021.xlsx'];
    para.acoust.filename.processed_data_age1=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\US&CAN_detailsa_2021_table1y+_ALL_11-9-2021.xlsx'];
    para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2021\Stratification_geographic_Lat_2021_final.xlsx'];
    %     para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2021\US&CAN_2021_US_transect_region_haul_age2+ auto_20210927.xlsx'];
    %     para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2021\US&CAN_2021_transect_region_haul_age1+ auto_20211029.xlsx'];
    para.bio_acoust.filename.Transect_region_haul_age2=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\US&CAN_2021_transect_region_haul_age1+ auto_final.xlsx'];
    para.bio_acoust.filename.Transect_region_haul_age1=['C:\Projects\EchoPro\EchoProGUI_Currrent\outputs\US&CAN_2021_transect_region_haul_age1+ auto_final.xlsx'];
    if para.proc.exclude_age1 == 1
        para.acoust.filename.processed_data=para.acoust.filename.processed_data_age2;
        para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age2;
    else
        para.acoust.filename.processed_data=para.acoust.filename.processed_data_age1;
        para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age1;
    end
    %% kriging related files
    para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2021\default_vario_krig_settings_2021_US.xlsx'];
    %     data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\US_only_krigedgrid2_5nm_2021_SHtransects_toChu.xlsx'];     % 2021 US cell res = 2.50 nmi with extended area coverage
    data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
elseif para.proc.source == 2   % CAN only
    para.bio.filename.trawl=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_haul_CAN.xlsx'];
    para.bio.filename.gear=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_gear_CAN.xlsx'];
    para.bio.filename.catch=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_catch_CAN.xlsx'];
    para.bio.filename.length=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_length_CAN.xlsx'];
    para.bio.filename.specimen=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_specimen_CAN_AGES.xlsx'];
    fprintf('******* No CAN stratification files *******\n')
    fprintf('******* No CAN Transsect-region-haul files *******\n')
    fprintf('******* No CAN NASC export data files *******\n')
    fprintf('******* No CAN vario-krig parameter file *******\n')
    fprintf('******* No CAN grid cell file *******\n\n')
elseif para.proc.source == 3 % US & CAN acoustic trawl data
    switch para.bio_data_type
        case 1   % Acoustic and Trawl survey data
            %% Biological trawl data
            para.bio.filename.trawl_US=['' para.data_root_dir 'Biological\2021\US\2021_biodata_haul.xlsx'];
            para.bio.filename.gear_US=['' para.data_root_dir 'Biological\2021\US\2021_biodata_gear.xlsx'];
            para.bio.filename.catch_US=['' para.data_root_dir 'Biological\2021\US\2021_biodata_catch.xlsx'];
            para.bio.filename.length_US=['' para.data_root_dir 'Biological\2021\US\2021_biodata_length.xlsx'];
            para.bio.filename.specimen_US=['' para.data_root_dir 'Biological\2021\US\2021_biodata_specimen.xlsx'];
            para.bio.filename.trawl_CAN=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_haul_CAN.xlsx'];
            para.bio.filename.gear_CAN=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_gear_CAN.xlsx'];
            para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_catch_CAN.xlsx'];
            para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_length_CAN.xlsx'];
            para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2021\CAN\2021_biodata_specimen_CAN_AGES.xlsx'];
            %% Stratification files
            para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2021\Stratification_geographic_Lat_2021_final.xlsx'];
            para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2021\US&CAN strata 2021_final.xlsx'];
            para.acoust.filename.processed_data_age2=['' para.data_root_dir 'Exports\US&CAN_detailsa_2021_table2y+_ALL_final.xlsx'];
            para.acoust.filename.processed_data_age1=['' para.data_root_dir 'Exports\US&CAN_detailsa_2021_table1y+_ALL_final.xlsx'];
%             para.acoust.filename.processed_data_age1=['outputs\US&CAN_detailsa_2021_table1y+_ALL_11-9-2021.xlsx'];
%             para.acoust.filename.processed_data_age2=['outputs\US&CAN_detailsa_2021_table1y+_ALL_11-9-2021.xlsx'];
            para.bio_acoust.filename.Transect_region_haul_age2=['' para.data_root_dir 'Stratification\2021\US&CAN_2021_transect_region_haul_age2+ auto_final.xlsx'];
            para.bio_acoust.filename.Transect_region_haul_age1=['' para.data_root_dir 'Stratification\2021\US&CAN_2021_transect_region_haul_age1+ auto_final.xlsx'];
            if para.proc.exclude_age1 == 1
                para.acoust.filename.processed_data=para.acoust.filename.processed_data_age2;
                para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age2;
            else
                para.acoust.filename.processed_data=para.acoust.filename.processed_data_age1;
                para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age1;
            end
            %% kriging related files
            para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2021\default_vario_krig_settings_2021_US&CAN.xlsx'];
            data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
        case 2   % US Bottom Trawl data
            fprintf('******* Bottom Trawl Data are not available *******\n\n')
        case 3   % US observer trawl data
            para.bio.filename.trawl=['' para.data_root_dir 'Observer Data\Hake_Trawl_Chu_2021.xlsx'];
            para.bio.filename.gear=[''];
            para.bio.filename.catch=['' para.data_root_dir 'Observer Data\Hake_Catch_Chu_2021.xlsx'];
            para.bio.filename.length=['' para.data_root_dir 'Observer Data\Hake_Length_Chu_2021.xlsx'];
            para.bio.filename.specimen=['' para.data_root_dir 'Observer Data\Hake_Age_Chu_2021.xlsx'];
            %% Stratification files
            para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2021\US&CAN strata 2021_final.xlsx'];
            para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2021_table2y+_ALL_final.xlsx'];
            para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2021\Stratification_geographic_Lat_2021_final.xlsx'];
            para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2021\US&CAN_2021_transect_region_haul_age2+ auto_final.xlsx'];
            %% kriging related files
            para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2021\default_vario_krig_settings_2021_US&CAN.xlsx'];
            data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
    end
end

%% kriging related files
data.in.filename.smoothed_contour=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Smoothing_EasyKrig.xlsx'];

para.proc.ST_BT_RT_ET_zero_removal_flag=0;      % 0 = not remove zeros before ST and after ET; 1 = remove zeros before ST and after ET
if para.bio_data_type ~= 1
    para.proc.stratification_index = 0;         % non- acoutical and trawl survey data only use INPFC stratification
else
para.proc.stratification_index=1;               % index for the chosen stratification
                                                % 1 = KS (trawl)-based, 2-7 = geographically based but close to trawl-based stratification
                                                % 0 = INPFC strata
                                                % 7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
                                                % 10 = one stratum for the whole survey 
end
para.proc.start_transect=1;                     % start transect number
para.proc.end_transect=122;                     % end transect number
para.proc.transect_offset=0;                    % transect offset added to the CAN transect when merge the uS and CAN data
para.proc.age1_haul=[0];                        % trawls to be excluded if age-1 is excluded
para.proc.KS_stratification=0;                  % 1 - stratification based on KS (or trawl) - based analysis
                                                % 0 - geographically defined strata
para.bio.haul_no_offset=200;                    % Canadian's trawl number offset
para.bio.CAN_strata_num0=[];                    % for combined satrta definiation file
para.bio.database_type='Oracle';                % biodata format: 'Oracle' or 'FSCS'
para.acoust.TS_station_num=2;                   % number of trawl sampling stations, whose data are used to compute the TS

                                             
