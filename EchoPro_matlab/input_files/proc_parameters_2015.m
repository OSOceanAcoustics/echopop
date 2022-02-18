function  proc_parameters_2015
%% this file is to provide all required input filenames & some process parameter settings
%%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/17/2015

global para data

% para.proc.source=3;
if para.proc.source ~= 3
%% Biological trawl data
    para.bio.filename.trawl=['' para.data_root_dir 'Biological\2015\at-sea US\2015_biodata_haul.xlsx'];
    para.bio.filename.gear=['' para.data_root_dir 'Biological\2015\at-sea US\2015_biodata_gear.xlsx'];
    para.bio.filename.catch=['' para.data_root_dir 'Biological\2015\at-sea US\2015_biodata_catch.xlsx'];
    para.bio.filename.length=['' para.data_root_dir 'Biological\2015\at-sea US\2015_biodata_length.xlsx'];
    para.bio.filename.specimen=['' para.data_root_dir 'Biological\2015\at-sea US\2015_biodata_specimen.xlsx'];
%% Stratification files
    para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2015\at-sea\US strata 2015 atseav2.xlsx'];
    para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\2015\US_atsea_detailsa_2015_table2y+_ALL_2015_10_14.xlsx'];
    para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2015\at-sea\US_atsea_2015_transect_region_haul.xlsx'];
    para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2015\at-sea\Stratification_geographic_Lat_2015_atsea.xlsx'];
%% kriging related files
    para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2015\at-sea US\default_vario_krig_settings.xlsx'];
elseif para.proc.source == 2   % CAN only
    para.bio.filename.trawl_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_haul_CAN.xlsx'];
    para.bio.filename.gear_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_gear_CAN.xlsx'];
    para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_catch_CAN.xlsx'];
    para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_length_CAN.xlsx'];
    para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_specimen_CAN.xlsx'];
    fprintf('******* No CAN stratification files *******\n')
    fprintf('******* No CAN Transsect-region-haul files *******\n')
    fprintf('******* No CAN NASC export data files *******\n')
    fprintf('******* No CAN vario-krig parameter file *******\n')
    fprintf('******* No CAN grid cell file *******\n\n')
elseif para.proc.source == 3 % US & CAN acoustic trawl data
    switch para.bio_data_type
        case 1   % Acoustic and Trawl survey data
            %% Biological trawl data
            para.bio.filename.trawl_US=['' para.data_root_dir 'Biological\2015\aged shoreside US\2015_biodata_haul.xlsx'];
            para.bio.filename.gear_US=['' para.data_root_dir 'Biological\2015\aged shoreside US\2015_biodata_gear.xlsx'];
            para.bio.filename.catch_US=['' para.data_root_dir 'Biological\2015\aged shoreside US\2015_biodata_catch.xlsx'];
            para.bio.filename.length_US=['' para.data_root_dir 'Biological\2015\aged shoreside US\2015_biodata_length.xlsx'];
            para.bio.filename.specimen_US=['' para.data_root_dir 'Biological\2015\aged shoreside US\2015_biodata_specimen.xlsx'];
            para.bio.filename.trawl_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_haul_CAN.xlsx'];
            para.bio.filename.gear_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_gear_CAN.xlsx'];
            para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_catch_CAN.xlsx'];
            para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_length_CAN.xlsx'];
            para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2015\aged shoreside CAN\2015_biodata_specimen_CAN.xlsx'];
            %% Stratification files
            para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\US&CAN strata 2015 shoreside.xlsx'];
            para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\US&CAN_snapshot1_2015_transect_region_haul_age2+ auto final.xlsx'];
            para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\Stratification_geographic_Lat_2015_shoreside.xlsx'];
%             para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_shoreside_snapshot1_detailsa_2015_table2y+_ALL_final_2016-02-09.xlsx'];
            para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2015_table2y+_shoreside_snapshot1_ALL_final.xlsx'];   % filename changed on 1/19/2021
            %% kriging related files
            para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2015\shoreside_snapshot-1 US&CAN\default_vario_krig_settings_final.xlsx'];
            para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2015\shoreside_snapshot-1 US&CAN\default_vario_krig_settings_KS_final.xlsx'];
        case 2   % US Bottom Trawl data
            fprintf('******* Bottom Trawl Data are not available *******\n\n')
        case 3   % US observer trawl data
            para.bio.filename.trawl=['' para.data_root_dir 'Observer Data\Hake_Trawl_Chu_2015.xlsx'];
            para.bio.filename.gear=[''];
            para.bio.filename.catch=['' para.data_root_dir 'Observer Data\Hake_Catch_Chu_2015.xlsx'];
            para.bio.filename.length=['' para.data_root_dir 'Observer Data\Hake_Length_Chu_2015.xlsx'];
            para.bio.filename.specimen=['' para.data_root_dir 'Observer Data\Hake_Age_Chu_2015.xlsx'];
            %% Stratification files
            para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\US&CAN strata 2015 shoreside.xlsx'];
            para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\US&CAN_snapshot1_2015_transect_region_haul_age2+ auto final.xlsx'];
            para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2015\shoreside_snapshot-1\Stratification_geographic_Lat_2015_shoreside.xlsx'];
%             para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_shoreside_snapshot1_detailsa_2015_table2y+_ALL_final_2016-02-09.xlsx'];
            para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2015_table2y+_shoreside_snapshot1_ALL_final.xlsx'];   % filename changed on 1/19/2021
            %% kriging related files
            para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2015\shoreside_snapshot-1 US&CAN\default_vario_krig_settings_final.xlsx'];
    end
end

% %% Stratification files
% para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2015\at-sea\Stratification_geographic_Lat_2015_atsea.xlsx'];

%% NASC data
para.proc.transect_info_filename=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Transect Bounds from 2015.xlsx'];               % ST, BT, RT,and ET of transect for removing extra zeros 

%% kriging related files
data.in.filename.smoothed_contour=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Smoothing_EasyKrig.xlsx'];
data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
% data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid1_25nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage

para.proc.ST_BT_RT_ET_zero_removal_flag=0;      % 0 = not remove zeros before ST and after ET; 1 = remove zeros before ST and after ET
para.proc.stratification_index=1;               % index for the chosen stratification
                                                % 1 = KS (trawl)-based, 2-7 = geographically based but close to trawl-based stratification
                                                % 0 = INPFC strata
                                                % 7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
                                                % 10 = one stratum for the whole survey 
para.proc.start_transect=1;                     % start transect number
para.proc.end_transect=1120;                     % end transect number
para.proc.transect_offset=0;                    % transect offset added to the CAN transect when merge the uS and CAN data
para.proc.age1_haul=[0];                        % trawls to be excluded if age-1 is excluded
para.proc.KS_stratification=1;                  % 1 - stratification based on KS (or trawl) - based analysis
                                                % 0 - geographically defined strata
para.bio.haul_no_offset=200;                    % Canadian's trawl number offset
para.bio.CAN_strata_num0=[];                    % for combined satrta definiation file
para.bio.database_type='Oracle';                % biodata format: 'Oracle' or 'FSCS'
para.acoust.TS_station_num=2;                   % number of trawl sampling stations, whose data are used to compute the TS

                                             
