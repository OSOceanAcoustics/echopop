function  proc_parameters_2012
%% this file is to provide all required input filenames & some process parameter settings
%%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       7/17/2015

global para data
% para.home_dir='C:\Projects\EchoPro\EchoProGUI_2012';
% para.proc.source=3;
if para.proc.source == 1   % US only
    para.bio.filename.trawl_US=['' para.data_root_dir 'Biological\2012\US\biodata_haul.xlsx'];
    para.bio.filename.gear_US=['' para.data_root_dir 'Biological\2012\US\biodata_gear.xlsx'];
    para.bio.filename.catch_US=['' para.data_root_dir 'Biological\2012\US\biodata_catch.xlsx'];
    para.bio.filename.length_US=['' para.data_root_dir 'Biological\2012\US\biodata_length.xlsx'];
    para.bio.filename.specimen_US=['' para.data_root_dir 'Biological\2012\US\biodata_specimen.xlsx'];
    fprintf('******* No US stratification files *******\n')
    fprintf('******* No US Transsect-region-haul files *******\n')
    fprintf('******* No US NASC export data files *******\n')
    fprintf('******* No US vario-krig parameter file *******\n')
    fprintf('******* No US grid cell file *******\n\n')
elseif para.proc.source == 2   % CAN only
    para.bio.filename.trawl_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_haul_Canada.xlsx'];
    para.bio.filename.gear_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_gear_Canada.xlsx'];
    para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_catch_Canada.xlsx'];
    para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_length_Canada.xlsx'];
    para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_specimen_Canada.xlsx'];
    fprintf('******* No CAN stratification files *******\n')
    fprintf('******* No CAN Transsect-region-haul files *******\n')
    fprintf('******* No CAN NASC export data files *******\n')
    fprintf('******* No CAN vario-krig parameter file *******\n')
    fprintf('******* No CAN grid cell file *******\n\n')
elseif para.proc.source == 3 % US & CAN acoustic trawl data
    switch para.bio_data_type
        case 1   % Acoustic and Trawl survey data
            %% Biological trawl data
            para.bio.filename.trawl_US=['' para.data_root_dir 'Biological\2012\US\biodata_haul.xlsx'];
            para.bio.filename.gear_US=['' para.data_root_dir 'Biological\2012\US\biodata_gear.xlsx'];
            para.bio.filename.catch_US=['' para.data_root_dir 'Biological\2012\US\biodata_catch.xlsx'];
            para.bio.filename.length_US=['' para.data_root_dir 'Biological\2012\US\biodata_length.xlsx'];
            para.bio.filename.specimen_US=['' para.data_root_dir 'Biological\2012\US\biodata_specimen.xlsx'];
            para.bio.filename.trawl_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_haul_Canada.xlsx'];
            para.bio.filename.gear_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_gear_Canada.xlsx'];
            para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_catch_Canada.xlsx'];
            para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_length_Canada.xlsx'];
            para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2012\CAN\biodata_specimen_Canada.xlsx'];
        case 2   % US Bottom Trawl data
            fprintf('******* Bottom Trawl Data are not available *******\n\n')
        case 3   % US observer trawl data
            para.bio.filename.trawl=['' para.data_root_dir 'Observer Data\Hake_Trawl_Chu_2012.xlsx'];
            para.bio.filename.gear=[''];
            para.bio.filename.catch=['' para.data_root_dir 'Observer Data\Hake_Catch_Chu_2012.xlsx'];
            para.bio.filename.length=['' para.data_root_dir 'Observer Data\Hake_Length_Chu_2012.xlsx'];
            para.bio.filename.specimen=['' para.data_root_dir 'Observer Data\Hake_Age_Chu_2012.xlsx'];
    end
end
%% Stratification files
para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2012\US&CAN_T_reg_haul_FINAL-12-13-2012.xlsx'];
para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2012\US&CAN strata 2012 11-30-2012.xlsx'];
para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2012\Stratification_geographic_Lat.xlsx'];

%% NASC data
para.proc.transect_info_filename=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Transect Bounds from 2012.xlsx'];               % ST, BT, RT,and ET of transect for removing extra zeros 
para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2012_table2y+_ALL_final.xlsx'];
%% compare 1500m with 35 nm west bound
% para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2012_table2y+_ALL_final - 1500m.xlsx'];

%% kriging related files
data.in.filename.smoothed_contour1=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Smoothing_EasyKrig.xlsx'];
data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2012\default_vario_krig_settings_final.xlsx'];
% para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2012\default_vario_krig_settings_orig.xlsx'];

%% old program
% data.in.filename.grid_cell=['' para.home_dir '\input_data\2012\krig_grid2_5nm_cut_centroids_2012_11Jan12.xlsx'];
% data.in.filename.smoothed_contour=['' para.home_dir '\kriging\mesh_cells\Smoothing_EasyKrig.xlsx'];

para.proc.ST_BT_RT_ET_zero_removal_flag=1;      % 0 = not remove zeros before ST and after ET; 1 = remove zeros before ST and after ET
para.proc.stratification_index=1;               % index for the chosen stratification
                                                % 1 = KS (trawl)-based, 2-7 = geographically based but close to trawl-based stratification
                                                % 0 = INPFC strata
                                                % 7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
para.proc.start_transect=1;                     % start transect number
para.proc.end_transect=144;                     % end transect number
para.proc.transect_offset=0;                    % transect offset added to the CAN transect when merge the uS and CAN data
para.proc.age1_haul=[0];                        % trawls to be excluded if age-1 is excluded
para.proc.KS_stratification=1;                  % 1 - stratification based on KS (or trawl) - based analysis
                                                % 0 - geographically defined strata
para.bio.haul_no_offset=100;                    % Canadian's trawl number offset
para.bio.CAN_strata_num0=[];                    % for combined satrta definiation file
para.bio.database_type='Oracle';                % biodata format: 'Oracle' or 'FSCS'
para.acoust.TS_station_num=2;                   % number of trawl sampling stations, whose data are used to compute the TS



