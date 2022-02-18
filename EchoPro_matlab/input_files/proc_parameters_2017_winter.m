function  proc_parameters_2017
%% this file is to provide all required input filenames & some process parameter settings
%%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/17/2015

global para data

para.proc.source=1;         % 1 = US; 2 CAN; 3 = all
if para.proc.source ~= 3
%% Biological trawl data
    para.bio.filename.trawl=['' para.data_root_dir 'Biological\2017w\2017W_biodata_haul.xlsx'];
    para.bio.filename.gear=['' para.data_root_dir 'Biological\2017w\2017W_biodata_gear.xlsx'];
    para.bio.filename.catch=['' para.data_root_dir 'Biological\2017w\2017W_biodata_catch.xlsx'];
    para.bio.filename.length=['' para.data_root_dir 'Biological\2017w\2017W_biodata_length.xlsx'];
    para.bio.filename.specimen=['' para.data_root_dir 'Biological\2017w\2017W_biodata_specimen_NOAGES.xlsx'];
%% Stratification files
    para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2017\US&CAN strata 2017 shoreside.xlsx'];
    para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2017_table1y+_ALL_draft2a.xlsx'];
    para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2017\Stratification_geographic_Lat_2017_shoreside.xlsx'];
%     para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2017\US&CAN_2017_transect_region_haul_age2+ auto draft1.xlsx'];
%     para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2017\US&CAN_2017_transect_region_haul_age1+ auto draft2.xlsx'];
    para.bio_acoust.filename.Transect_region_haul=['' para.data_root_dir 'Stratification\2017\US&CAN_2017_transect_region_haul_age1+ auto draft2a.xlsx'];
%% kriging related files
    para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2017w\default_vario_krig_settings_INPFC_final.xlsx'];
%     para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2017\default_vario_krig_settings_INPFC_10-27-2017.xlsx'];

end

%% kriging related files
data.in.filename.smoothed_contour=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Smoothing_EasyKrig.xlsx'];
data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage

para.proc.ST_BT_RT_ET_zero_removal_flag=0;      % 0 = not remove zeros before ST and after ET; 1 = remove zeros before ST and after ET
para.proc.stratification_index=1;               % index for the chosen stratification
                                                % 1 = KS (trawl)-based, 2-7 = geographically based but close to trawl-based stratification
                                                % 0 = INPFC strata
                                                % 7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
                                                % 10 = one stratum for the whole survey 
para.proc.start_transect=1;                     % start transect number
para.proc.end_transect=200;                     % end transect number
para.proc.transect_offset=0;                    % transect offset added to the CAN transect when merge the uS and CAN data
para.proc.age1_haul=[0];                        % trawls to be excluded if age-1 is excluded
para.proc.KS_stratification=0;                  % 1 - stratification based on KS (or trawl) - based analysis
                                                % 0 - geographically defined strata
para.bio.haul_no_offset=200;                    % Canadian's trawl number offset
para.bio.CAN_strata_num0=[];                    % for combined satrta definiation file
para.bio.database_type='Oracle';                % biodata format: 'Oracle' or 'FSCS'
para.acoust.TS_station_num=2;                   % number of trawl sampling stations, whose data are used to compute the TS

                                             
