function  proc_parameters_2001
% initialization operation to read all default psettings and files

global para data hdl

%% Biological trawl data
para.bio.filename.trawl_US=['' para.data_root_dir 'Biological\2001\US\biodata_haul.xls'];
para.bio.filename.gear_US=['' para.data_root_dir 'Biological\2001\US\biodata_gear.xls'];
para.bio.filename.catch_US=['' para.data_root_dir 'Biological\2001\US\biodata_catch.xls'];
para.bio.filename.length_US=['' para.data_root_dir 'Biological\2001\US\biodata_length.xls'];
para.bio.filename.specimen_US=['' para.data_root_dir 'Biological\2001\US\biodata_specimen.xls'];
para.bio.filename.trawl_CAN='';
para.bio.filename.gear_CAN='';
para.bio.filename.catch_CAN=['' para.data_root_dir 'Biological\2001\CAN\biodata_catch.xls'];
para.bio.filename.length_CAN=['' para.data_root_dir 'Biological\2001\CAN\biodata_length.xls'];
para.bio.filename.specimen_CAN=['' para.data_root_dir 'Biological\2001\CAN\biodata_specimen.xls'];

%% Stratification files
para.bio_acoust.filename.Transect_region_haul='';
para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2001\US&CAN strata 2001.xlsx'];
para.proc.stratification_filename=['' para.data_root_dir 'Stratification\2001\Stratification_geographic_Lat_2001.xlsx'];

%% NASC data
para.proc.transect_info_filename=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Transect Bounds to 2011.xlsx'];                           % ST, BT, RT,and ET of transect for removing extra zeros 
para.acoust.filename.processed_data=['' para.data_root_dir 'Exports\US&CAN_detailsa_2001_table2y+_ALL_final.xlsx'];

%% kriging related files
data.in.filename.smoothed_contour=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\Smoothing_EasyKrig.xlsx'];
data.in.filename.grid_cell=['' para.data_root_dir 'Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_2013.xlsx'];                         % 2013 cell res = 2.50 nmi with extended area coverage
para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2001\default_vario_krig_settings_final.xlsx'];
% para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2001\default_vario_krig_settings_orig.xlsx'];

% old parameter files
% para.acoust.filename.strata=['' para.data_root_dir 'Stratification\2009\US&CAN strata 2009_old.xlsx'];
% para.krig.vario_krig_para_filename=['' para.data_root_dir 'Kriging files & parameters\2009\default_vario_krig_settings_orig.xlsx'];


para.proc.ST_BT_RT_ET_zero_removal_flag=1;      % 0 = not remove zeros before ST and after ET; 1 = remove zeros before ST and after ET
para.proc.stratification_index=1;               % index for the chosen stratification
                                                % 1 = KS (trawl)-based, 2-3 = geographically based but close to trawl-based stratification
                                                % 0 = INPFC strata
para.proc.start_transect=1;                     % start transect number
para.proc.CAN_start_transect=0;                % start transect number for transects on Canadian Ship
para.proc.end_transect=200;                     % end transect number
para.proc.transect_offset=0;                    % transect offset added to the CAN transect when merge the uS and CAN data
para.proc.age1_haul=[ ];                        % trawls to be excluded if age-1 is excluded
para.proc.KS_stratification=1;                  % 1 - stratification based on KS (or trawl) - based analysis
                                                % 0 - geographically defined strata
para.bio.haul_no_offset=200;                    % Canadian's trawl number offset
para.bio.CAN_strata_num0=[];                    % for combined satrta definiation file
para.bio.database_type='Oracle';                % biodata format: 'Oracle' or 'FSCS'
para.acoust.TS_station_num=1;                   % number of trawl sampling stations, whose data are used to compute the TS
