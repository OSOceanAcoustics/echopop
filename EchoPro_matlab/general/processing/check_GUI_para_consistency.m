function para=check_GUI_para_consistency(hdl,para)
% check whether the processing parameter (para.proc) are consistent with the corresponding 
% variables set on the "Processing" GUI window 

% extrapolation
if get(hdl.proc.radio_no_extrapolation,'value') == 1
    para.proc.extrapolation=0;
else
    para.proc.extrapolation=1;
end

% age data
if get(hdl.proc.radio_agedata,'value') == 1
    para.proc.age_data_status=1;
else
    para.proc.age_data_status=2;
end

%% biological data type
para.bio_data_type = get(hdl.proc.popup_trawl_data,'value');

% stratification
switch get(hdl.proc.popup_proc_stratification,'value')
    case 1
        para.proc.stratification_index=1;
        para.proc.KS_stratification=1;
    case 2
        para.proc.stratification_index=0;
        para.proc.KS_stratification=0;
    case 3
        para.proc.stratification_index=2;
        para.proc.KS_stratification=1;
    otherwise
        disp('not a valid option!')
end
     
% Biomass
para.proc.kriging_input=get(hdl.proc.popup_proc_variable,'value');

% Transect & Bootstrapping parameters
para.proc.start_transect=str2num( get(hdl.proc.edit_start_transect_num,'string'));                 % start transect number
para.proc.end_transect=str2num( get(hdl.proc.edit_end_transect_num,'string'));                     % end transect number
para.proc.bootstrap.limit=str2num( get(hdl.proc.edit_bootstrap_limit,'string'));
para.proc.transect_reduction_fraction=num2str(get(hdl.proc.edit_transect_reduction_fraction,'string'));          % spacing reduction factor

% kriging setup
para.proc.kriging=get(hdl.proc.radio_kriging_proc,'value');                         % 1= kriging processing, 0 = no kriging
para.proc.default_para=get(hdl.proc.radio_default_para,'value');                 % 1 = default variogram & kriging parameter file; 0 = manually determines variogram & kriging parameters
