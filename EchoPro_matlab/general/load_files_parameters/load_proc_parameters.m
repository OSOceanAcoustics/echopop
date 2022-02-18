function load_proc_parameters(survey_year)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para

% construct parameter filename
if isfield(para, 'platform_name') 
    if strcmp(para.platform_name, 'SD')  % SD data
%%     load SD datafile names
        para_setting_filename=['proc_parameters_' survey_year '_SD'];
    else   % SFV
        para_setting_filename=['proc_parameters_' survey_year];
    end
else
%%     load FSV datafile names
    para_setting_filename=['proc_parameters_' survey_year];
end

% para_setting_filename=['proc_parameters_' survey_year];

eval(para_setting_filename);
return