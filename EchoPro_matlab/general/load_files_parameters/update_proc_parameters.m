function   update_proc_parameters(hObject)
%% call appropriate functions to update process settings based on the settings in 
%% 'Reload Parameter Files' and/or 'Processing' windows
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global hdl para data

if nargin < 1
% update all parameters based on the new settings 
else
    switch get(hObject,'Tag')
        case 'pb_load_update'
            update_filename_window_settings
        case 'pb_processing'
            update_processing_window_settings
    end
end


return