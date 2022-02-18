function    proc_acoustic_data
% process the acoustic echogram data to get biomass estimates
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global para data

if para.acoust.file_type == 1
    % process processed acoustic data (NASC data) exported from EchoView
    proc_NASC_data
else
    % processed echogram data --> echo-integrated NASC data
    proc_acoustic_raw_data
end

