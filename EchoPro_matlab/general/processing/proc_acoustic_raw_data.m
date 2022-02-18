function proc_acoustic_raw_data
%% process raw acoustic echogram data, either EK60 or EK500
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global hdl para data

if para.acoust.file_type ~= 1
%% load transect-haul file
    out=load_transect_trawl_info(para.bio.filename.transect_haul);
    date.bio.TX=out.TX;
end

for i=para.proc.start_transect:para.proc.end_transect
    
    load_transect_files           % VL log, evr file, EV bot file    
end