function      stratum_id=get_geographic_strata(lat,strata_def_filename,stratum_indx)
% the function get_geographic_strata will return the stratum_id based on
% the strata definition specified geographically in file whose name is
% provided
% INPUTS:
%    lat = latitude array
%    strata_def_filename = filename where the stratification is provided
%    stratum_indx = index of the stratification method
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       11/29/2013

global para data

if para.proc.stratification_index == 10
    dat=xlsread(strata_def_filename,'length strata byhaul_1stratum');
else
    dat=xlsread(strata_def_filename,para.proc.stratification_index+1);
end
stratum_id=ones(size(lat));

%% number of strata
n=size(dat,1);
%% added the following two lines on 12-18-2019
%% ----------------------------------------------------
ind=find(lat <= dat(1,2));
stratum_id(ind)=dat(1,1);
%% ----------------------------------------------------

for i=2:n
    if dat(i-1,2) < dat(i,2)
       ind=find(lat >= dat(i-1,2) & lat <= dat(i,2));
    elseif i > 2
       ind=find(lat >= dat(i-2,2) & lat <= dat(i,2));
    end
    stratum_id(ind)=dat(i,1);
end
% % ind = find(stratum_id > length(data.bio.strata));
% % if ~isempty(ind)
% %     stratum_id(ind) = length(data.bio.strata);
% % end
% %% replaced the above 4 lines with the following 16 lines on 12-18-2019
% %% -------------------------------------------------------------------------------
% ind=find(lat > dat(end,2));
% stratum_id(ind)=dat(end,1)+1;
% 
% %% fill-in strata without trawls that are outside the strata with those of the cloest tratum that has trawls
% n_min = min(data.bio.haul_strata);
% if n_min > 1
%     for i = 1:n_min - 1
%         data.bio.strata(i) = data.bio.strata(n_min);
%     end
% end
% n_max = max(data.bio.haul_strata);
% if n_max <= size(dat,1)+1
%     for i = max(data.bio.haul_strata)+1:dat(end,1)+1
%         data.bio.strata(i) = data.bio.strata(n_max);
%     end
% end
% %% -----------------------------------------------------------------------------

%% add the following 16 lines on 12-19-2019
%% -------------------------------------------------------------------------------
if para.bio_data_type == 3 | para.proc.stratification_index == 0
    ind=find(lat > dat(end,2));
%     stratum_id(ind)=dat(end,1)+1;
    stratum_id(ind)=dat(end,1);
    
    % fill-in strata without trawls that are outside the strata with those of the cloest tratum that has trawls
    n_min = min(data.bio.haul_strata);
    if n_min > 1
        for i = 1:n_min - 1
            data.bio.strata(i) = data.bio.strata(n_min);
        end
    end
    n_max = max(data.bio.haul_strata);
    if n_max <= size(dat,1)+1
        for i = max(data.bio.haul_strata)+1:dat(end,1)+1
            data.bio.strata(i) = data.bio.strata(n_max);
        end
    end
end
%% -----------------------------------------------------------------------------

return