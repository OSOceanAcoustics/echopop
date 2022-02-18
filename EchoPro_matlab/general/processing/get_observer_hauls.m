function hauls = get_observer_hauls(VL_BiomassInt, stratification_filename)
%% find corresponding haul for each cell based on geographically stratified strata for observer data
% Dezhang Chu
% 10-25-2018

global para data

lat = VL_BiomassInt(:,5);
n = length(lat);

%% read geographic stratification file
dat=xlsread(stratification_filename,para.proc.stratification_index+1);

for i=1:n
    dif = lat(i) - dat(:,2);
    ind = find(dif > 0);
    if isempty(ind)
        ind = 0;
    else
        switch length(ind)
            case 1
               m=length(ind);
            case 2
               m=length(ind);
            case 3
               m=length(ind);
            case 4
               m=length(ind);
        end
    end
    hauls(i) = data.bio.strata(ind(end) + 1).trawls(1);    % only the first trawl since it will fall in the same stratum
end

end