function construct_catch_table
% construct a catch table
% [trawl_no catch expanded_number expanded_weight];
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global para data

k=0;
for i=1:length(data.bio.catch)
    for j=1:length(data.bio.catch(i).species)
        if data.bio.catch(i).species(j).ID == para.bio.species_code_ID
            k=k+1;
            data.final.table.catch(k,1)=data.bio.catch(i).trawl_no;
            data.final.table.catch(k,2)=data.bio.catch(i).species(j).exp_cnt;
            data.final.table.catch(k,3)=data.bio.catch(i).species(j).exp_wgt;
        end
    end
end

return