function construct_observer_catch_trawl_output_matrices
% construct the hake (or taget species) catch output matrix for
% visualization using observer data
% & generate report tables
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para data

%% Catch Table
%% Columns   1-3:   'Trawl number'   'Abundance'     'Biomass'          

n_tot_trawl=length(data.bio.catch);  
hake_trawl_no=[];
k=0;
for i=1:n_tot_trawl
    for j=1:length(data.bio.catch(i).species)
        if para.bio.species_code_ID == data.bio.catch(i).species(j).ID
            k=k+1;
            data.final.table.catch(k,1)=data.bio.catch(i).trawl_no;
            data.final.table.catch(k,2)=data.bio.catch(i).species(j).exp_cnt;   % number
            data.final.table.catch(k,3)=data.bio.catch(i).species(j).exp_wgt;   % weight in kg           
        end
    end
end
            
n_hake_len_gender_trawl=length(data.bio.hake_length_sex);
n_hake_len_wgt_gender_age_trawl=length(data.bio.hake_length_weight_sex_age);
acc_k=0;
%% Trawl Table
%% Columns   1-6:   'Trawl number'   'Transect Number'   'EQ Lat'   'EQ Lon'  'Stratum  id'  'Depth' 
%% Columns   7-13:  'surface temperature' ' gear depth temperature'  'hake length'   'hake gender'  'hake age'   'hake weight' 'Gear depth'

% dat=xlsread(para.acoust.filename.strata);

if para.bio_data_type == 1 % acoustical survey trawls
    if para.proc.stratification_index == 10
        dat=xlsread(para.acoust.filename.strata,'length strata byhaul_1stratum');
    else
        dat=xlsread(para.acoust.filename.strata,para.proc.stratification_index+1);
    end
    haul_strata=[dat(:,4) dat(:,2)];
elseif para.bio_data_type == 3  % observer
    % only allows pre-stratification
    dat = xlsread(para.proc.stratification_filename,1);
    %% construct a haul_strata array that is consistent with the acoustic trawl format
    n_strata = size(dat,1);
    haul_strata = [];
    for i = 1:n_strata
        if i == 1
            ind = find(data.bio.trawl.EQlat <= dat(i,2));
        elseif i == n_strata
            ind = find(data.bio.trawl.EQlat > dat(i,2));
        else
            ind = find(data.bio.trawl.EQlat > dat(i,2) & data.bio.trawl.EQlat <= dat(i+1,2));
        end
        haul_strata = [haul_strata; dat(i)*ones(length(ind),1)  data.bio.trawl.trawl_no(ind)];
    end
end

data.bio.haul_strata = haul_strata;

for i=1:length(data.bio.trawl.trawl_no)
    %% length-sex samples (station #1)
    n_len_gender=[];
    for j=1:n_hake_len_gender_trawl
        if data.bio.trawl.trawl_no(i) == data.bio.hake_length_sex(j).trawl_no
            n_len_gender=length(data.bio.hake_length_sex(j).length);
            break
        end
    end
    catch_weight(i)=0;
    for k1=1:n_tot_trawl
        if data.bio.trawl.trawl_no(i) == data.bio.catch(k1).trawl_no 
            for k2=1:length(data.bio.catch(k1).species)
              if data.bio.catch(k1).species(k2).ID ==  para.bio.species_code_ID
                  catch_weight(i)=data.bio.catch(k1).species(k2).exp_wgt;
%                   fprintf('trawl_seq = %d\t trawl no = %d\t unaged weight = %7.3f\n',i,data.bio.catch(k1).trawl_no,catch_weight(i));
                 break
              end
            end
        end
    end
%     data.final.table.trawl(1,5)=0;
    if ~isempty(n_len_gender)
        indx1=acc_k+(1:n_len_gender);
        data.final.table.trawl(indx1,1)=data.bio.trawl.trawl_no(i);    % trawl number
        if ~isempty(para.bio.filename.gear)
            data.final.table.trawl(indx1,2)=data.bio.gear.transect(i);    % transect number
        end
        %% corresponding index in the trawl structure
        trawl_ind=find(data.bio.trawl.trawl_no == data.bio.trawl.trawl_no(i));
        data.final.table.trawl(indx1,3)=data.bio.trawl.EQlat(trawl_ind);      % EQ Lat
        data.final.table.trawl(indx1,4)=data.bio.trawl.EQlon(trawl_ind);      % EQ Lon
        indx_stratum=find(haul_strata(:,1) == data.bio.trawl.trawl_no(i));
        if ~isempty(indx_stratum)
            data.final.table.trawl(indx1,5)=haul_strata(indx_stratum,2);    % Stratum
        else
            data.final.table.trawl(indx1,5)=nan;
            if para.proc.bio_info_disp == 1
                fprintf('Undefined stratum for trawl_no = %d !!!\n',data.bio.trawl_trawl_no(i));
            end
        end
        data.final.table.trawl(indx1,6)=data.bio.trawl.ave_dep(trawl_ind);    % Depth
        if ~isempty(para.bio.filename.gear)
            data.final.table.trawl(indx1,7)=data.bio.gear.Tsurf(i);       % surface temperature
            data.final.table.trawl(indx1,8)=data.bio.gear.Tdep(i);        % gear depth temperature
            data.final.table.trawl(indx1,13)=data.bio.gear.footropeD(i)-10;    % Gear Depth (assume 20 m openning)
        end
        data.final.table.trawl(indx1,9)=data.bio.hake_length_sex(j).length(:);       % hake length
        data.final.table.trawl(indx1,10)=data.bio.hake_length_sex(j).Gender(:);       % hake gender
        data.final.table.trawl(indx1,11)=nan;                                         % hake age
        data.final.table.trawl(indx1,12)=nan;                                         % hake weight (aged)
        data.final.table.trawl(indx1,14)=catch_weight(i);                                % hake weight (un-aged)
        acc_k=indx1(end);
    end
%     ind=find(data.final.table.trawl(:,5) == 0);
%     if ~isempty(ind) & para.proc.bio_info_disp == 1
%         fprintf('not empty: %d\n',length(ind))
%     end
    %% length-weight-age-sex samples (station #2)
    n_len_wgt_gender_age=[];
    for j=1:n_hake_len_wgt_gender_age_trawl
        if data.bio.trawl.trawl_no(i) == data.bio.hake_length_weight_sex_age(j).trawl_no
            n_len_wgt_gender_age=length(data.bio.hake_length_weight_sex_age(j).length);
            break
        end
    end
    if ~isempty(n_len_wgt_gender_age) & n_len_wgt_gender_age ~= 0
        indx2=acc_k+(1:n_len_wgt_gender_age);
        data.final.table.trawl(indx2,1)=data.bio.trawl.trawl_no(i);    % trawl number
        if ~isempty(para.bio.filename.gear)
            data.final.table.trawl(indx2,2)=data.bio.gear.transect(i);    % transect number
        end
        %% corresponding index in the trawl structure
        trawl_ind=find(data.bio.trawl.trawl_no == data.bio.trawl.trawl_no(i));
        data.final.table.trawl(indx2,3)=data.bio.trawl.EQlat(trawl_ind);      % EQ Lat
        data.final.table.trawl(indx2,4)=data.bio.trawl.EQlon(trawl_ind);      % EQ Lon
        indx_stratum=find(haul_strata(:,1) == data.bio.trawl.trawl_no(i));
        if ~isempty(indx_stratum)
            data.final.table.trawl(indx2,5)=haul_strata(indx_stratum,2);    % Stratum
        else
            data.final.table.trawl(indx2,5)=nan;
        end
        data.final.table.trawl(indx2,6)=data.bio.trawl.ave_dep(trawl_ind);    % Depth
        if ~isempty(para.bio.filename.gear)
            data.final.table.trawl(indx2,7)=data.bio.gear.Tsurf(i);       % surface temperature
            data.final.table.trawl(indx2,8)=data.bio.gear.Tdep(i);        % gear depth temperature
            data.final.table.trawl(indx2,13)=data.bio.gear.footropeD(i)-10;    % Gear Depth (assume 20 m openning)
        end
        data.final.table.trawl(indx2,9)=data.bio.hake_length_weight_sex_age(j).length(:);       % hake length
        data.final.table.trawl(indx2,10)=data.bio.hake_length_weight_sex_age(j).sex(:);      % hake gender
        data.final.table.trawl(indx2,11)=data.bio.hake_length_weight_sex_age(j).age(:);         % hake age
        data.final.table.trawl(indx2,12)=data.bio.hake_length_weight_sex_age(j).wgt(:);         % hake weight (aged)
        acc_k=indx2(end);
    end
%     if ~isempty(ind) & para.proc.bio_info_disp == 1
%         fprintf('not empty: %d\n',length(ind))
%     end
end  
    
return