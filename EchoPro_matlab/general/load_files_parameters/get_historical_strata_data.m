function      get_historical_strata_data
% get quantities and keys associated with strata for:
% 1. trawl information
% 2. Length - key 
% 3. Length-weight key 
% 4. Length-age - key for abundance & biomass
% 5. etc
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Modification:            3/24/2013
%% Modification:       12/09/2020  to change the excel file reading lines 47 & 50, and commented out lines 23-30, and line 44
%% Modification:       05/20/2021  to change the extroplated strata from max(dat0(:,1)) to max(dat0(:,2)) on line 152
%% Modification:       12/11/2021  to change the extroplated strata from max(dat0(:,1)) if-loop and added line 153-160
%% Modification:       12/14/2021  to fill-in parameters in the strata where no fish samples (age2+ case) in there by using the parameters of the adjacient strata


global data para

t0 =clock;

fid=fopen(para.acoust.filename.strata);
if fid < 0
    return
end

% i=0;
% while ~feof(fid)
%    i=i+1;
%    line=str2num(fgetl(fid));
% %    line=fgetl(fid);
% %    fprintf('%d\t %s\n',i, line);
% end
% fclose(fid);

%%% read haul-strata data file
data.bio.strata=[];
data.bio.len_haul_M=[];
data.bio.len_haul_F=[];
data.bio.len_haul_ALL=[];
data.bio.aged_len_haul_M=[];
data.bio.aged_len_haul_F=[];
data.bio.aged_len_haul_ALL=[];
data.bio.compact_len_haul_M=[];
data.bio.compact_len_haul_F=[];
data.bio.compact_len_haul_ALL=[];

% file_length=i;
if para.proc.stratification_index == 10
%     dat0 = xlsread(para.acoust.filename.strata,'length strata byhaul_1stratum',['A2:E' num2str(file_length)]);
    dat0 = xlsread(para.acoust.filename.strata,'length strata byhaul_1stratum');   % changed on 12/09/2020
else
%     dat0 = xlsread(para.acoust.filename.strata,para.proc.stratification_index+1,['A2:E' num2str(file_length)]);
    dat0 = xlsread(para.acoust.filename.strata,para.proc.stratification_index+1);  % changed on 12/09/2020
end

if para.proc.stratification_index == 0 & para.proc.KS_stratification == 0 % INPFC with 
    if dat0(1,1) == 2015
        %% remove age1 hauls but not those used in transect_region_hual files, to be consistent with K-S stratification
        para.proc.age1_haul = [4 9 10 13 15 18 23 30 35 39 41 48 60 69 227];
    end
end

% [junk, strata_name] = xlsread(para.acoust.filename.strata,'',['C2:C' num2str(file_length)]);
% indx=find(str2num(para.survey_year ) == dat0(:,1));
% dat=dat0(indx,[2 4 5]);                         % strata-haul-weight factor

if para.bio_data_type ~= 3
    indx = find(str2num(para.survey_year ) == dat0(:,1));
    dat = dat0(indx,[2 4 5]);                         % strata-haul-weight factor    
    stratum_id = dat(:,1);
% %     %% 2019 OR & WA trawls
% %     included_hauls = [42   43   44   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64   65   66   67   69   70   71];
% %     ind = find(dat(:,2) >= included_hauls(1) & dat(:,2) <= included_hauls(end));
% %     dat = dat(ind,:);
% %     stratum_id = stratum_id(ind);
elseif para.bio_data_type == 3  % observer data
    dat0 = xlsread(para.proc.stratification_filename,para.proc.stratification_index+1);    
    n=size(dat0,1);    
    lat = data.bio.trawl.EQlat;
    ind = find(lat <= dat0(1,2));
    stratum_id(ind) = dat0(1,1);
    for i=2:n
        ind=find(lat > dat0(i-1,2) & lat <= dat0(i,2));
        stratum_id(ind)=dat0(i,1);
    end
    ind = find(lat > dat0(end,2));
    stratum_id(ind) = dat0(end,1) + 1;
    
    dat = [stratum_id(:) data.bio.trawl.trawl_no(:) ones(length(data.bio.trawl.trawl_no),1)];
end

data.bio_acoust.haul_wgt_tbl=dat(:,2:3);      
data.bio.haul_strata=unique(dat(:,1));
%% remove strata index 0 -- YOY
ind = find(data.bio.haul_strata == 0);
data.bio.haul_strata(ind)=[];
n=length(data.bio.haul_strata);                        % number of strata

if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
        fprintf('----- gether trawl info first for all strata ...\n')
    end
end

% n=data.bio.haul_strata(end);
for i=1:n
    if data.bio.haul_strata(i) ~= 0                    % not process stratum 0
        indx=find(dat(:,1) == data.bio.haul_strata(i));
        [uniques_trawls indx1]=unique(dat(indx,2));
        if ~isempty(para.bio.CAN_strata_num0) & data.bio.haul_strata(i) >= para.bio.CAN_strata_num0
            %% combining CAN hauls and add haul number offset is haul numbers are
            %% overlapping with those from US
            uniques_trawls=uniques_trawls+para.bio.haul_no_offset;
        end
        if para.proc.exclude_age1 == 1
            %% exclude age-1 hauls
            [intersect_hauls,IA,IB]= intersect(uniques_trawls,para.proc.age1_haul);
            if ~isempty(intersect_hauls)
                [selected_hauls,IA,IB]= setxor(uniques_trawls,intersect_hauls);
                ind=[];
                for j=1:length(IA)
                    ind0=find(uniques_trawls == selected_hauls(j));
                    ind=[ind; ind0];
                end
                uniques_trawls=uniques_trawls(ind);
            else
                ind=1:length(indx1);
            end
        else
            ind = 1:length(indx1);
        end
        if para.proc.transect_reduction_fraction ~= 0   % reduction factor for simulations of transect spacing reduction 
           intersect_uniques_trawls = intersect(uniques_trawls,data.tmp.TX_selected_trawl_no_indx);
        else
           intersect_uniques_trawls = uniques_trawls;
        end
        
        data.bio.strata(data.bio.haul_strata(i)).trawls=intersect_uniques_trawls;
        data.bio.strata(data.bio.haul_strata(i)).wgt=dat(indx(indx1(ind)),3);
    end
end
%% fill-in trawl information for stratum number < first non-empty stratum 
k1 = 1;
if isempty(data.bio.strata(1).trawls)
    k = 1;
    %% find the first non-empty stratum
    while isempty(data.bio.strata(k).trawls)
        k = k + 1;
        if k > n
            fprintf('++++++++ No non-empty stratum, Processing is stopped !! ++++++++\n');
            break
        end
    end
    k1 = k;   % first non-zero stratum index
    %% fill-in zero-trawl strata
    for i = 1:k-1
        data.bio.strata(i).trawls=data.bio.strata(data.bio.haul_strata(k)).trawls;
        data.bio.strata(i).wgt=data.bio.strata(data.bio.haul_strata(k)).wgt;
        data.bio.haul_strata(i) = i;
    end
% %     %% to deal with additional stratum
% %     %% for example: "Conception" when there were no trawls conducted in the region but there were hake regions
% %     data.bio.haul_strata(data.bio.haul_strata(1:n))=data.bio.haul_strata(1:n);
% %     for i = 1:data.bio.haul_strata(1)-1
% %         data.bio.strata(i).trawls=data.bio.strata(data.bio.haul_strata(1)).trawls;
% %         data.bio.strata(i).wgt=data.bio.strata(data.bio.haul_strata(1)).wgt;
% %         data.bio.haul_strata(i) = i;
% %     end
end

%% interpolate trawl information to stranum (no tralws) by using the parameters of the adjacient strata (stratum)
for i = k1:n
    i1 = i-1;
    k = i;
    while isempty(data.bio.strata(k).trawls)
        k = k + 1;
        if k > n
            break
        end
    end
    i2 = k;
    if i2 <= n   % parameters are from the averaged values from two adjacient strata
        for k = i1+1:i2-1
            %% union of two closest strata (i1 and i2)
            data.bio.strata(k).trawls =[data.bio.strata(i1).trawls ; data.bio.strata(i2).trawls];
            data.bio.strata(k).wgt = [data.bio.strata(i1).wgt; data.bio.strata(i2).wgt];
            data.bio.haul_strata(k) = i;
        end
    else  % using values from non-zero trawls stratum, i.e. stratum i1 = i-1
        data.bio.strata(i).trawls=data.bio.strata(i1).trawls;
        data.bio.strata(i).wgt=data.bio.strata(i1).wgt;
        data.bio.haul_strata(i) = i;
    end
end
%% extrpolate trawl information to stranum (no tralws) > last non-empty stratum (situation in using the Observer data or A-SHOP data)
%% automatically determine the maximum number of stratum   11/2/2021 (SD, and/or Observer uses geographic_stratification xlsx file)
if para.proc.KS_stratification == 1
    n_strata_max =  max(dat0(:,2));     % change max(dat0(:,1)) to max(dat0(:,2)) on 5/29/2021
                                        % change back from max(dat0(:,2)) to max(dat0(:,1)) on 10/18/2021
                                        % change to if loop for determining what which column corresponding strata index
else
    n_strata_max =  max(dat0(:,1));
    if n_strata_max > 20
        n_strata_max =  max(dat0(:,2));
    end
end

for i = data.bio.haul_strata(end)+1:n_strata_max   % change max(dat0(:,1)) to max(dat0(:,2)) on 5/29/2021
                                                   % change back from max(dat0(:,2)) to max(dat0(:,1)) on 10/18/2021
    data.bio.strata(i).trawls=data.bio.strata(data.bio.haul_strata(end)).trawls;
    data.bio.strata(i).wgt=data.bio.strata(data.bio.haul_strata(end)).wgt;
    data.bio.haul_strata(i) = i;
    if length(data.bio.strata(i).trawls) < 1
        disp(data.bio.strata(i).trawls)
    end
end
n=length(data.bio.strata);
% construct hake-haul number array
data.bio.hake_trawl_num=[];
nhaul1=length(data.bio.hake_length_sex);
nhaul2=length(data.bio.hake_length_weight_sex_age);
if nhaul2 > nhaul1
    for i=1:nhaul2
        if i <= nhaul1
             data.bio.hake_trawl_num=union(data.bio.hake_trawl_num,data.bio.hake_length_sex(i).trawl_no);
        end
        data.bio.hake_trawl_num=union(data.bio.hake_trawl_num,data.bio.hake_length_weight_sex_age(i).trawl_no);
    end
else
    for i=1:nhaul1
        if i <= nhaul2
             data.bio.hake_trawl_num=union(data.bio.hake_trawl_num,data.bio.hake_length_weight_sex_age(i).trawl_no);
        end
        data.bio.hake_trawl_num=union(data.bio.hake_trawl_num,data.bio.hake_length_sex(i).trawl_no);
    end    
end
    
data.bio.n_hake_trawl=length(data.bio.hake_trawl_num);
data.bio.len_haul_M=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
data.bio.len_haul_F=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
data.bio.len_haul_ALL=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
data.bio.aged_len_haul_M=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
data.bio.aged_len_haul_F=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
% data.bio.aged_len_haul_N=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
data.bio.aged_len_haul_ALL=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
%% construct a sortable arry len_sex_trawl_array(i) from a non-sortable: data.bio.hake_length_sex(i).trawl_no array 
for i=1:nhaul1
   len_sex_trawl_array(i)=data.bio.hake_length_sex(i).trawl_no;
end
%% construct a sortable arry len_wgt_sex_age_trawl_arry(i) from a non-sortable: hake_length_weight_sex_age(i).trawl_no array 
for i=1:nhaul2
   len_wgt_sex_age_trawl_arry(i)=data.bio.hake_length_weight_sex_age(i).trawl_no;
end

data.bio.unique_aged_haul_no=int64([]);
data.bio.unique_haul_no=int64([]);
%% populate length-haul matrix
for i=1:data.bio.n_hake_trawl
   L=[];L2=[];
   Lmale=[];Lmale2=[];
   Lfemale=[];Lfemale2=[];
   Lunsexed=[];Lunsexed2=[];
%% find trawl index corresponding to the length array at station #1
   ind1=find(len_sex_trawl_array == data.bio.hake_trawl_num(i));
   if ~isempty(ind1)
      L=[L data.bio.hake_length_sex(ind1).length];
      Lmale=[Lmale data.bio.hake_length_sex(ind1).length(data.bio.hake_length_sex(ind1).Gender == 1)];
      Lfemale=[Lfemale data.bio.hake_length_sex(ind1).length(data.bio.hake_length_sex(ind1).Gender == 2)];
%       Lunsexed=[Lunsexed data.bio.hake_length_sex(ind1).length(data.bio.hake_length_sex(ind1).Gender == 3)];
   else
       if para.proc.bio_info_disp == 1
           fprintf('------ no corresponding trawl #%d found in the FSCS log file at Station 1 !!!\n',data.bio.hake_trawl_num(i));
       end
   end
    
%% find trawl index corresponding to the length array at station #2
   ind2=find(len_wgt_sex_age_trawl_arry == data.bio.hake_trawl_num(i));
   data.bio.unique_aged_haul_no(i)=data.bio.hake_trawl_num(i);
   data.bio.unique_haul_no(i)=data.bio.hake_trawl_num(i);
   if ~isempty(ind2)
      ind2_age=find(~isnan(data.bio.hake_length_weight_sex_age(ind2).age) == 1);
      L2=data.bio.hake_length_weight_sex_age(ind2).length(ind2_age);
      if ~isempty(L2)
          L=[L L2];
          Lmale2=data.bio.hake_length_weight_sex_age(ind2).length(data.bio.hake_length_weight_sex_age(ind2).sex(ind2_age) == 1);
          Lmale=[Lmale Lmale2];
          Lfemale2=data.bio.hake_length_weight_sex_age(ind2).length(data.bio.hake_length_weight_sex_age(ind2).sex(ind2_age) == 2);
          Lfemale=[Lfemale Lfemale2];         
      end
   else
       if para.proc.bio_info_disp == 1
           fprintf('------ no corresponding trawl #%d found in the FSCS log file at Station 2 !!!\n',data.bio.hake_trawl_num(i));
       end
   end
   binned_L=hist(L,para.bio.hake_len_bin);
   data.bio.len_haul_ALL(:,i)= binned_L(:);
   binned_Lmale=hist(Lmale,para.bio.hake_len_bin);
   data.bio.len_haul_M(:,i)= binned_Lmale(:);
   binned_Lfemale=hist(Lfemale,para.bio.hake_len_bin);
   data.bio.len_haul_F(:,i)= binned_Lfemale(:);
   binned_Lunsexed=hist(Lunsexed,para.bio.hake_len_bin);
   data.bio.len_haul_N(:,i)= binned_Lunsexed(:);
   
%% length-haul matrix with aged hake only
   if ~isempty(L2)
       binned_L2=hist(L2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_ALL(:,i)= binned_L2(:);
   end
   if ~isempty(Lmale2)
       binned_Lmale2=hist(Lmale2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_M(:,i)= binned_Lmale2(:);
   end
   if ~isempty(Lfemale2)
       binned_Lfemale2=hist(Lfemale2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_F(:,i)= binned_Lfemale2(:);
   end
   if ~isempty(Lunsexed2)
       binned_Lunsexed2=hist(Lunsexed2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_N(:,i)= binned_Lfemale2(:);
   end

end

% remove hauls with no age but only length values in biodata_specimen file
ind_nan=find(sum(data.bio.aged_len_haul_ALL)==0);
data.bio.unique_aged_haul_no(ind_nan)=[];
data.bio.aged_len_haul_M(:,ind_nan)=[];
data.bio.aged_len_haul_F(:,ind_nan)=[];
data.bio.aged_len_haul_ALL(:,ind_nan)=[];

% remove hauls with no length values in biodata_specimen file
ind_nan=find(sum(data.bio.len_haul_ALL)==0);
% data.bio.unique_haul_no(ind_nan)=int64([]); % remove int64() on 2019-12-18 
data.bio.unique_haul_no(ind_nan)= [];
data.bio.len_haul_M(:,ind_nan)=[];
data.bio.len_haul_F(:,ind_nan)=[];
data.bio.len_haul_ALL(:,ind_nan)=[];


for i=1:length(para.bio.hake_len_bin)
    indx=find(data.bio.len_haul_M(i,:) ~= 0);
    data.bio.compact_len_haul_M(i,1:length(indx))=data.bio.hake_trawl_num(indx);
    indx=find(data.bio.len_haul_F(i,:) ~= 0);
    data.bio.compact_len_haul_F(i,1:length(indx))=data.bio.hake_trawl_num(indx);
    indx=find(data.bio.len_haul_ALL(i,:) ~= 0);
    data.bio.compact_len_haul_ALL(i,1:length(indx))=data.bio.hake_trawl_num(indx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create and fill-in strata structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
        fprintf('----- create and fill-in strata structure ...\n')
    end
end

for ii = 1:n   % number of strata
    if isfield(para, 'platform_name')
        if strcmp(para.platform_name, 'SD')  % SD data
            fprintf('  Stratum %d -> elapse time = %4.0f (sec) \n', ii, etime(clock, t0))
        end
    end
    if para.bio_data_type == 3 & data.bio.haul_strata(ii) == 100
        i = data.bio.haul_strata(ii);   %  stratum number
    else
        i = ii;
        %     if i == 5
        %         disp(i)
        %     end
        L1=[];
        Wgt1=0;
        Gender1=[];
        L2=[];
        Wgt2=[];
        Gender2=[];
        Age2=[];
        nt=length(data.bio.strata(i).trawls);
        data.bio.strata(i).Len_key_Mx=zeros(1,length(para.bio.hake_len_bin));
        data.bio.strata(i).Len_key_Fx=zeros(1,length(para.bio.hake_len_bin));
        data.bio.strata(i).Len_key_Nx=zeros(1,length(para.bio.hake_len_bin));
        data.bio.strata(i).Len_key_Mn=zeros(1,length(para.bio.hake_len_bin));
        data.bio.strata(i).Len_key_Fn=zeros(1,length(para.bio.hake_len_bin));
        data.bio.strata(i).Len_key_Nn=zeros(1,length(para.bio.hake_len_bin));
        for j=1:nt   % loop through trawls
            %%%%% station #1 -> length - gender
            L1j=[];
            for k=1:length(data.bio.catch)
                if data.bio.catch(k).trawl_no == data.bio.strata(i).trawls(j)
                    catch_trawl_no=k;
                    break
                end
            end
            for k=1:length(data.bio.hake_length_sex)  % find the coresponding trawl in the stratum
                if data.bio.strata(i).trawls(j) == data.bio.hake_length_sex(k).trawl_no
                    data.bio.strata(i).num_of_fish(j)=data.bio.hake_length_sex(k).n;
                    species_ind=0;
                    for l=1:length(data.bio.catch(catch_trawl_no).species)
                        if data.bio.catch(catch_trawl_no).species(l).ID == para.bio.species_code_ID
                            species_ind=l;
                        end
                    end
                    if species_ind == 0
                        if para.proc.bio_info_disp == 1
                            fprintf('---- cannot find hake in this chosen trawl!!! ----')
                            fprintf('trawl = %d\n', catch_trawl_no);
                        end
                    else
                        Wgt1=Wgt1+data.bio.catch(catch_trawl_no).species(species_ind).exp_wgt;
                    end
                    
                    L1j=round(data.bio.hake_length_sex(k).length);
                    % accumulative quantities over trawls within the stratum i from staiotn #1
                    L1=[L1 L1j];
                    Gender1=[Gender1 data.bio.hake_length_sex(k).Gender];
                    data.bio.hake_length_sex(k).stratum=i;
                end
            end
            
            L2j=[];
            Wgt2j=[];
            Gender2j=[];
            Age2j=[];
            
            %%%%% station #2 -> length - weight - gender - age
            for k=1:length(data.bio.hake_length_weight_sex_age)  % find the coresponding trawl in the stratum
                if data.bio.strata(i).trawls(j) == data.bio.hake_length_weight_sex_age(k).trawl_no
                    L2j=round(data.bio.hake_length_weight_sex_age(k).length);
                    Gender2j=data.bio.hake_length_weight_sex_age(k).sex;
                    Wgt2j=data.bio.hake_length_weight_sex_age(k).wgt;
                    Age2j=data.bio.hake_length_weight_sex_age(k).age;
                    indx_nan=find(isnan(Wgt2j) ==1);
                    Wgt2j(indx_nan)=[];
                    L2j(indx_nan)=[];
                    Age2j(indx_nan)=[];
                    Gender2j(indx_nan)=[];
                    data.bio.hake_length_weight_sex_age(k).stratum=i;
                    L2=[L2 L2j];
                    Gender2=[Gender2 Gender2j];
                    Wgt2=[Wgt2 Wgt2j];
                    Age2=[Age2 Age2j];
                    indx_age_nan=find(isnan(Age2j) ==1);
                    if ~isempty(indx_age_nan) &  para.proc.bio_info_disp == 1
                        fprintf('trawl # = %d\t  # of Nan Age = %d !!\n',data.bio.hake_length_weight_sex_age(k).trawl_no,length(indx_age_nan))
                    end
                end
            end
            if str2num(para.survey_year) >= 2003 & para.acoust.TS_station_num == 1
                TS0j=20*log10(L1j)-68;
                TSj=10*log10(nanmean(10.^(TS0j/10)));
                data.bio.strata(i).Lave_j(j)=nanmean(L1j);
                data.bio.strata(i).TS_lin_j(j)=TSj;                   % 10*log10(<sv>
                data.bio.strata(i).sig_bs_j(j)=10.^(TSj/10);
            else
                TS0j=20*log10([L1j L2j])-68;
                TSj=10*log10(nanmean(10.^(TS0j/10)));
                data.bio.strata(i).Lave_j(j)=nanmean([L1j L2j]);
                data.bio.strata(i).TS_lin_j(j)=TSj;                   % 10*log10(<sv>
                data.bio.strata(i).sig_bs_j(j)=10.^(TSj/10);
            end
            if isempty(Age2) &  para.proc.bio_info_disp == 1
                fprintf('No Age2+ samples at station #1: Stratum = %d\t Trawl = %d\n',i,data.bio.strata(i).trawls(j));
            end
        end                                 % end trawl loop within the stratum
        
        nM=length(find(Gender1 == 1));
        nF=length(find(Gender1 == 2));
        nTot=length(L1);
        
        if nTot+ length(L2) > 0  % ignore un-defined stratum
            %    fprintf('stratum %d\t   M: %d\t   F: %d\t  Ratio: %4.2f\t Total: %d\n',i,nM,nF,nF/nM,nTot);
            %% obtain biological inofrmation within each stratum
            %%  biological sampling station #1
            data.bio.strata(i).length1=L1;
            data.bio.strata(i).gender1=Gender1;
            data.bio.strata(i).meanLen1=mean(L1);
            data.bio.strata(i).stdLen1=std(L1);
            data.bio.strata(i).wgt1=Wgt1;
            
            %%  biological sampling station #2
            data.bio.strata(i).length2=L2;
            data.bio.strata(i).gender2=Gender2;
            data.bio.strata(i).weight2=Wgt2;               % weight array of individual fish weight
            data.bio.strata(i).wgt2=nansum(Wgt2);
            if all(isnan(Age2))
                fprintf('No Age2+ samples at station #2: Stratum = %d\t Trawl = %d\n',i,data.bio.strata(i).trawls(j));
                fprintf('Solution:  find age from length ... \n');
                if isempty(L2)  % station 2 samples for stratum i
                    [L1, Gender1, L2, Age2, Wgt2, Gender2, data] = get_stratumB_from_stratumA(data, i, i-1);
                else
                    Age2 = find_age_from_length(data, para, L2);
                end
            end
            data.bio.strata(i).age2=Age2;
            
            if  para.acoust.TS_station_num == 1 & nTot > 0  % 2015-12-15 edited to take into account the case when nTot = 0, and para.acoust.TS_station_num = 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   using station 1 data only for TS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                L=round(L1);
                Gender=Gender1;
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   using both station 1 and 2 data for TS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                L=round([L1 L2]);
                Gender=[Gender1 Gender2];
            end
            % normalized within each trawl
            data.bio.strata(i).sig_bs=nanmean(data.bio.strata(i).sig_bs_j);
            data.bio.strata(i).sig_b=4*pi*data.bio.strata(i).sig_bs;
            data.bio.strata(i).gender=Gender;
            data.bio.strata(i).L=L;
            data.bio.strata(i).Lave=nanmean(data.bio.strata(i).Lave_j);
            data.bio.strata(i).Male_ind=find(data.bio.strata(i).gender == 1);
            data.bio.strata(i).Female_ind=find(data.bio.strata(i).gender == 2);
            data.bio.strata(i).Gender(data.bio.strata(i).Male_ind)=1;
            data.bio.strata(i).Gender(data.bio.strata(i).Female_ind)=2;
            data.bio.strata(i).LaveM=nanmean(data.bio.strata(i).L(data.bio.strata(i).Male_ind));
            data.bio.strata(i).LaveF=nanmean(data.bio.strata(i).L(data.bio.strata(i).Female_ind));
            data.bio.strata(i).LaveALL=nanmean(L);
            data.bio.strata(i).nM=length(data.bio.strata(i).Male_ind);
            data.bio.strata(i).nF=length(data.bio.strata(i).Female_ind);
            data.bio.strata(i).nALL=length(L);
            %% quantities associated with sample station #1
            data.bio.strata(i).n=length(L1);
            data.bio.strata(i).L1=L1;
            male_indx=find(data.bio.strata(i).gender1 == 1);
            female_indx=find(data.bio.strata(i).gender1 == 2);
            data.bio.strata(i).nM1=length(male_indx);
            data.bio.strata(i).nF1=length(female_indx);
            if data.bio.strata(i).nM1 == 0
                male_wgt=0;
            else
                L1_M=L1(male_indx);
                male_wgt=nansum(interp1(para.bio.hake_len_bin,data.bio.len_wgt_M,L1_M));
            end
            if data.bio.strata(i).nF1 == 0
                female_wgt=0;
            else
                L1_F=L1(female_indx);
                female_wgt=nansum(interp1(para.bio.hake_len_bin,data.bio.len_wgt_F,L1_F));
            end
            if  data.bio.strata(i).nF1 == 0 &  data.bio.strata(i).nF1 == 0
                data.bio.strata(i).nM_wgt1=0;
                data.bio.strata(i).nF_wgt1=0;
            else
                data.bio.strata(i).nM_wgt1=Wgt1*male_wgt/(male_wgt+female_wgt);
                data.bio.strata(i).nF_wgt1=Wgt1*female_wgt/(male_wgt+female_wgt);
            end
            
            data.bio.strata(i).n1=length(L1);
            if data.bio.strata(i).nM+data.bio.strata(i).nF ~= data.bio.strata(i).n
                if data.bio.strata(i).n ~= 0  &  para.proc.bio_info_disp == 1
                    fprintf('Station 1: stratum: %d\t Unsexed samples = %f %% \n',i, 1-(data.bio.strata(i).nM+data.bio.strata(i).nF)/data.bio.strata(i).n);
                end
            end
            
            % create sex-length key (station 1 samples)
            % 2015-12-15 edited to take into account the case when nTot (from station 1 only) = 0, and para.acoust.TS_station_num = 1
            if ~isempty(Gender1) | (isempty(Gender1) & para.acoust.TS_station_num == 1)
                %% 2015-12-15 revisions
                if ~isempty(Gender1)
                    LenM_ind=find(data.bio.strata(i).gender1 == 1);
                    LenF_ind=find(data.bio.strata(i).gender1 == 2);
                    data.bio.strata(i).Len_key_Mx=hist(data.bio.strata(i).length1(LenM_ind),para.bio.hake_len_bin);
                    data.bio.strata(i).Len_key_Fx=hist(data.bio.strata(i).length1(LenF_ind),para.bio.hake_len_bin);
                    data.bio.strata(i).Len_key_Nx=hist(data.bio.strata(i).length1,para.bio.hake_len_bin);               % un_aged samples
                    num_len=data.bio.strata(i).n1;
                else
                    LenM_ind=find(data.bio.strata(i).gender2 == 1);
                    LenF_ind=find(data.bio.strata(i).gender2 == 2);
                    data.bio.strata(i).Len_key_Mx=hist(data.bio.strata(i).length2(LenM_ind),para.bio.hake_len_bin);
                    data.bio.strata(i).Len_key_Fx=hist(data.bio.strata(i).length2(LenF_ind),para.bio.hake_len_bin);
                    data.bio.strata(i).Len_key_Nx=hist(data.bio.strata(i).length2,para.bio.hake_len_bin);               % un_aged samples
                    num_len=data.bio.strata(i).nALL;                                % number of samples in stations 1 & 2
                end
                if length(LenM_ind) > 0
                    data.bio.strata(i).Len_key_Mn=data.bio.strata(i).Len_key_Mx/length(LenM_ind);
                else
                    data.bio.strata(i).Len_key_Mn=zeros(size(data.bio.strata(i).Len_key_Mx));
                end
                if length(LenF_ind) > 0
                    data.bio.strata(i).Len_key_Fn=data.bio.strata(i).Len_key_Fx/length(LenF_ind);
                else
                    data.bio.strata(i).Len_key_Fn==zeros(size(data.bio.strata(i).Len_key_Fx));
                end
                %        data.bio.strata(i).Len_key_Nn=data.bio.strata(i).Len_key_Nx/(length(LenM_ind)+length(LenF_ind));
                %             data.bio.strata(i).Len_key_Nn=data.bio.strata(i).Len_key_Nx/data.bio.strata(i).n1;     % modified on 10/25/2010
                data.bio.strata(i).Len_key_Nn=data.bio.strata(i).Len_key_Nx/num_len;                                 % modified on 12/15/2015
            else
                if  para.proc.bio_info_disp == 1
                    fprintf('---- no sexed data in stratum  %d\n',i);
                end
            end
            
            % create sex-length-age matrix  (station 2 samples)
            indx2=[];
            for j=1:length(para.bio.hake_age_bin)
                for k=1:length(para.bio.hake_len_bin)
                    % Male
                    if k == length(para.bio.hake_len_bin)
                        L2jkM_ind=find(Gender2 == 1 & Age2 == para.bio.hake_age_bin(j) & L2 >= para.bio.hake_len_bin(k));
                    else
                        L2jkM_ind=find(Gender2 == 1 & Age2 == para.bio.hake_age_bin(j) & (L2 >= para.bio.hake_len_bin(k) & L2 < para.bio.hake_len_bin(k+1)));
                    end
                    % number key
                    data.bio.strata(i).Len_Age_key_M(k,j)=length(L2jkM_ind);
                    % weight key
                    data.bio.strata(i).Len_Age_key_wgt_M(k,j)=sum(Wgt2(L2jkM_ind));
                    indx2=union(indx2,L2jkM_ind);
                    % Female
                    if k == length(para.bio.hake_len_bin)
                        L2jkF_ind=find(Gender2 == 2 & Age2 == para.bio.hake_age_bin(j) & L2 >= para.bio.hake_len_bin(k));
                    else
                        L2jkF_ind=find(Gender2 == 2 & Age2 == para.bio.hake_age_bin(j) & (L2 >= para.bio.hake_len_bin(k) & L2 < para.bio.hake_len_bin(k+1)));
                    end
                    % number key
                    data.bio.strata(i).Len_Age_key_F(k,j)=length(L2jkF_ind);
                    % weight key
                    data.bio.strata(i).Len_Age_key_wgt_F(k,j)=sum(Wgt2(L2jkF_ind));
                    indx2=union(indx2,L2jkF_ind);
                    % ALL
                    % number key
                    data.bio.strata(i).Len_Age_key_ALL(k,j)=data.bio.strata(i).Len_Age_key_M(k,j)+data.bio.strata(i).Len_Age_key_F(k,j);
                    % weight key
                    data.bio.strata(i).Len_Age_key_wgt_ALL(k,j)=data.bio.strata(i).Len_Age_key_wgt_M(k,j)+data.bio.strata(i).Len_Age_key_wgt_F(k,j);
                end
            end
            total_matrix_wgt2=nansum(nansum(data.bio.strata(i).Len_Age_key_wgt_ALL));
            total_wgt2=sum(Wgt2);
            if abs(total_matrix_wgt2-total_wgt2) > 1e-2 &  para.proc.bio_info_disp == 1
                fprintf('stratum = %d\t  Matrix wgt2 = %f\t  wgt2 = %f\n',i,total_matrix_wgt2,total_wgt2);
            end
            %% normalized the length-age matrices
            %% length-age number keys
            data.bio.strata(i).matrix_NtotM=nansum(nansum(data.bio.strata(i).Len_Age_key_M));
            data.bio.strata(i).matrix_NtotF=nansum(nansum(data.bio.strata(i).Len_Age_key_F));
            data.bio.strata(i).matrix_NtotALL=data.bio.strata(i).matrix_NtotM+data.bio.strata(i).matrix_NtotF;
            data.bio.strata(i).Len_Age_key_Mn=data.bio.strata(i).Len_Age_key_M/data.bio.strata(i).matrix_NtotM;
            data.bio.strata(i).Len_Age_key_Fn=data.bio.strata(i).Len_Age_key_F/data.bio.strata(i).matrix_NtotF;
            data.bio.strata(i).Len_Age_key_ALLn=data.bio.strata(i).Len_Age_key_ALL/data.bio.strata(i).matrix_NtotALL;
            
            %         if i == 5
            %             disp(i)
            %         end
            %% length-age weight keys
            data.bio.strata(i).matrix_NtotM_wgt=nansum(nansum(data.bio.strata(i).Len_Age_key_wgt_M));
            data.bio.strata(i).matrix_NtotF_wgt=nansum(nansum(data.bio.strata(i).Len_Age_key_wgt_F));
            data.bio.strata(i).matrix_NtotALL_wgt=data.bio.strata(i).matrix_NtotM_wgt+data.bio.strata(i).matrix_NtotF_wgt;
            data.bio.strata(i).Len_Age_key_wgt_Mn=data.bio.strata(i).Len_Age_key_wgt_M/data.bio.strata(i).matrix_NtotM_wgt;
            data.bio.strata(i).Len_Age_key_wgt_Fn=data.bio.strata(i).Len_Age_key_wgt_F/data.bio.strata(i).matrix_NtotF_wgt;
            data.bio.strata(i).Len_Age_key_wgt_ALLn=data.bio.strata(i).Len_Age_key_wgt_ALL/data.bio.strata(i).matrix_NtotALL_wgt;
            
            data.bio.strata(i).total_N=data.bio.strata(i).n+data.bio.strata(i).matrix_NtotALL;   % station #1 + station #2
            data.bio.strata(i).total_wgt=data.bio.strata(i).wgt1+data.bio.strata(i).wgt2;   % station #1 + station #2
            % len and len-Age number key proportions
            data.bio.strata(i).Len_Age_M_proportion=data.bio.strata(i).matrix_NtotM/data.bio.strata(i).total_N;
            data.bio.strata(i).Len_Age_F_proportion=data.bio.strata(i).matrix_NtotF/data.bio.strata(i).total_N;
            data.bio.strata(i).Len_M_proportion=data.bio.strata(i).nM1/data.bio.strata(i).total_N;
            data.bio.strata(i).Len_F_proportion=data.bio.strata(i).nF1/data.bio.strata(i).total_N;
            
            % len and len-Age weight key proprotions
            data.bio.strata(i).Len_Age_M_wgt_proportion=data.bio.strata(i).matrix_NtotM_wgt/data.bio.strata(i).total_wgt;
            data.bio.strata(i).Len_Age_F_wgt_proportion=data.bio.strata(i).matrix_NtotF_wgt/data.bio.strata(i).total_wgt;
            data.bio.strata(i).Len_M_wgt_proportion=data.bio.strata(i).nM_wgt1/data.bio.strata(i).total_wgt;
            data.bio.strata(i).Len_F_wgt_proportion=data.bio.strata(i).nF_wgt1/data.bio.strata(i).total_wgt;
            
            check1=data.bio.strata(i).Len_Age_M_proportion+data.bio.strata(i).Len_Age_F_proportion+data.bio.strata(i).Len_M_proportion+data.bio.strata(i).Len_F_proportion;
            check2=data.bio.strata(i).Len_Age_M_wgt_proportion+data.bio.strata(i).Len_Age_F_wgt_proportion+data.bio.strata(i).Len_M_wgt_proportion+data.bio.strata(i).Len_F_wgt_proportion;
            if abs(check1+check2 - 2) > 0.001 &  para.proc.bio_info_disp == 1
                fprintf('un-sexed samples are found in stratum: %d\t number check = %f\t weight check = %f\n',i,check1,check2);
            end
            %% generate len-wgt-key for each stratum (only use the data from sample sation #2)
            indx2=[];
            for k=1:length(para.bio.hake_len_bin)
                % Male
                if k == length(para.bio.hake_len_bin)
                    L2kM_ind=find(Gender2 == 1 & data.bio.strata(i).length2 >= para.bio.hake_len_bin(k));
                else
                    L2kM_ind=find(Gender2 == 1 & (data.bio.strata(i).length2 >= para.bio.hake_len_bin(k) & data.bio.strata(i).length2 < para.bio.hake_len_bin(k+1)));
                end
                if length(L2kM_ind) >= 5
                    data.bio.strata(i).Len_wgt_key_M(k)=nanmean(data.bio.strata(i).weight2(L2kM_ind));
                else
                    data.bio.strata(i).Len_wgt_key_M(k)=data.bio.len_wgt_all.reg_w0M*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_pM;
                end
                % Female
                if k == length(para.bio.hake_len_bin)
                    L2kF_ind=find(Gender2 == 2 & data.bio.strata(i).length2 >= para.bio.hake_len_bin(k));
                else
                    L2kF_ind=find(Gender2 == 2 & (data.bio.strata(i).length2 >= para.bio.hake_len_bin(k) & data.bio.strata(i).length2 < para.bio.hake_len_bin(k+1)));
                end
                if length(L2kF_ind) >= 5
                    data.bio.strata(i).Len_wgt_key_F(k)=nanmean(data.bio.strata(i).weight2(L2kF_ind));
                else
                    data.bio.strata(i).Len_wgt_key_F(k)=data.bio.len_wgt_all.reg_w0F*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_pF;
                end
                data.bio.strata(i).Len_wgt_key_F_ind(k)=length(L2kF_ind);
                % ALL
                ALL_nk=length(L2kM_ind)+length(L2kF_ind);
                if length(ALL_nk) >= 5
                    data.bio.strata(i).Len_wgt_key_ALL(k)=(data.bio.strata(i).Len_wgt_key_M(k)*length(L2kM_ind)+data.bio.strata(i).Len_wgt_key_F(k)*length(L2kF_ind))/ALL_nk;
                else
                    data.bio.strata(i).Len_wgt_key_ALL(k)=data.bio.len_wgt_all.reg_w0*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_p;
                end
                data.bio.strata(i).Len_wgt_key_ALL_ind(k)=ALL_nk;
            end
            if isnan(sum(sum(data.bio.strata(i).Len_Age_key_ALLn)))
                disp([ i ii])
            end
        else
            fprintf('Stratum %d has no trawl samples !!\n',i);
        end  % end of ignore un-defined stratum
    end   %% end of if loop for stratum index == 0
%     if length(data.bio.strata(i).Len_key_Mn) == 1
%         disp(i)
%     end
end                             % end of strata loop

if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
        fprintf('  Stratum %d -> elapse time = %4.0f (sec) \n', ii, etime(clock, t0))
    end
end

if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
        fprintf('=== Start copyiny stratum information to adjacient stratum -> elapse time = %4.0f (sec) \n', etime(clock, t0))
    end
end
%% modifiy data.bio.strata struct if it is empty due to removed trawls 
if para.proc.transect_reduction_fraction ~= 0
    stratum_id_i = 1:max(stratum_id);
    for i = stratum_id_i
        if isempty(data.bio.strata(i).sig_b) | (para.proc.exclude_age1 == 1 & max(data.bio.strata(i).age2) < 2)
            fprintf('no sig_b for stratum = %d\n', i)
            if i == 1 
                j = 2;
                while (isempty(data.bio.strata(j).sig_b) | max(data.bio.strata(j).age2) < 2)
                    j = j +  1;
                end
                %% copy contents of stratum j to stratum i
                copy_stratum_contents(i, j);   
            elseif i == stratum_id_i(end)
                j = stratum_id_i(end) - 1;
                while (isempty(data.bio.strata(j).sig_b) | max(data.bio.strata(j).age2) < 2)
                    j = j -  1;
                end
                %% copy contents of stratum j to stratum i
                copy_stratum_contents(i, j);   
            else
                j1 = i - 1;
                while (isempty(data.bio.strata(j1).sig_b) | max(data.bio.strata(j1).age2) < 2)
                    j1 = j1 -  1;
                end        
                found_nonempty_stratum = 0;
                j2 = i + 1;
                while (isempty(data.bio.strata(j2).sig_b) | max(data.bio.strata(j2).age2) < 2)
                    j2 = j2 +  1;
                    if j2 > stratum_id_i(end)
                        j2 = j2 - 1;
                        break
                    end
                end
                if (isempty(data.bio.strata(j2).sig_b) | max(data.bio.strata(j2).age2) < 2)
                %% copy contents of stratum j1 to stratum i
                    copy_stratum_contents(i, j1);
                else
                    if i-j1 < j2 - i   % stratum i close to stratum j1
                %% copy contents of stratum j1 to stratum i
                        copy_stratum_contents(i, j1);
                    elseif i-j1 > j2 - i  % stratum i close to stratum j2
                %% copy contents of stratum j2 to stratum i
                        copy_stratum_contents(i, j2);
                    else
                        if length(data.bio.strata(j1).trawls) >= length(data.bio.strata(j2).trawls)
                        % stratum i has same "stratum distance" to stratum j1 and to stratum j2, but there are more trawls in stratum j1 than in stratum j2
                        %% copy contents of stratum j1 to stratum i
                            copy_stratum_contents(i, j1);
                        else
                        % stratum i has same "stratum distance" to stratum j1 and to stratum j2, but there are more trawls in stratum j2 than in stratum j1
                        %% copy contents of stratum j2 to stratum i
                            copy_stratum_contents(i, j2);
                        end
                    end
                %% copy contents of averages of strata j1 and j2 to stratum i
% %                     copy_stratum_contents(i, j1, j2);
                end
            end
        end
%         if length(data.bio.strata(i).Len_key_Mn) == 1
%             disp(i)
%         end
    end  % end of strata llop
end   % end of "if transect_reduction"
%% double check if there are zero-len-sge-key
for  i = 1:max(stratum_id)
    %% age1+
    if sum(sum(data.bio.strata(i).Len_Age_key_ALLn)) == 0
        fprintf('Age1+ zero-key: Stratum = %d\n', i)
    else
        % age2+
        if sum(sum(data.bio.strata(i).Len_Age_key_ALLn(:,2:end))) == 0
            fprintf('Age2+ zero-key: Stratum = %d\n', i)
        end
    end
end
if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
        fprintf('*****************************  Done with Strata  ***************************** \n', ii, etime(clock, t0))
    end
end
return