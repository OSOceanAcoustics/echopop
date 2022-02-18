function      get_historical_strata_data_exclude1
% get parameters associated with strata for historical data (1995 -2001):
% trawl information
% Length - key
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013


global data para

fid=fopen(para.acoust.filename.strata);
if fid < 0
    return
end

i=0;
while ~feof(fid)
   i=i+1;
   line=str2num(fgetl(fid));
end
fclose(fid);

%%% read haul-strata data file
para.bio.strata=[];
file_length=i;
dat0 = xlsread(para.acoust.filename.strata,'',['A2:E' num2str(file_length)]);
[junk, strata_name] = xlsread(para.acoust.filename.strata,'',['C2:C' num2str(file_length)]);
indx=find(str2num(para.survey_year ) == dat0(:,1));
dat=dat0(indx,[2 4 5]);
data.bio_acoust.haul_wgt_tbl=dat(:,2:3);
data.bio.haul_strata=unique(dat(:,1));
n=length(data.bio.haul_strata);                        % number of strata
for i=1:n
    indx=find(dat(:,1) == data.bio.haul_strata(i));
    [uniques_trawls indx1]=unique(dat(indx,2));
    if ~isempty(para.bio.CAN_strata_num0) & data.bio.haul_strata(i) >= para.bio.CAN_strata_num0
   %% combining CAN hauls and add haul number offset is haul numbers are
   %% overlapping with those from US 
        uniques_trawls=uniques_trawls+para.bio.haul_no_offset;
    end
    para.bio.strata(data.bio.haul_strata(i)).trawls=uniques_trawls;
    para.bio.strata(data.bio.haul_strata(i)).wgt=dat(indx(indx1),3);
    for j=1:length(uniques_trawls)     % # of hauls in stratum i
        strata_trawls=para.bio.strata(data.bio.haul_strata(i)).trawls;
    end
end

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
data.bio.aged_len_haul_ALL=zeros(length(para.bio.hake_len_bin),data.bio.n_hake_trawl);
%% construct a sortable arry len_sex_trawl_array(i) from a non-sortable: data.bio.hake_length_sex(i).trawl_no array 
for i=1:nhaul1
   len_sex_trawl_array(i)=data.bio.hake_length_sex(i).trawl_no;
end
%% construct a sortable arry len_wgt_sex_age_trawl_arry(i) from a non-sortable: hake_length_weight_sex_age(i).trawl_no array 
for i=1:nhaul2
   len_wgt_sex_age_trawl_arry(i)=data.bio.hake_length_weight_sex_age(i).trawl_no;
end

%% populate length-haul matrix
for i=1:data.bio.n_hake_trawl
   L=[];L2=[];
   Lmale=[];Lmale2=[];
   Lfemale=[];Lfemale2=[];
%% find trawl index corresponding to the length array at station #1
   ind1=find(len_sex_trawl_array == data.bio.hake_trawl_num(i));
   if ~isempty(ind1)
      L=[L data.bio.hake_length_sex(ind1).length];
      Lmale=[Lmale data.bio.hake_length_sex(ind1).length(data.bio.hake_length_sex(ind1).Gender == 1)];
      Lfemale=[Lfemale data.bio.hake_length_sex(ind1).length(data.bio.hake_length_sex(ind1).Gender == 2)];
   else
      fprintf('222 ------ no corresponding trawl #%d found in the FSCS log file at Station 1 !!!\n',data.bio.hake_trawl_num(i));
   end
    
%% find trawl index corresponding to the length array at station #2
   ind2=find(len_wgt_sex_age_trawl_arry == data.bio.hake_trawl_num(i));
   if ~isempty(ind2)
      L2=data.bio.hake_length_weight_sex_age(ind2).length;
      L=[L L2];
      Lmale2=data.bio.hake_length_weight_sex_age(ind2).length(data.bio.hake_length_weight_sex_age(ind2).sex == 1);
      Lmale=[Lmale Lmale2];
      Lfemale2=data.bio.hake_length_weight_sex_age(ind2).length(data.bio.hake_length_weight_sex_age(ind2).sex == 2);
      Lfemale=[Lfemale Lfemale2];
   else
      fprintf('------ no corresponding trawl #%d found in the FSCS log file at Station 2 !!!\n',data.bio.hake_trawl_num(i));
   end
   binned_L=hist(L,para.bio.hake_len_bin);
   data.bio.len_haul_ALL(:,i)= binned_L(:);
   binned_Lmale=hist(Lmale,para.bio.hake_len_bin);
   data.bio.len_haul_M(:,i)= binned_Lmale(:);
   binned_Lfemale=hist(Lfemale,para.bio.hake_len_bin);
   data.bio.len_haul_F(:,i)= binned_Lfemale(:);
   
%% length-haul matrix with aged hake only
   if ~isempty(L2)
       binned_L2=hist(L2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_ALL(:,i)= binned_L2(:);
   else
   end
   if ~isempty(Lmale2)
       binned_Lmale2=hist(Lmale2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_M(:,i)= binned_Lmale2(:);
   end
   if ~isempty(Lfemale2)
       binned_Lfemale2=hist(Lfemale2,para.bio.hake_len_bin);
       data.bio.aged_len_haul_F(:,i)= binned_Lfemale2(:);
   end

end

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
for i=1:n   % number of strata
    L1=[];
    Wgt1=0;
    Gender1=[];
    L2=[];
    Wgt2=[];
    Gender2=[];
    Age2=[];
    nt=length(para.bio.strata(i).trawls);
    para.bio.strata(i).Len_key_Mx=zeros(1,length(para.bio.hake_len_bin));
    para.bio.strata(i).Len_key_Fx=zeros(1,length(para.bio.hake_len_bin));
    para.bio.strata(i).Len_key_Nx=zeros(1,length(para.bio.hake_len_bin));
    para.bio.strata(i).Len_key_Mn=zeros(1,length(para.bio.hake_len_bin));
    para.bio.strata(i).Len_key_Fn=zeros(1,length(para.bio.hake_len_bin));
    para.bio.strata(i).Len_key_Nn=zeros(1,length(para.bio.hake_len_bin));
    for j=1:nt   % loop through trawls
      %%%%% station #1 -> length - gender
        L1j=[];
        for k=1:length(data.bio.catch)
            if data.bio.catch(k).trawl_no == para.bio.strata(i).trawls(j)
                catch_trawl_no=k;
                break
            end
        end
        for k=1:length(data.bio.hake_length_sex)  % find the coresponding trawl in the stratum
           if para.bio.strata(i).trawls(j) == data.bio.hake_length_sex(k).trawl_no
               para.bio.strata(i).num_of_fish(j)=data.bio.hake_length_sex(k).n;
               species_ind=0;
               for l=1:length(data.bio.catch(catch_trawl_no).species)
                  if data.bio.catch(catch_trawl_no).species(l).ID == para.bio.species_code_ID
                      species_ind=l;
                  end
               end
               if species_ind ==0
                   fprintf('---- cannot find hake in this chosen trael!!! ----')
               end
               Wgt1=Wgt1+data.bio.catch(catch_trawl_no).species(species_ind).exp_wgt;
               
               L1j=data.bio.hake_length_sex(k).length;
             % accumulative quantities over trawls within the stratum i from staiotn #1
               L1=[L1 L1j];
               Gender1=[Gender1 data.bio.hake_length_sex(k).Gender];
               data.bio.hake_length_sex(k).stratum=i;
%                if i == 12
%                    disp(i)
%                end
           end
        end      

        L2j=[];
        Wgt2j=[];
        Gender2j=[];
        Age2j=[];
        
      %%%%% station #2 -> length - weight - gender - age
        for k=1:length(data.bio.hake_length_weight_sex_age)  % find the coresponding trawl in the stratum
           if para.bio.strata(i).trawls(j) == data.bio.hake_length_weight_sex_age(k).trawl_no
               L2j=data.bio.hake_length_weight_sex_age(k).length;
               Gender2j=data.bio.hake_length_weight_sex_age(k).sex;
               Wgt2j=data.bio.hake_length_weight_sex_age(k).wgt;
               Age2j=data.bio.hake_length_weight_sex_age(k).age;
               if para.proc.exclude_age1 == 1 
                   indx_nan=find(isnan(Wgt2j) ==1 | Age2j == 1);
               else
                   indx_nan=find(isnan(Wgt2j) ==1);
               end
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
               if ~isempty(indx_age_nan)
                   fprintf('trawl # = %d\t  # of Nan Age = %d !!\n',data.bio.hake_length_weight_sex_age(k).trawl_no,length(indx_age_nan))
               end
           end
        end
        if str2num(para.survey_year) >= 2003 & para.acoust.TS_station_num == 1
            TS0j=20*log10(L1j)-68;
            TSj=10*log10(nanmean(10.^(TS0j/10)));
            para.bio.strata(i).Lave_j(j)=nanmean(L1j);
            para.bio.strata(i).TS_lin_j(j)=TSj;                   % 10*log10(<sv>
            para.bio.strata(i).sig_bs_j(j)=10.^(TSj/10);
%             if i == 12
%                 fprintf('Stratum = %d\t Trawl # = %d\t n = %d\t  %6.4e \n',i,para.bio.strata(i).trawls(j),length(L1j),para.bio.strata(i).sig_bs_j(j));
%             end
        else
            TS0j=20*log10([L1j L2j])-68;
            TSj=10*log10(nanmean(10.^(TS0j/10)));
            para.bio.strata(i).Lave_j(j)=nanmean([L1j L2j]);
            para.bio.strata(i).TS_lin_j(j)=TSj;                   % 10*log10(<sv>
            para.bio.strata(i).sig_bs_j(j)=10.^(TSj/10);
            if j == 8
  %              disp(L1j)
            end
        end
        if isempty(Age2)
 %           fprintf('No samples at station #1: Stratum = %d\t Trawl = %d\n',i,para.bio.strata(i).trawls(j));
        end
    end                                 % end trawl loop within the stratum
 
    nM=length(find(Gender1 == 1));
    nF=length(find(Gender1 == 2));
    nTot=length(L1);
%    fprintf('stratum %d\t   M: %d\t   F: %d\t  Ratio: %4.2f\t Total: %d\n',i,nM,nF,nF/nM,nTot);
%% obtain biological inofrmation within each stratum
%%  biological sampling station #1
    para.bio.strata(i).length1=L1;
    para.bio.strata(i).gender1=Gender1;
    para.bio.strata(i).meanLen1=mean(L1);
    para.bio.strata(i).stdLen1=std(L1);
    para.bio.strata(i).wgt1=Wgt1;
    
%%  biological sampling station #2
    para.bio.strata(i).length2=L2;
    para.bio.strata(i).gender2=Gender2;
    para.bio.strata(i).weight2=Wgt2;               % weight array of individual fish weight
    para.bio.strata(i).wgt2=nansum(Wgt2);
    para.bio.strata(i).age2=Age2;

    if str2num(para.survey_year) >= 2003  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   using station 1 data only for TS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L=L1;
        Gender=Gender1;
    else
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   using both station 1 and 2 data for TS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       L=[L1 L2];
       Gender=[Gender1 Gender2];
    end
   
% normalized within each trawl
    para.bio.strata(i).sig_bs=nanmean(para.bio.strata(i).sig_bs_j);
    para.bio.strata(i).sig_b=4*pi*para.bio.strata(i).sig_bs;
    para.bio.strata(i).gender=Gender;
    para.bio.strata(i).L=L;
    para.bio.strata(i).Lave=nanmean(para.bio.strata(i).Lave_j);
    para.bio.strata(i).Male_ind=find(para.bio.strata(i).gender == 1); 
    para.bio.strata(i).Female_ind=find(para.bio.strata(i).gender == 2);
    para.bio.strata(i).Gender(para.bio.strata(i).Male_ind)=1;
    para.bio.strata(i).Gender(para.bio.strata(i).Female_ind)=2;
    para.bio.strata(i).LaveM=nanmean(para.bio.strata(i).L(para.bio.strata(i).Male_ind));
    para.bio.strata(i).LaveF=nanmean(para.bio.strata(i).L(para.bio.strata(i).Female_ind));
    para.bio.strata(i).nM=length(para.bio.strata(i).Male_ind);
    para.bio.strata(i).nF=length(para.bio.strata(i).Female_ind);
    %% quantities associated with sample station #
    para.bio.strata(i).n=length(L1);
    para.bio.strata(i).L1=L1;
    male_indx=find(para.bio.strata(i).gender1 == 1);
    female_indx=find(para.bio.strata(i).gender1 == 2);
    para.bio.strata(i).nM1=length(male_indx);
    para.bio.strata(i).nF1=length(female_indx);
    L1_M=L1(male_indx);
    L1_F=L1(female_indx);
    male_wgt=nansum(interp1(para.bio.hake_len_bin,data.bio.len_wgt_M,L1_M));
    female_wgt=nansum(interp1(para.bio.hake_len_bin,data.bio.len_wgt_F,L1_F));
    para.bio.strata(i).nM_wgt1=Wgt1*male_wgt/(male_wgt+female_wgt);
    para.bio.strata(i).nF_wgt1=Wgt1*female_wgt/(male_wgt+female_wgt);
    
    para.bio.strata(i).n1=length(L1);
    if para.bio.strata(i).nM+para.bio.strata(i).nF ~= para.bio.strata(i).n
        fprintf('Station 1: stratum: %d\t Unsexed samples = %f\n',i, 1-(para.bio.strata(i).nM+para.bio.strata(i).nF)/para.bio.strata(i).n);
    end
    
% create sex-length key (station 1 samples)
    if ~isempty(Gender1)
%%%%%%%%% commented out on 10-21-2011  but put it back in on 10-25-2011 to make this part only fo rsample station #1 %%%%%%%%%%%%%%%%%%
        LenM_ind=find(para.bio.strata(i).gender1 == 1);
        LenF_ind=find(para.bio.strata(i).gender1 == 2);
        para.bio.strata(i).Len_key_Mx=hist(para.bio.strata(i).length1(LenM_ind),para.bio.hake_len_bin);
        para.bio.strata(i).Len_key_Fx=hist(para.bio.strata(i).length1(LenF_ind),para.bio.hake_len_bin);
        para.bio.strata(i).Len_key_Nx=hist(para.bio.strata(i).length1,para.bio.hake_len_bin);
%%%%%%%%% commented out on 10-21-2011   %%%%%%%%%%%%%%%%%%

% %% modified on 10-21-2011
%         LenM_ind=find(para.bio.strata(i).Gender == 1);
%         LenF_ind=find(para.bio.strata(i).Gender == 2);
%         para.bio.strata(i).Len_key_Mx=hist(para.bio.strata(i).L(LenM_ind),para.bio.hake_len_bin);
%         para.bio.strata(i).Len_key_Fx=hist(para.bio.strata(i).L(LenF_ind),para.bio.hake_len_bin);
%         para.bio.strata(i).Len_key_Nx=hist(para.bio.strata(i).L,para.bio.hake_len_bin);
% %%%% end of modification made on 10-21-2011

        if length(LenM_ind) > 0
            para.bio.strata(i).Len_key_Mn=para.bio.strata(i).Len_key_Mx/length(LenM_ind);
        else
            para.bio.strata(i).Len_key_Mn=0;
        end
        if length(LenF_ind) > 0
            para.bio.strata(i).Len_key_Fn=para.bio.strata(i).Len_key_Fx/length(LenF_ind);
        else
            para.bio.strata(i).Len_key_Fn=0;
        end
%        para.bio.strata(i).Len_key_Nn=para.bio.strata(i).Len_key_Nx/(length(LenM_ind)+length(LenF_ind));     
        para.bio.strata(i).Len_key_Nn=para.bio.strata(i).Len_key_Nx/para.bio.strata(i).n1;     % modified on 10/25/2010
    else
        fprintf('---- no sexed data in stratum  %d\n',i);
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
            para.bio.strata(i).Len_Age_key_M(k,j)=length(L2jkM_ind);
            % weight key
            para.bio.strata(i).Len_Age_key_wgt_M(k,j)=sum(Wgt2(L2jkM_ind));
            indx2=union(indx2,L2jkM_ind);
        % Female
            if k == length(para.bio.hake_len_bin)
               L2jkF_ind=find(Gender2 == 2 & Age2 == para.bio.hake_age_bin(j) & L2 >= para.bio.hake_len_bin(k));
            else
               L2jkF_ind=find(Gender2 == 2 & Age2 == para.bio.hake_age_bin(j) & (L2 >= para.bio.hake_len_bin(k) & L2 < para.bio.hake_len_bin(k+1)));
            end
            % number key
            para.bio.strata(i).Len_Age_key_F(k,j)=length(L2jkF_ind);
             % weight key
            para.bio.strata(i).Len_Age_key_wgt_F(k,j)=sum(Wgt2(L2jkF_ind));
            indx2=union(indx2,L2jkF_ind);
       % ALL 
             % number key
           para.bio.strata(i).Len_Age_key_ALL(k,j)=para.bio.strata(i).Len_Age_key_M(k,j)+para.bio.strata(i).Len_Age_key_F(k,j);
             % weight key
           para.bio.strata(i).Len_Age_key_wgt_ALL(k,j)=para.bio.strata(i).Len_Age_key_wgt_M(k,j)+para.bio.strata(i).Len_Age_key_wgt_F(k,j);
        end
        if  isnan(para.bio.strata(i).Len_Age_key_wgt_M(k,j))
            disp([i,j,k])
        end
    end
    total_matrix_wgt2=nansum(nansum(para.bio.strata(i).Len_Age_key_wgt_ALL));
    total_wgt2=sum(Wgt2);
    if abs(total_matrix_wgt2-total_wgt2) > 1e-2
        fprintf('strata = %d\t  Matrix wgt2 = %f\t  wgt2 = %f\n',i,total_matrix_wgt2,total_wgt2);
    end
 %% normalized the length-age matrices
 %% length-age number keys
    para.bio.strata(i).matrix_NtotM=nansum(nansum(para.bio.strata(i).Len_Age_key_M));
    para.bio.strata(i).matrix_NtotF=nansum(nansum(para.bio.strata(i).Len_Age_key_F));
    para.bio.strata(i).matrix_NtotALL=para.bio.strata(i).matrix_NtotM+para.bio.strata(i).matrix_NtotF;
    para.bio.strata(i).Len_Age_key_Mn=para.bio.strata(i).Len_Age_key_M/para.bio.strata(i).matrix_NtotM;
    para.bio.strata(i).Len_Age_key_Fn=para.bio.strata(i).Len_Age_key_F/para.bio.strata(i).matrix_NtotF;    
    para.bio.strata(i).Len_Age_key_ALLn=para.bio.strata(i).Len_Age_key_ALL/para.bio.strata(i).matrix_NtotALL;
 %% length-age weight keys
    para.bio.strata(i).matrix_NtotM_wgt=max(eps,nansum(nansum(para.bio.strata(i).Len_Age_key_wgt_M)));
    para.bio.strata(i).matrix_NtotF_wgt=max(eps,nansum(nansum(para.bio.strata(i).Len_Age_key_wgt_F)));
    para.bio.strata(i).matrix_NtotALL_wgt=para.bio.strata(i).matrix_NtotM_wgt+para.bio.strata(i).matrix_NtotF_wgt;
    para.bio.strata(i).Len_Age_key_wgt_Mn=para.bio.strata(i).Len_Age_key_wgt_M/para.bio.strata(i).matrix_NtotM_wgt;
    para.bio.strata(i).Len_Age_key_wgt_Fn=para.bio.strata(i).Len_Age_key_wgt_F/para.bio.strata(i).matrix_NtotF_wgt;    
    para.bio.strata(i).Len_Age_key_wgt_ALLn=para.bio.strata(i).Len_Age_key_wgt_ALL/para.bio.strata(i).matrix_NtotALL_wgt;
%     if para.bio.strata(i).matrix_NtotALL ~= length(L2)
%         fprintf('Station 2: stratum: %d\t Unsexed samples ratio = %f\n',i, 1-(para.bio.strata(i).matrix_NtotM+para.bio.strata(i).matrix_NtotF)/length(L2));
%     end

%     figure(1);imagesc(para.bio.strata(i).Len_Age_Mn);colorbar;title(sprintf('stratum = %d',i))
%     figure(2);imagesc(para.bio.strata(i).Len_Age_Fn);colorbar;title(sprintf('stratum = %d',i));pause(1)

   para.bio.strata(i).total_N=para.bio.strata(i).n+para.bio.strata(i).matrix_NtotALL;   % station #1 + station #2
   para.bio.strata(i).total_wgt=para.bio.strata(i).wgt1+para.bio.strata(i).wgt2;   % station #1 + station #2
   % len and len-Age number key proportions
   para.bio.strata(i).Len_Age_M_proportion=para.bio.strata(i).matrix_NtotM/para.bio.strata(i).total_N;
   para.bio.strata(i).Len_Age_F_proportion=para.bio.strata(i).matrix_NtotF/para.bio.strata(i).total_N;
   para.bio.strata(i).Len_M_proportion=para.bio.strata(i).nM1/para.bio.strata(i).total_N;
   para.bio.strata(i).Len_F_proportion=para.bio.strata(i).nF1/para.bio.strata(i).total_N;
   
   % len and len-Age weight key proprotions
   para.bio.strata(i).Len_Age_M_wgt_proportion=para.bio.strata(i).matrix_NtotM_wgt/para.bio.strata(i).total_wgt;
   para.bio.strata(i).Len_Age_F_wgt_proportion=para.bio.strata(i).matrix_NtotF_wgt/para.bio.strata(i).total_wgt;
   para.bio.strata(i).Len_M_wgt_proportion=para.bio.strata(i).nM_wgt1/para.bio.strata(i).total_wgt;
   para.bio.strata(i).Len_F_wgt_proportion=para.bio.strata(i).nF_wgt1/para.bio.strata(i).total_wgt;
   
   check1=para.bio.strata(i).Len_Age_M_proportion+para.bio.strata(i).Len_Age_F_proportion+para.bio.strata(i).Len_M_proportion+para.bio.strata(i).Len_F_proportion;
   check2=para.bio.strata(i).Len_Age_M_wgt_proportion+para.bio.strata(i).Len_Age_F_wgt_proportion+para.bio.strata(i).Len_M_wgt_proportion+para.bio.strata(i).Len_F_wgt_proportion;
   if abs(check1+check2 - 2) > 0.001
       fprintf('un-sexed samples are found in stratum: %d\t number check = %f\t weight check = %f\n',i,check1,check2);
   end
   %% generate len-wgt-key for each stratum (only use the data from sample sation #2)
   indx2=[];
   for k=1:length(para.bio.hake_len_bin)
       % Male
       if k == length(para.bio.hake_len_bin)
           L2kM_ind=find(Gender2 == 1 & para.bio.strata(i).length2 >= para.bio.hake_len_bin(k));
       else
           L2kM_ind=find(Gender2 == 1 & (para.bio.strata(i).length2 >= para.bio.hake_len_bin(k) & para.bio.strata(i).length2 < para.bio.hake_len_bin(k+1)));
       end
       if length(L2kM_ind) >= 5
            para.bio.strata(i).Len_wgt_key_M(k)=nanmean(para.bio.strata(i).weight2(L2kM_ind));
       else
           para.bio.strata(i).Len_wgt_key_M(k)=data.bio.len_wgt_all.reg_w0M*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_pM;
       end
  %     para.bio.strata(i).Len_wgt_key_M_ind(k)=length(L2kM_ind);
       % Female
       if k == length(para.bio.hake_len_bin)
           L2kF_ind=find(Gender2 == 2 & para.bio.strata(i).length2 >= para.bio.hake_len_bin(k));
       else
           L2kF_ind=find(Gender2 == 2 & (para.bio.strata(i).length2 >= para.bio.hake_len_bin(k) & para.bio.strata(i).length2 < para.bio.hake_len_bin(k+1)));
       end
       if length(L2kF_ind) >= 5
            para.bio.strata(i).Len_wgt_key_F(k)=nanmean(para.bio.strata(i).weight2(L2kF_ind));
       else
           para.bio.strata(i).Len_wgt_key_F(k)=data.bio.len_wgt_all.reg_w0F*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_pF;
       end
       para.bio.strata(i).Len_wgt_key_F_ind(k)=length(L2kF_ind);
       % ALL
       ALL_nk=length(L2kM_ind)+length(L2kF_ind);
       if length(ALL_nk) >= 5
           para.bio.strata(i).Len_wgt_key_ALL(k)=(para.bio.strata(i).Len_wgt_key_M(k)*length(L2kM_ind)+para.bio.strata(i).Len_wgt_key_F(k)*length(L2kF_ind))/ALL_nk;
       else
           para.bio.strata(i).Len_wgt_key_ALL(k)=data.bio.len_wgt_all.reg_w0*para.bio.hake_len_bin(k)^data.bio.len_wgt_all.reg_p;
       end       
       para.bio.strata(i).Len_wgt_key_ALL_ind(k)=ALL_nk;       
   end
end                             % end of strata loop

return