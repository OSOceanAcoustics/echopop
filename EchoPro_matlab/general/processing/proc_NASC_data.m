function proc_NASC_data
% process echo-integration NASC data to obtain length-age-sex structured biomass estimate and all required and useful tables & results
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Modification:       12/12/2013
%% Modification:       11/17/2016   handle 1995 data since the NASC export has 13 columns but no haul assignment (old export format)
%% Modification:       12/10/2019   including flexibility of using observer biological data
%% Modification:       03/10/2020   transect reduction will include trawl reduction
%% Modification:       05/29/2021   1995 VL values were missing --> interval --> NaN, change interval = 0.5 nmi (lines 169-171)

global para data hdl  

tic
fprintf('\n\n================== Biomass Estimate of %s Data ===================\n\n',para.survey_year)

fprintf('Start Time = %s\n', datetime)
%% delete the kriging information popup window
if isfield(hdl,'krig') && isfield(hdl.krig,'msg')
    if ishandle(hdl.krig.msg)
       delete(hdl.krig.msg)
    end
end

%% make processing parameters consistent with whatever shown on the window
para=check_GUI_para_consistency(hdl,para);

%% load & process biological (trawl) data
data_in_filename=data.in.filename;
data=[];
data.in.filename=data_in_filename;
data_in=data.in;
data.final.table=[];
data.in=data_in;
para.bio.strata=[];
data.bio=[];
if para.bio_data_type ~= 1 | (para.bio_data_type == 1 & ~strcmp(para.survey_year, '2011'))   % for acoustic trawls in 2011, there are age1 trawls need to be excluded if chosen
    para.proc.age1_haul = [];
end

% para.bio_data_type=1;
fprintf('Load Biological trawl Data ...\n')
    switch para.bio_data_type
        case 1  % acoustic trawl
            %% load trawl biological data
            get_combined_biological_data;
% %             %% process biological data
% %             % get length-weight key (an array of weight per fish at legnth)
% %             proc_len_wgt_data;
% %             % construct the trawl information table
% %             construct_catch_trawl_output_matrices
        case 2    % trawl survey
        case 3   % observer
            %% load trawl biological data
            get_combined_observer_biological_data;
% %             %% process biological data
% %             % get length-weight key (an array of weight per fish at legnth)
% %             proc_len_wgt_data;
% %             % construct the trawl information table
% %             construct_observer_catch_trawl_output_matrices
    end


%% Biomass Estimate Data Processing
while para.proc.bootstrap.cnt < para.proc.bootstrap.limit
    para.proc.bootstrap.cnt=para.proc.bootstrap.cnt+1; 
    if para.proc.bootstrap.limit > 1 & para.proc.bootstrap.cnt <= para.proc.bootstrap.limit
        fprintf('\n\n************ Bootstrapping Count = %d *********** \n\n',para.proc.bootstrap.cnt);
    end
    
    if isempty(data.bio.hake_length_weight_sex_age)
        f = errordlg('No Age Data. de-select Age Data!', 'Error Dialog', 'modal');
        return
    end
    %% load and process echo sounder data
    %% create stratum-based length-age and length-weight matirices with gender structure
    if para.proc.bootstrap.cnt == 1             % first bootstrap run
        data0 = data;                           % save original data struct
        data.final.bootstrap.out=[];
        data.proc.start_time=clock;
        %% load VL-interval based NASC data that were exported from EchoView previously
        VL_BiomassInt=get_NASC_data(para.acoust.filename.processed_data,para.survey_year,para.proc.transect_offset);
        if para.proc.ST_BT_RT_ET_zero_removal_flag == 1  % remove zeros before ST, after ET, and between BT & RT along transect line
            VL_BiomassInt=remove_unwanted_zeros(para.proc.transect_info_filename,VL_BiomassInt);
        end
        data.tmp.VL_BiomassInt=VL_BiomassInt;
        para.proc.transect_reduction_fraction=str2num(get(hdl.proc.edit_transect_reduction_fraction,'string'));
    else
        VL_BiomassInt=data.tmp.VL_BiomassInt;
    end
    if para.proc.transect_reduction_fraction ~= 0
        [TX_selected, selected_TX_reg_ind]=get_resampled_sonar_data(data,para,VL_BiomassInt);
        %% update VL_BiomassInt
        VL_BiomassInt=VL_BiomassInt(selected_TX_reg_ind,:);     
        %% update len_wgt_all
        update_data_bio_struct(TX_selected, data0);
        if size(VL_BiomassInt,1) > size(data.tmp.VL_BiomassInt,1)
            disp('----------------- wrong reduction!! ---------------------')
        end
    end
    
    switch para.bio_data_type
        case 1  % acoustic trawl
            %% process biological data
            % get length-weight key (an array of weight per fish at legnth)
            proc_len_wgt_data;
            % construct the trawl information table
            construct_catch_trawl_output_matrices
        case 2    % trawl survey
        case 3   % observer
            %% process biological data
            % get length-weight key (an array of weight per fish at legnth)
            proc_len_wgt_data;
            % construct the trawl information table
            construct_observer_catch_trawl_output_matrices
    end
    data.bio.strata=[];
    get_historical_strata_data;
    transect_num=VL_BiomassInt(:,1);
    ind=find(transect_num >= para.proc.start_transect & transect_num <= para.proc.end_transect);
    VL_BiomassInt=VL_BiomassInt(ind,:);
    transect_num=VL_BiomassInt(:,1);
    NASC_int=VL_BiomassInt(:,12);
    transect_spacing = VL_BiomassInt(:,8);
    %% change NaN transect spacing to 10 nmi for 2001
    if ~isempty(strmatch(para.survey_year,{'2001'})) | ~isempty(strmatch(para.survey_year,{'1995'}))
        ind=find(isnan(transect_spacing) ==1);
        transect_spacing(ind)=10;
    end
    %% get stratification information
    haul_stratum_id=get_geographic_strata(VL_BiomassInt(:,5),para.proc.stratification_filename,para.proc.stratification_index+1);
    para.bio.strata=data.bio.strata;
    %         haul_stratum_id=VL_BiomassInt(:,7);
    layer_dep=VL_BiomassInt(:,9);
    layer_height=VL_BiomassInt(:,10);
    bot_dep=VL_BiomassInt(:,11);
    
    if str2num(char(para.survey_year)) < 2003
        layer_dep= nan;
        layer_height= nan;
        hauls=nan;
    else  % later than 2003
        if  size(VL_BiomassInt,2) > 12
            if para.bio_data_type == 3   % only for >= 2003
                % observer trawl data
                hauls = get_observer_hauls(VL_BiomassInt, para.proc.stratification_filename);
                VL_BiomassInt(:,13) = hauls;
            elseif para.bio_data_type == 1
                % acoustic survey trawls
                hauls = VL_BiomassInt(:,13);
            end
        end
    end
    
    ind = find(haul_stratum_id > length(data.bio.strata));
    if ~isempty(ind)
        haul_stratum_id(ind) = length(data.bio.strata);
    end
    
        
    n=size(VL_BiomassInt,1);
    if isfield(para, 'platform_name')    % sail drone has no Vessel Log distance
        if strcmp(para.platform_name, 'SD')  % SD data
            interval = 0.5*ones(n,1);     % nominal interval distance 0.5 nm based on GPS log distance
        else  % FSV Vessel Log distance in nm
            interval(1:n-1) = VL_BiomassInt(2:n,3)-VL_BiomassInt(1:n-1,3);
            interval(n)=VL_BiomassInt(n,4)-VL_BiomassInt(n,3);
        end
    else     % FSV Vessel Log distance in nm
        interval(1:n-1) = VL_BiomassInt(2:n,3)-VL_BiomassInt(1:n-1,3);
        interval(n)=VL_BiomassInt(n,4)-VL_BiomassInt(n,3);
    end
    %% fix NaN for VL of 1995 (column 2)
    if ~isempty(strmatch(para.survey_year,{'1995'}))
       nan_VL_ind = find(isnan(interval) == 1);
       interval(nan_VL_ind) = 0.5;   % set to 0.5 nmi
    end
    interval(n+1:end)=[];
    
    median_VL=nanmedian(interval);
    ind=find(abs(interval-median_VL) > 0.05);  % remove outliers at the end of transect
    interval(ind)=VL_BiomassInt(ind,4)-VL_BiomassInt(ind,3);
    % the following two lines were added on 12/4/2013 to remove some
    % outliers in VL data
    ind=find(abs(interval)> 1);
    interval(ind)=median_VL;
    
    %%%%%%%% load updated strata table: Transect-Region-Haul --> re-assigne haul-based region strata
    if para.proc.KS_stratification == 1 % | para.proc.stratification_index > 1
        haul_stratum_id=ones(size(haul_stratum_id));
        if para.proc.stratification_index == 10
            haul_strata=xlsread(para.acoust.filename.strata,'length strata byhaul_1stratum');
        else
            haul_strata=xlsread(para.acoust.filename.strata,para.proc.stratification_index+1);
        end
        %% Regions with Reg Index = 999 or no haul assignment will be used
        ind=find(VL_BiomassInt(:,2) < 900 & hauls > 0);
        
        if size(VL_BiomassInt,2) > 12 & strcmp(para.survey_year, '1995') ~= 1  % newer exported format
            haul_num=hauls(ind);

            for i=1:length(ind)
%                 if haul_num(i) == 42
%                     disp([i haul_num(i)])
%                 end
%                 if VL_BiomassInt(ind(i),1) == 38 & VL_BiomassInt(ind(i),2) == 45
%                     disp([i ind(i)])
%                 end
                stratum_indx=find(haul_num(i) == haul_strata(:,4));
                if ~isempty(stratum_indx)
                    haul_stratum_id(ind(i))=haul_strata(stratum_indx,2);
                end
            end
        else
            if ~isempty(para.bio_acoust.filename.Transect_region_haul)
                T_reg_haul=xlsread(para.bio_acoust.filename.Transect_region_haul);
                for i=1:length(T_reg_haul)
                    indx=find(VL_BiomassInt(:,1)== T_reg_haul(i,1) & VL_BiomassInt(:,2) == T_reg_haul(i,2));
                    indx1=find(T_reg_haul(i,3) == haul_strata(:,4));
                    haul_stratum_id(indx)=haul_strata(indx1,4);
                end
            else
                haul_stratum_id=VL_BiomassInt(:,7);
            end
        end
    end
    %%%%%%%%  end of re-assigne strata
    fprintf('Process un-kriged Acoustic & Biological Data ...\n')
    
    %% create acoustically weighted length-age-sex structured output matrices
    Final_out.Len_Age_Matrix_AcoustM=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for male hake
    Final_out.Len_Age_Matrix_AcoustF=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for female hake
    Final_out.Len_Age_Matrix_AcoustALL=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for all hake
    Final_out.Len_Age_AcoustM_age1_from_unaged=zeros(length(para.bio.hake_len_bin),1);
    Final_out.Len_Age_AcoustF_age1_from_unaged=zeros(length(para.bio.hake_len_bin),1);
    Final_out.Len_Age_AcoustALL_age1_from_unaged=zeros(length(para.bio.hake_len_bin),1);
    
    
    %% initialize weight-key as a function of length and age
    Wgt_len_age=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
    Final_Wgt_len_age_M=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
    Final_Wgt_len_age_F=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
    Final_Wgt_len_age_ALL=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
    
    %% initialize other variables
    ntk_male=[]; ntk_female=[]; ntk=[];  Wgt_male_int=[]; Wgt_female_int=[]; Wgt_ALL_int=[];
    nntk_male=[]; nntk_female=[]; nntk=[];  nWgt_male_int=[]; nWgt_female_int=[]; nWgt_ALL_int=[];
    mix_sa_ratio=[];sig_b=[];len_wgt_ALL_ii=[];
    
    nInt=size(VL_BiomassInt,1);
    Wgt_ALL_int_mix=0;
    
    ntk_age1=0;
    nasc_age1=0;
    wgt_age1=0;
    for ii=1:nInt
        if rem(ii,500) == 0
            fprintf('*******************************  Processing VL Int %d **************************************\n',ii);
        end
        if haul_stratum_id(ii) ~= 0 & ~isnan(haul_stratum_id(ii))
            if data.bio.strata(haul_stratum_id(ii)).LaveALL < para.bio.Lmin  % exclude YOY hake
                NASC_int(ii)=0;
            end
            sig_b(ii)=data.bio.strata(haul_stratum_id(ii)).sig_b;   % 4*pi* diff. backscattering cross section sig_bs
%                 if haul_stratum_id(ii) == 10
%                     disp([ii VL_BiomassInt(ii,1) VL_BiomassInt(ii,2) haul_stratum_id(ii) ])
%                 end
        else
            sig_b(ii)=1e10;
            haul_stratum_id(ii)=nan;
            NASC_int(ii)=0;
        end
%         if ~isnan(NASC_int(ii)) & NASC_int(ii) > eps %& transect_num(ii) > 85% & haul_stratum_id(ii) == 9
%             fprintf('Index = %d\t Stratum = %d\t NASC_int = %6.2f\n',ii, haul_stratum_id(ii),NASC_int(ii))
%         end
        mix_sa_ratio(ii)=1;
        %% find proportion coefficient for mix species
        if size(VL_BiomassInt,2) > 12  & strcmp(para.survey_year, '1995') ~= 1 
            haul_ind=find(data.bio_acoust.haul_wgt_tbl(:,1) == hauls(ii));
            if haul_ind ~= 0
                mix_sa_ratio(ii)=data.bio_acoust.haul_wgt_tbl(haul_ind,2);
            else
                mix_sa_ratio(ii)=0;
            end
        else
            if nanmean(data.bio.strata(haul_stratum_id(ii)).wgt) < 0.99
                mix_sa_ratio(ii)=nanmean(data.bio.strata(haul_stratum_id(ii)).wgt);
            end
        end
%         if ii > 10615
%             disp(ii)
%         end
        if ~isnan(haul_stratum_id(ii)) % & haul_stratum_id(ii) > 0
            if NASC_int(ii) >= 1   & interval(ii) > 1
                fprintf('Interval > 1 & NASC_{int} > 1.0 nm: indx = %d\t interval = %3.1f\t NASC_{int} = %10.2f\n', ii, interval(ii),NASC_int(ii))
            end
            
            if transect_spacing(ii) ~= 0
                nntk(ii)=round(mix_sa_ratio(ii)*NASC_int(ii)/sig_b(ii));     % normalized # hake/colume within the interval and spacing #/area
            else
                nntk(ii)=0;     % normalized # hake/colume within the interval and spacing #/area
            end
            ntk(ii)=(mix_sa_ratio(ii)*interval(ii)*transect_spacing(ii)*NASC_int(ii)/sig_b(ii));               % total # hake/colume within the interval and spacing
            if para.proc.bootstrap.limit > 1 & para.proc.transect_reduction_fraction ~= 0 & ntk(ii) ~= 0
                ntk(ii) = ntk(ii)/(1-para.proc.transect_reduction_fraction/100);
            end
            ntk1(ii)=ntk(ii);                   % age1 is included
%             if ntk(ii) > 0 & nntk(ii) > 0
%                 disp(ii)
%             end
            if para.proc.exclude_age1 == 1
                %% number density
                age1_num_proportion=nansum(data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn(:,1));
                %% NASC
                sig_bs_aged_ave=nansum((para.acoust.sig_b_coef*para.bio.hake_len_bin.^2)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn);
                age1_nasc_proportion=(para.acoust.sig_b_coef*para.bio.hake_len_bin.^2)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn(:,1)/sig_bs_aged_ave;
                if sig_bs_aged_ave == 0
                    fprintf('sigma_{bs} = 0:  indx = %d\t, stratum = %d\t, sigma_{bs} = %8.2f\n', ii, haul_stratum_id(ii), sig_bs_aged_ave)
                end
                %% weight or biomass density
                % length distribution of age-1 fish
                indx1=find(para.bio.age1_min_len == para.bio.hake_len_bin);
                len_dist_age1=data.bio.strata(haul_stratum_id(ii)).Len_key_Nn(indx1:end).*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn(indx1:end,1)';
                if nansum(len_dist_age1) > 1e-10  | nansum(data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn(indx1:end,1)) > 1e-10
                    % modified on 12/30/2015 to handle and situation that
                    % age1 mixed
                    age1_wgt_proportion=nansum(data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn(:,1))/nansum(nansum(data.bio.strata(haul_stratum_id(ii)).Len_Age_key_ALLn));
                else
                    age1_wgt_proportion=0;
                end
            else        % age-1 will not be excluded
                age1_num_proportion=0;
                age1_nasc_proportion=0;
                age1_wgt_proportion=0;
            end
            age2y_num_proportion=1-age1_num_proportion;
            age2y_nasc_proportion=1-age1_nasc_proportion;
            age2y_wgt_proportion=1-age1_wgt_proportion;
            
            
            if ~isnan(ntk(ii))
                %% proportion the numbers+
                %% Length station (#1)
                N_len_M=ntk(ii)*data.bio.strata(haul_stratum_id(ii)).Len_M_proportion;
                N_len_F=ntk(ii)*data.bio.strata(haul_stratum_id(ii)).Len_F_proportion;
                %% length and age station (#2)
                N_len_age_M=ntk(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion;
                N_len_age_F=ntk(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion;
                
                %% update Len-Age matrix for unaged hake
                Len_M_ii=(N_len_M*data.bio.strata(haul_stratum_id(ii)).Len_key_Mn(:));
                Len_F_ii=(N_len_F*data.bio.strata(haul_stratum_id(ii)).Len_key_Fn(:));
                Final_out.Len_Age_Matrix_AcoustM(:,end)=Final_out.Len_Age_Matrix_AcoustM(:,end)+Len_M_ii;
                Final_out.Len_Age_Matrix_AcoustF(:,end)=Final_out.Len_Age_Matrix_AcoustF(:,end)+Len_F_ii;
                Final_out.Len_Age_Matrix_AcoustALL(:,end)=Final_out.Len_Age_Matrix_AcoustALL(:,end)+ (Len_M_ii + Len_F_ii);
                %% update Len-Age matrix for aged hake
                Len_Age_M_ii=nan*ones(1,length(para.bio.hake_age_bin));
                Len_Age_F_ii=nan*ones(1,length(para.bio.hake_age_bin));
                if N_len_age_M ~=0
                    %% add male number to the matrix
                    Len_Age_M_ii=(N_len_age_M*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_Mn(:,1:end));
                    Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)+Len_Age_M_ii;
                    Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)+Len_Age_M_ii;
                end
                if N_len_age_F ~=0
                    %% add female number to the matrix
                    Len_Age_F_ii=(N_len_age_F*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_Fn(:,1:end));
                    Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)+Len_Age_F_ii;
                    Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)+Len_Age_F_ii;
                end
            end
            [len_wgt_M, len_wgt_F, len_wgt_ALL]=get_wgt_len_key(haul_stratum_id(ii));
            
            %%  sum of the proportions of abundance of FSCS station1 (length-sex) and FSCS station 2 (length-weight-sex-age)
            ntk_male(ii)=(ntk(ii)*(data.bio.strata(haul_stratum_id(ii)).Len_M_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion));
            ntk_female(ii)=(ntk(ii)*(data.bio.strata(haul_stratum_id(ii)).Len_F_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion));
            var1=(data.bio.strata(haul_stratum_id(ii)).Len_M_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion);
            var2=(data.bio.strata(haul_stratum_id(ii)).Len_F_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion);
            var=var1+var2;
%             if abs(ntk_male(ii)+ntk_female(ii) - (ntk(ii))) > 1e-2 
%                 fprintf('ii = %d\t, stratum # = %d\t var = %10.5f\t m + f = %10.1f\t n_tot = %10.1f\n',ii,haul_stratum_id(ii),var,ntk_male(ii)+ntk_female(ii),ntk(ii))
%             end
            
            %% convert abundance (numbers) to biomass (weight)
            Wgt_male_int(ii)=nansum(ntk_male(ii)*len_wgt_M);                                % in kg
            Wgt_female_int(ii)=nansum(ntk_female(ii)*len_wgt_F);                            % in kg
            Wgt_unsexed_int(ii)=nansum((ntk(ii)-ntk_male(ii)-ntk_female(ii))*len_wgt_ALL);  % in kg for unsexed hake
            Wgt_ALL_int(ii)=Wgt_male_int(ii)+Wgt_female_int(ii)+Wgt_unsexed_int(ii);
            len_wgt_ALL_ii(ii)=len_wgt_ALL;
            
            if ~isnan(ntk(ii))
                if ~isnan(sum(sum(data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn)))
                    %% 2D Matrix for both aged and unaged fish
                    Wgt_len_age=Wgt_len_age+Wgt_ALL_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
                    
                    %% acumulative 2D Matrix for both aged and unaged biomass
                    Final_Wgt_len_age_M=Final_Wgt_len_age_M+Wgt_male_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
                    Final_Wgt_len_age_F=Final_Wgt_len_age_F+Wgt_female_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
                    Final_Wgt_len_age_ALL=Final_Wgt_len_age_ALL+Wgt_ALL_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
 %% Check age-1
%                     age1_array=[para.bio.hake_len_bin(12:16)' Final_out.Len_Age_Matrix_AcoustALL(12:16,1)  Final_Wgt_len_age_ALL(12:16,1)];
%                     age1_array=[age1_array Final_Wgt_len_age_ALL(12:16,1)./(Final_Wgt_len_age_ALL(12:16,1)+eps)];
%                     if nanmean(age1_array(:,4)) > 0
%                         disp(age1_array)
%                     end
                end
                
                if mix_sa_ratio(ii) < 0.99
%                     if isnan(Wgt_ALL_int(ii))
%                         disp(Wgt_ALL_int(ii))
%                     end
                    Wgt_ALL_int_mix=Wgt_ALL_int_mix+nansum(Wgt_ALL_int(ii));
                end
            end
            
            %%%%%% normalized quantities
            nntk_male(ii)=round(nntk(ii)*(data.bio.strata(haul_stratum_id(ii)).Len_M_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion));
            nntk_female(ii)=round(nntk(ii)*(data.bio.strata(haul_stratum_id(ii)).Len_F_proportion+data.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion));
            %% convert the averaged and normalized number of fish to the averaged and normalized weight in kg/nm^2
            nWgt_male_int(ii)=nntk_male(ii)*len_wgt_M;                                  % in kg
            nWgt_female_int(ii)=nntk_female(ii)*len_wgt_F;                              % in kg
            nWgt_unsexed_int(ii)=(nntk(ii)-nntk_male(ii)-nntk_female(ii))*len_wgt_ALL;  % in kg for unsexed hake
            nWgt_ALL_int(ii)=nWgt_male_int(ii)+nWgt_female_int(ii)+nWgt_unsexed_int(ii);
            
%             if ii == 2213
%                 disp([ii age2y_nasc_proportion])
%             end
            %% exclude age-1 portion of the kriging input quantitie: abundance (number) density and biomass (weight) density
            nntk(ii)=nntk(ii)*age2y_num_proportion;
            NASC_int(ii)=NASC_int(ii)*age2y_nasc_proportion;
            nWgt_ALL_int(ii)=nWgt_ALL_int(ii)*age2y_wgt_proportion;
            ntk0(ii)=ntk(ii);
            
            if isnan(NASC_int(ii))
                fprintf('indx = %d\t age2+ prop = %3.2f\t NASC_{int} = %6.2f\n', ii, age2y_nasc_proportion, NASC_int(ii))
            end
            %% other parameteres: numbers (abundance) and biomass (weight) of each VL-interval based area (Interval x Transect spacing)
            %% Abundance (number)
            ntk(ii)=ntk(ii)*age2y_num_proportion;
            ntk_male(ii)=ntk_male(ii)*age2y_num_proportion;
            ntk_female(ii)=ntk_female(ii)*age2y_num_proportion;
            if ~isnan(ntk(ii))
                Final_out.Len_Age_AcoustM_age1_from_unaged=Final_out.Len_Age_AcoustM_age1_from_unaged+age1_num_proportion*Len_M_ii;
                Final_out.Len_Age_AcoustF_age1_from_unaged=Final_out.Len_Age_AcoustF_age1_from_unaged+age1_num_proportion*Len_F_ii;
                Final_out.Len_Age_AcoustALL_age1_from_unaged=Final_out.Len_Age_AcoustALL_age1_from_unaged+age1_num_proportion*(Len_M_ii+Len_F_ii);
            end
            %% Biomass (weight)
            Wgt_male_int(ii)=Wgt_male_int(ii)*age2y_wgt_proportion;
            Wgt_female_int(ii)=Wgt_female_int(ii)*age2y_wgt_proportion;
            added_all_age_biomass0=Wgt_ALL_int(ii);
            added_all_age_biomass1=nansum(nansum(Wgt_ALL_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn));
            added_age2_biomass1=nansum(nansum(Wgt_ALL_int(ii)*data.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn(:,2:end)));
            Wgt_ALL_int(ii)=Wgt_ALL_int(ii)*age2y_wgt_proportion;
            added_age2_biomass0=Wgt_ALL_int(ii);
        end  % end of if -  ~isnan(haul_stratum_id(ii)) & haul_stratum_id(ii) > 0
        if abs(1-added_all_age_biomass0/added_all_age_biomass1) > 1e-5
            fprintf('inconsistency in all_age biomass portitioning: %d\t  %8.3f\t %8.3f\n',ii,added_all_age_biomass0,added_all_age_biomass1)
        end
        if abs(1-added_age2_biomass0/added_age2_biomass1) > 1e-5 & para.proc.exclude_age1 == 1
            fprintf('inconsistency in age2+ biomass portitioning: %d\t  %8.3f\t %8.3f\n',ii,added_age2_biomass0,added_age2_biomass1)
        end
        if abs(ntk(ii)) == 0 &  Wgt_ALL_int(ii) > eps
            fprintf('inconsistency in abundance and biomass: %d\t %d\t %10.3f\n',ii,ntk(ii),Wgt_ALL_int(ii));
        end
%         if  ii > 11615
%             disp(num2str([ii  ntk(ii)  Wgt_ALL_int(ii)]))
%         end
        
    end  % end of interval loop
    
    %% construct un-kriged VL-interval based NASC-Biomass table
    if str2num(char(para.survey_year)) < 2003
        %% Columns   1-9:    'Transect'    'Region no.'    'VL_start'    'VL_end'     'Lat'      'Lon'      'stratum'      'Depth'      'NASC'
        
        %% Columns  10-15:   'ntk_male'    'ntk_female'    'ntk_total'   'wgt_male'    'wgt_female'    'wgt_total'          -- absolute values
        
        %% Columns 16-21:    'nntk_male'   'nntk_female'   'nntk_total'  'nwgt_male'   'nwgt_female'    'nwgt_total'        -- normalized values
        
        VL_BiomassInt(:,[7:21])=[haul_stratum_id(:) bot_dep(:) NASC_int(:) ntk_male(:) ntk_female(:) ntk(:)  Wgt_male_int(:) Wgt_female_int(:) Wgt_ALL_int(:) ...
            nntk_male(:) nntk_female(:) nntk(:)  nWgt_male_int(:) nWgt_female_int(:) nWgt_ALL_int(:)];
        
    else
        %% Columns   1-9:    'Transect'    'Region no.'    'VL_start'    'VL_end'     'Lat'      'Lon'      'stratum'      'Depth'      'NASC'
        
        %% Columns  10-15:   'ntk_male'    'ntk_female'    'ntk_total'   'wgt_male'    'wgt_female'    'wgt_total'          -- absolute values
        
        %% Columns 16-21:    'nntk_male'   'nntk_female'   'nntk_total'  'nwgt_male'   'nwgt_female'    'nwgt_total'        -- normalized values
        
        %% Columns 22-26:    'Layer Depth'  'Layer height'    'Spacing'   'Interval'   'mix_sa_ratio'  'sig_b' 'wgt_per_fish'
        VL_BiomassInt(:,[7:28])=[haul_stratum_id(:) bot_dep(:) NASC_int(:) ntk_male(:) ntk_female(:) ntk(:)  Wgt_male_int(:) Wgt_female_int(:) Wgt_ALL_int(:) ...
            nntk_male(:) nntk_female(:) nntk(:)  nWgt_male_int(:) nWgt_female_int(:) nWgt_ALL_int(:) layer_dep(:) layer_height(:) transect_spacing(:) interval(:) mix_sa_ratio(:)  sig_b(:) len_wgt_ALL_ii(:)];
    end
    
    %% find age-1 hake related quantities
    data.final.grand_tot_biomass=nansum(nansum(Final_Wgt_len_age_ALL))*1e-6;
    %% construct final age-1 biomass table from age-1 abundance table
    if para.proc.survey_len_wgt_key_sex == 1
        biomass_age1_M0=(data.bio.len_wgt_M*Final_out.Len_Age_Matrix_AcoustM(:,1));
        biomass_age1_F0=(data.bio.len_wgt_F*Final_out.Len_Age_Matrix_AcoustF(:,1));
        biomass_age1_ALL0=(data.bio.len_wgt_ALL*Final_out.Len_Age_Matrix_AcoustALL(:,1));
    else
        biomass_age1_M0=(data.bio.len_wgt_ALL*Final_out.Len_Age_Matrix_AcoustM(:,1));
        biomass_age1_F0=(data.bio.len_wgt_ALL*Final_out.Len_Age_Matrix_AcoustF(:,1));
        biomass_age1_ALL0=(data.bio.len_wgt_ALL*Final_out.Len_Age_Matrix_AcoustALL(:,1));
    end
    if biomass_age1_ALL0 ~= 0 & ~isnan(biomass_age1_ALL0)
        Mcoef=biomass_age1_M0/biomass_age1_ALL0;
        Fcoef=biomass_age1_F0/biomass_age1_ALL0;
    else
        Mcoef=0;
        Fcoef=0;
    end
    biomass_age1_M=Mcoef*nansum(Final_Wgt_len_age_ALL(:,1))*1e-6;
    biomass_age1_F=Fcoef*nansum(Final_Wgt_len_age_ALL(:,1))*1e-6;
    biomass_age1_ALL=nansum(Final_Wgt_len_age_ALL(:,1))*1e-6;
    biomass_age1_ALL=max(biomass_age1_M+biomass_age1_F,biomass_age1_ALL);
    biomass_age1_mix=biomass_age1_ALL-biomass_age1_M-biomass_age1_F;
    data.final.grand_age1_male_biomass=biomass_age1_M;
    data.final.grand_age1_female_biomass=biomass_age1_F;
    data.final.grand_age1_ALL_biomass=biomass_age1_ALL;
    data.final.grand_biomass_age1_mix=biomass_age1_mix;
    
    %% construct overall length-age-sex structured un-kriged biomass table and quantities
    %% total age-1 biomass
    data.final.grand_tot_biomass1y0=biomass_age1_ALL;
    
    if para.proc.exclude_age1 == 1
        data.final.grand_male_biomass=nansum(nansum(Final_Wgt_len_age_M))*1e-6-biomass_age1_M;
        data.final.grand_female_biomass=nansum(nansum(Final_Wgt_len_age_F))*1e-6-biomass_age1_F;
        data.final.grand_tot_biomass2y=data.final.grand_tot_biomass-biomass_age1_ALL;
        data.final.grand_tot_mix_biomass=Wgt_ALL_int_mix*1e-6-biomass_age1_mix;
    else
        data.final.grand_male_biomass=nansum(Wgt_male_int)*1e-6;
        data.final.grand_female_biomass=nansum(Wgt_female_int)*1e-6;
        data.final.grand_tot_biomass1y=nansum(Wgt_ALL_int)*1e-6;
        data.final.grand_tot_biomass2y=data.final.grand_tot_biomass1y - data.final.grand_tot_biomass1y0;
        data.final.grand_tot_mix_biomass=Wgt_ALL_int_mix*1e-6;
    end
    
    if isfield(data.final,'table') & isfield(data.final.table,'trawl')
        trawl=data.final.table.trawl;
        data.final.table=Final_out;
        data.final.table.trawl=trawl;
    else
        data.final.table=Final_out;
    end
    %% create length-haul (length-age-sex structured) un-kriged tables
    %% find all trawls
    t=data.bio.compact_len_haul_ALL(:);
    unique_haul=unique(t)';
    unique_haul=int64([0 data.bio.unique_haul_no]);
    unique_aged_haul=[0 data.bio.unique_aged_haul_no];
    data.final.table.len_haul_M=[unique_haul; para.bio.hake_len_bin(:) data.bio.len_haul_M];
    data.final.table.len_haul_F=[unique_haul; para.bio.hake_len_bin(:) data.bio.len_haul_F];
    data.final.table.len_haul_ALL=[unique_haul; para.bio.hake_len_bin(:) data.bio.len_haul_ALL];
    data.final.table.aged_len_haul_M=[unique_aged_haul; para.bio.hake_len_bin(:) data.bio.aged_len_haul_M];
    data.final.table.aged_len_haul_F=[unique_aged_haul; para.bio.hake_len_bin(:) data.bio.aged_len_haul_F];
    data.final.table.aged_len_haul_ALL=[unique_aged_haul; para.bio.hake_len_bin(:) data.bio.aged_len_haul_ALL];
    data.final.table.compact_len_haul_M=[para.bio.hake_len_bin(:) data.bio.compact_len_haul_M];
    data.final.table.compact_len_haul_F=[para.bio.hake_len_bin(:) data.bio.compact_len_haul_F];
    data.final.table.compact_len_haul_ALL=[para.bio.hake_len_bin(:) data.bio.compact_len_haul_ALL];
    %% create final length-age-sex structured tables with len and age listed
    %% abundance tables
    num_Male=round(data.final.table.Len_Age_Matrix_AcoustM);
    data.final.table.Len_Age_Matrix_AcoustM=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) num_Male];
    num_Female=round(data.final.table.Len_Age_Matrix_AcoustF);
    data.final.table.Len_Age_Matrix_AcoustF=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) num_Female];
    num_ALL=num_Male+num_Female;
    data.final.table.Len_Age_Matrix_AcoustALL=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) num_ALL];
    
    %% biomass table with header
    data.final.table.Wgt_Len_Age_Matrix_AcoustM=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) Final_Wgt_len_age_M*1e-9];         % kg --> mmt
    data.final.table.Wgt_Len_Age_Matrix_AcoustF=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) Final_Wgt_len_age_F*1e-9];
    data.final.table.Wgt_Len_Age_Matrix_AcoustALL=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) Final_Wgt_len_age_ALL*1e-9];
    data.final.table.Wgt_len_age=[0  para.bio.hake_age_bin ; para.bio.hake_len_bin(:) Final_Wgt_len_age_ALL*1e-6];   % kg -> mt
    
    %% unsexed
    unsex_total_biomass=data.final.grand_tot_biomass2y-(data.final.grand_male_biomass+data.final.grand_female_biomass);
    data.final.grand.unsexed_biomass=nansum(Wgt_unsexed_int)*1e-6;
    
    %% save the un-kriged VL-interval based NASC-Biomass table to the global variable 'data'
    data.final.table.biomass=VL_BiomassInt;
    
    %% data structure item names & display biomass results
    if str2num(char(para.survey_year)) < 2003
        data.final.table.biomass_description={'Transect','Region no.','VL_start','VL_end','Lat','Lon','stratum','Depth','NASC','ntk_male','ntk_female','ntk_total','wgt_male','wgt_female','wgt_total', ...
            'nntk_male','nntk_female','nntk_total','nwgt_male','nwgt_female','nwgt_total'};
    else
        data.final.table.biomass_description={'Transect','Region no.','VL_start','VL_end','Lat','Lon','stratum','Depth','NASC','ntk_male','ntk_female','ntk_total','wgt_male','wgt_female','wgt_total', ...
            'nntk_male','nntk_female','nntk_total','nwgt_male','nwgt_female','nwgt_total','Layer Depth','Layer height','Transect spacing','interval','mix_coef','sig_b','wgt_per_fish'};
        get_fish_school_characteristics
    end
    data.final.table.trawl_description={'Trawl_no','Transect','Lat','Lon','stratum','Depth','Surface T','TD Depth T','Length','Gender','Age','Aged Weight','Layer Depth','Unaged Weight'};
    
    
    fprintf('\n ******************** All Survey Transects (before Kriging):  Estimated Biomass *********** \n **************  Male: %5.3f (KMT)\n ************* Female: %5.3f (KMT)\n ************ Un-sexed: %5.3f (KMT)\n ******* Age-1 Biomass: %5.3f (KMT)\n ************ Hake Mix: %5.3f (KMT) \n\n',...
        data.final.grand_male_biomass,data.final.grand_female_biomass,data.final.grand.unsexed_biomass,data.final.grand_tot_biomass1y0,data.final.grand_tot_mix_biomass)
    
    fprintf('----------------------------------------\n');
    if para.proc.exclude_age1 == 1
        fprintf('     Total (Age-2+): %8.3f (KMT)  \n\n',data.final.grand_tot_biomass2y);
    else
        fprintf('     Total Biomass (Age1+): %8.3f (KMT)  \n\n',data.final.grand_tot_biomass1y);
        fprintf('     Total Biomass (Age2+): %8.3f (KMT)  \n\n',data.final.grand_tot_biomass2y);
    end
    
    toc
    
    if para.proc.kriging == 1
        para.proc.kriging = 0;
        krig_flag = 1;
    else
        krig_flag = 0;
    end
    CVjh = jolly_hampton_CV(para, data);
    CVjh_mean = nanmean(CVjh);
    fprintf('Jolly-Hampton (un-kriged): CV = %6.4f \n\n',CVjh_mean)
    para.proc.kriging = krig_flag;
%     end
    
    %% check whether kriging will be performed
    if para.proc.kriging == 1
        fprintf('\n Estimating kriged biomass using Geostatistics ...\n')
        para.krig.proc_opt=1;
        para.dataprep.fileID='Current Data';
        %%%%% add zeros at the end of west bound transects - number of zeros to be
        %%%%% added is defined in EchoPro_Root_Dir\general\initialization.m file (line 51)
        %          add_zeros([31:32 41:42])
        %          add_zeros([11:12 31:32 41:42])
        %          add_zeros([11:12 17 31:32 41:42 51 53])   % 2013 hake survey - produce the similar biomass as non-kriged biomass estimate.
        kriging_proc
        if isfield(data.final.table,'kriged_biomass0')
            %         data.final.table.biomass=data.final.table.biomass0;   % recover the original
            fprintf('Final Check:\n  Un-kriged biomass = %10.2f (kmt)\n     Kriged biomass = %10.2f (kmt)\n', ...
                nansum(data.final.table.biomass(:,15))*1e-6,nansum(data.final.table.kriged_biomass0(:,10))*1e-6);
        end
        if para.proc.bootstrap.cnt == para.proc.bootstrap.limit & para.proc.bootstrap.cnt > 1
            fprintf(['Final Check:\nMean # Transects = %4d\t std(# Transects) = %4d\n', ...
                    'Mean Kriged biomass (US+CAN)   = %10.2f (kmt)\t   std(Kriged_Biomass) (US+CAN) = %4.2f (kmt)\n', ...
                    'Mean Kriged biomass (CAN)      = %10.2f (kmt)\t   std(Kriged_Biomass) (UCAN) = %4.2f (kmt)\n', ...
                    'Mean Kriged CV (US+CAN)        = %10.2f%% \t std(Kriged CV) (US+CAN)        = %4.2f%%\n', ...
                    'Mean Kriged CV (Jolly-Hampton) = %10.2f%% \t std(Kriged CV) (Jolly-Hampton) = %4.2f%%\n'], ...
                    round(mean(data.final.bootstrap.out(:,2))), round(std(data.final.bootstrap.out(:,2))), ...
                    mean(data.final.bootstrap.out(:,3)), std(data.final.bootstrap.out(:,3)), ...
                    mean(data.final.bootstrap.out(:,4)), std(data.final.bootstrap.out(:,4)), ...
                    mean(data.final.bootstrap.out(:,5))*100, std(data.final.bootstrap.out(:,5))*100, ...
                    mean(data.final.bootstrap.out(:,6))*100, std(data.final.bootstrap.out(:,6))*100);
        end
    else     
         fprintf('Final Check:\n  Un-kriged biomass = %10.2f (kmt)\n',nansum(data.final.table.biomass(:,15))*1e-6);
   end
   fprintf('               NASC = %10.2f x 10^3 (m^2/nmi^2)\n', sum(data.final.table.biomass(:,9)*1e-3))
end

fprintf('End Time = %s\n', datetime)

