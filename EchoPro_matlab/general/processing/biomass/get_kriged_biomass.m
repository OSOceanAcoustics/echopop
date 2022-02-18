function [var,var1,cv_s,C0,data]=get_kriged_biomass(data,para)
%% function get_kriged_biomass.m is to obtain the biomass density, variance and
%% acoustically weight length-age-sex strctured abundance & biomass tables
%% based on the kriged input variables (para.proc.kriging_input):
%% 1 = biomass density
%% 2 = NASC
%% 3 = number density
%% OUTPUTS:
%% var = kriged value;
%% var2 =kriging varaince
%% cv_s = local sample CV within the kriging search redius
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       12/12/2013


%% initialize acoustically weighted length-age-sex structured output matrices (over write the non-kriged tables)
Final_out.Len_Age_Matrix_AcoustM=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for male hake
Final_out.Len_Age_Matrix_AcoustF=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for female hake
Final_out.Len_Age_Matrix_AcoustALL=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for all hake

%% initialize weight-key as a function of length and age
Final_Wgt_len_age_M=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
Final_Wgt_len_age_F=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
Final_Wgt_len_age_ALL=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));

Final_Wgt_len_M=zeros(length(para.bio.hake_len_bin),1);
Final_Wgt_len_F=zeros(length(para.bio.hake_len_bin),1);
Final_Wgt_len_ALL=zeros(length(para.bio.hake_len_bin),1);


%% Kriged NASC or Number density
haul_stratum_id=data.in.US_CAN_Mesh_stratum;
% ind = find(haul_stratum_id > length(data.bio.strata));
% if ~isempty(ind)
%     haul_stratum_id(ind) = length(data.bio.strata);
% end

strata_lat=unique(haul_stratum_id);

% d=xlsread(data.in.filename.grid_cell);
A0=para.proc.kriging_A0;                % nominal area of each grid
% Area=A0*d(:,5);            % actual cell area in nmi^2
Area=data.in.US_CAN_Mesh_area;
Nc=length(Area);

%% 'stratum','NASC','ntk_male','ntk_female','ntk_total','wgt_male','wgt_female','wgt_total'
stratum=nan*ones(Nc,1);ntk=nan*ones(Nc,1);ntk_male=nan*ones(Nc,1);ntk_female=nan*ones(Nc,1);
NASC=nan*ones(Nc,1);wgt_int_M=nan*ones(Nc,1);wgt_int_F=nan*ones(Nc,1);wgt_int_ALL=nan*ones(Nc,1);
n_len_bin=length(para.bio.hake_len_bin);
n_age_bin=length(para.bio.hake_age_bin);
if para.proc.kriging_input == 1
    %% Kriged Biomass density
   data.out.krig.gv(data.out.krig.gv<0 | isnan(data.out.krig.gv) == 1)=0;
    var=data.out.krig.gv;
    var1=data.out.krig.ge;
    C0=data.out.vario.c0;
    cv_s=data.out.krig.ge_s;            % sample variance
    %% create weight-based length, age, and sex structured 2D keys (matrices)
    %% from abundance-based length, age, and sex structured 2D keys (matrices)
    
    %% ------ separation of aged and unaged quantities ----
    % create weight-based len-age key from number- or abundance-based length-age keys
    for i=1:size(para.bio.strata,2)   % loop through number of strata
        if ~isempty(para.bio.strata(i).trawls)
            %% aged and unaged len-age distributed biomass proportion within each stratum
            aged_proportion(i)=para.bio.strata(i).Len_Age_M_wgt_proportion+para.bio.strata(i).Len_Age_F_wgt_proportion;
            unaged_proportion(i)=1-aged_proportion(i);
            %% aged male and female weight proportion within each stratum
            Dist_Wgt(i).M=para.bio.strata(i).Len_Age_M_wgt_proportion*para.bio.strata(i).Len_Age_key_wgt_Mn;
            Dist_Wgt(i).F=para.bio.strata(i).Len_Age_F_wgt_proportion*para.bio.strata(i).Len_Age_key_wgt_Fn;
            while isnan(para.bio.strata(i).Len_Age_key_wgt_Mn(40,20))
                % this is caused by only sampling at sation 1 but no samples at station 2
                % find a stratum has similar mean length distribution to replace the Len_Age_key_wgt table
                Dist_Wgt(i).M=zeros(n_len_bin,n_age_bin);
                Dist_Wgt(i).F=zeros(n_len_bin,n_age_bin);
                para.bio.strata(i).Len_Age_key_wgt_Mn(40,20)=0;
            end
            %% weight per unit length distribution for unaged hake
            W_Ln_ALL_array=data.bio.len_wgt_ALL.*para.bio.strata(i).Len_key_Nn;
            W_Ln_ALL_array_sum(i)=sum(W_Ln_ALL_array);
            % average weight in kg per fish
            if W_Ln_ALL_array_sum(i) ~= 0  % there are fish sampled in biological sample station #1
                W_Ln_ALL(i).N=W_Ln_ALL_array(:)/W_Ln_ALL_array_sum(i);    % normalized weight per unit length distribution
                if (para.bio.strata(i).Len_M_wgt_proportion+para.bio.strata(i).Len_F_wgt_proportion) == 0
                    %% if the station #1 are all zeros, using 50% for male and female ratio
                    M_proportion=0.5;
                    F_proportion=0.5;
                else
                    M_proportion=para.bio.strata(i).Len_M_wgt_proportion/(para.bio.strata(i).Len_M_wgt_proportion+para.bio.strata(i).Len_F_wgt_proportion);
                    F_proportion=para.bio.strata(i).Len_F_wgt_proportion/(para.bio.strata(i).Len_M_wgt_proportion+para.bio.strata(i).Len_F_wgt_proportion);
                end
                unaged_M_wgt_proportion(i)=unaged_proportion(i)*M_proportion;
                unaged_F_wgt_proportion(i)=unaged_proportion(i)*F_proportion;
            else  % there are no fish sampled in biological sample station #1
                W_Ln_ALL(i).N=0;
                unaged_M_wgt_proportion(i)=0;
                unaged_F_wgt_proportion(i)=0;
            end
        end
    end
%     data.chk.ind_stratum_age12_age13=[];
    for ii=1:length(var)        % loop through grid cells
%         if ii > 480 & ii < 483
%             disp(ii)
%         end
        if ~isnan(var(ii))
            strata_indx=data.in.US_CAN_Mesh_stratum(ii);
            if strata_indx == 0
                strata_indx=1;
                fprintf('strat index =0: ii = %d\n',ii);
            end
            if data.bio.strata(haul_stratum_id(ii)).LaveALL < para.bio.Lmin  % remove YOY hake
                var(ii)=0;
            end
            if haul_stratum_id(ii) ~= 0 & ~isnan(haul_stratum_id(ii))
                sig_b(ii)=data.bio.strata(strata_indx).sig_b;
            else
                sig_b(ii)=1e10;
            end

            %% aged and non-aged biomass and abundance
            %% aged biomass
            Wgt_len_age_M_ii=var(ii)*Dist_Wgt(strata_indx).M*Area(ii);
            Wgt_len_age_F_ii=var(ii)*Dist_Wgt(strata_indx).F*Area(ii);
            Final_Wgt_len_age_M=Final_Wgt_len_age_M+Wgt_len_age_M_ii;
            Final_Wgt_len_age_F=Final_Wgt_len_age_F+Wgt_len_age_F_ii;
            Wgt_len_age_ALL_ii=var(ii)*aged_proportion(strata_indx)*para.bio.strata(strata_indx).Len_Age_key_wgt_ALLn*Area(ii);
            Final_Wgt_len_age_ALL=Final_Wgt_len_age_ALL+Wgt_len_age_ALL_ii;
            %% unaged biomass
            %% increament biomass within the current grid cell
            Wgt_len_M_ii=var(ii)*unaged_M_wgt_proportion(strata_indx)*Area(ii)*W_Ln_ALL(strata_indx).N;
            Wgt_len_F_ii=var(ii)*unaged_F_wgt_proportion(strata_indx)*Area(ii)*W_Ln_ALL(strata_indx).N;
            Wgt_len_ALL_ii=var(ii)*unaged_proportion(strata_indx)*Area(ii)*W_Ln_ALL(strata_indx).N;
            %% accumulative biomass 2-D matrix
            Final_Wgt_len_M=Final_Wgt_len_M+Wgt_len_M_ii;
            Final_Wgt_len_F=Final_Wgt_len_F+Wgt_len_F_ii;
            Final_Wgt_len_ALL=Final_Wgt_len_ALL+var(ii)*unaged_proportion(strata_indx)*Area(ii)*W_Ln_ALL(strata_indx).N;
            if any(isnan( Final_Wgt_len_age_M)) & var(ii)~= 0
                disp(ii)
            end
            %% 1D table
            %% aged and un-aged biomass within the grid cell
            wgt_int_M(ii)=nansum(nansum(Wgt_len_age_M_ii,2)+Wgt_len_M_ii);
            wgt_int_F(ii)=nansum(nansum(Wgt_len_age_F_ii,2)+Wgt_len_F_ii);
            wgt_int_ALL(ii)=nansum(nansum(Wgt_len_age_ALL_ii,2)+Wgt_len_ALL_ii);
%             if ii == 1913
%                 disp(ii)
%             end
            %% aged and un-aged abundance & NASC within the grid cell
            [len_wgt_M, len_wgt_F, len_wgt_ALL]=get_wgt_len_key(haul_stratum_id(ii));
            %% use all fish regardless of sex and convert weight to numbers
            ntk_male(ii)=wgt_int_M(ii)/len_wgt_ALL;
            ntk_female(ii)=wgt_int_F(ii)/len_wgt_ALL;
%             ntk(ii)=ntk_male(ii)+ntk_female(ii);
            ntk(ii)=wgt_int_ALL(ii)/len_wgt_ALL;                % modified on 12/30/2015 to include un-sexed biomass
            NASC(ii)=ntk(ii)*data.bio.strata(strata_indx).sig_b;
            stratum(ii)=strata_indx;
            if abs(ntk(ii)) == 0 &  wgt_int_ALL(ii) > eps
                fprintf('inconsistency in abundance and biomass: %d\t %d\t %10.3f\n',ii,ntk(ii),wgt_int_ALL(ii));
            end
            
        else
            NASC(ii)=0;
            wgt_int_M(ii)=0;                                % in kg
            wgt_int_F(ii)=0;                                % in kg
            ntk(ii)=0;  % in kg for unsexed hake
            wgt_int_ALL(ii)=0;
            ntk_male(ii)=0;
            ntk_female(ii)=0;          
            stratum(ii)=haul_stratum_id(ii);
        end
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of 2/9/2013 revision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isnan(ntk(ii))
            %% proportion the numbers+
            %% Length station (#1) -unaged abundance
            N_len_M=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_M_proportion;
            N_len_F=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_F_proportion;
            %% length and age station (#2) -aged abundance
            N_len_age_M=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion;
            N_len_age_F=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion;
            %% update Len-Age matrix for unaged hake
            Final_out.Len_Age_Matrix_AcoustM(:,end)=Final_out.Len_Age_Matrix_AcoustM(:,end)+(N_len_M*para.bio.strata(haul_stratum_id(ii)).Len_key_Mn(:));
            Final_out.Len_Age_Matrix_AcoustF(:,end)=Final_out.Len_Age_Matrix_AcoustF(:,end)+(N_len_F*para.bio.strata(haul_stratum_id(ii)).Len_key_Fn(:));
            Final_out.Len_Age_Matrix_AcoustALL(:,end)=Final_out.Len_Age_Matrix_AcoustALL(:,end) ...
                + (N_len_M*para.bio.strata(haul_stratum_id(ii)).Len_key_Mn(:) + N_len_F*para.bio.strata(haul_stratum_id(ii)).Len_key_Fn(:));
            %% update Len-Age matrix for aged hake
            if N_len_age_M ~=0
                Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)+(N_len_age_M*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Mn(:,1:end));
                Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1) ...
                    +(N_len_age_M*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Mn(:,1:end));
            end
            if N_len_age_F ~=0
                Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)+(N_len_age_F*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Fn(:,1:end));
                Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1) ...
                    +(N_len_age_F*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Fn(:,1:end));
            end
            coef(ii)=para.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion+para.bio.strata(haul_stratum_id(ii)).Len_F_proportion ...
                   + para.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion +para.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion;
            %% mean weight per hake in the specified stratum
%             if isnan(sum(sum(Final_out.Len_Age_Matrix_AcoustM)))
%                 disp(ii)
%             end
%             if sum(Final_out.Len_Age_Matrix_AcoustALL(:,13)) > sum(Final_out.Len_Age_Matrix_AcoustALL(:,12))
%                 disp([ii haul_stratum_id(ii) sum(Final_out.Len_Age_Matrix_AcoustALL(:,12:13))])
%                 data.chk.ind_stratum_age12_age13=[data.chk.ind_stratum_age12_age13; [ii haul_stratum_id(ii) sum(Final_out.Len_Age_Matrix_AcoustALL(:,12:13))]];
%             end
        end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of 2/9/2013 revision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
       

    end  % end grid cell loop
%% combine aged and unaged fish to form a 2D Matrix  %%%%%
    Final_Wgt_len_age_M=[Final_Wgt_len_age_M Final_Wgt_len_M]*1e-6;         % kg -> 1000 tons
    Final_Wgt_len_age_F=[Final_Wgt_len_age_F Final_Wgt_len_F]*1e-6;         % kg -> 1000 tons
    Final_Wgt_len_age_ALL=[Final_Wgt_len_age_ALL Final_Wgt_len_ALL]*1e-6;   % kg -> 1000 tons

    %% create acoustically weighted kriged matrix
    nA=length(para.bio.hake_age_bin);
    nL=length(para.bio.hake_len_bin);
    
    %% distribution of unaged biomass to aged biomass
    Wgt_len_age_M_norm=zeros(nL,nA);
    Wgt_len_age_F_norm=zeros(nL,nA);
    Uaged2aged_mat_M=Wgt_len_age_M_norm;
    Uaged2aged_mat_F=Wgt_len_age_F_norm;
    threshold=1e-10;
    for  i=1:nL
        if sum(Final_Wgt_len_age_M(i,1:nA)) < threshold & sum(Final_Wgt_len_age_F(i,1:nA)) < threshold 
            %% neither male and nor female hake at ith length bin is found for aged hake (bio-sample station 2) 
            if Final_Wgt_len_age_M(i,end) < threshold  
            %% no unaged hake found at bio-sample station 1
                Uaged2aged_mat_M(i,:)=0;
                Uaged2aged_mat_F(i,:)=0;            
            else 
            %% unaged hake found at bio-sample station 1
                sum_over_ageM=sum(Final_Wgt_len_age_M(:,1:nA),2);
                ind_M=find(sum_over_ageM > eps);
                [val,ind_sel_M]=min(abs(ind_M-i));
                ind_sel_M=ind_sel_M+min(ind_M)-1;
                sum_over_ageF=sum(Final_Wgt_len_age_F(:,1:nA),2);
                ind_F=find(sum_over_ageF > eps);
                [val,ind_sel_F]=min(abs(ind_F-i));
                ind_sel_F=ind_sel_F+min(ind_F)-1;
                if ind_sel_F < ind_sel_M
                    %% closet length bin has no aged female, using the smaller length bin aged female data 
                    while sum(Final_Wgt_len_age_F(ind_sel_F,1:nA)) == 0
                        ind_sel_F=ind_sel_F-1;
                    end
                    Uaged2aged_mat_M(i,:)=Final_Wgt_len_age_M(i,end)*Final_Wgt_len_age_F(ind_sel_F,1:nA)/sum(Final_Wgt_len_age_F(ind_sel_F,1:nA));
                    Uaged2aged_mat_F(i,:)=Final_Wgt_len_age_F(i,end)*Final_Wgt_len_age_F(ind_sel_F,1:nA)/sum(Final_Wgt_len_age_F(ind_sel_F,1:nA));
                else
                    %% closet length bin has no aged male, using the smaller length bin aged male data 
                    while sum(Final_Wgt_len_age_M(ind_sel_M,1:nA)) == 0
                        ind_sel_M=ind_sel_M-1;
                    end
                    Uaged2aged_mat_M(i,:)=Final_Wgt_len_age_M(i,end)*Final_Wgt_len_age_M(ind_sel_M,1:nA)/sum(Final_Wgt_len_age_M(ind_sel_M,1:nA));
                    Uaged2aged_mat_F(i,:)=Final_Wgt_len_age_F(i,end)*Final_Wgt_len_age_M(ind_sel_M,1:nA)/sum(Final_Wgt_len_age_M(ind_sel_M,1:nA));
                end
            end
        elseif sum(Final_Wgt_len_age_M(i,1:nA)) < threshold & Final_Wgt_len_age_M(i,end) > threshold
            %% no male hake at ith length bin for aged hake (bio-sample station 2) but has for unaged hake (bio-sample station 1)
            Uaged2aged_mat_M(i,:)=Final_Wgt_len_age_M(i,end)*Final_Wgt_len_age_F(i,1:nA)/sum(Final_Wgt_len_age_F(i,1:nA));
            Uaged2aged_mat_F(i,:)=Final_Wgt_len_age_F(i,end)*Final_Wgt_len_age_F(i,1:nA)/sum(Final_Wgt_len_age_F(i,1:nA));
        elseif sum(Final_Wgt_len_age_F(i,1:nA)) < threshold & Final_Wgt_len_age_F(i,end) > threshold
            %% no female hake at ith length bin for aged hake (bio-sample station 2) but has for unaged hake (bio-sample station 1)
            Uaged2aged_mat_M(i,:)=Final_Wgt_len_age_M(i,end)*Final_Wgt_len_age_M(i,1:nA)/sum(Final_Wgt_len_age_M(i,1:nA));
            Uaged2aged_mat_F(i,:)=Final_Wgt_len_age_F(i,end)*Final_Wgt_len_age_M(i,1:nA)/sum(Final_Wgt_len_age_M(i,1:nA));
        elseif sum(Final_Wgt_len_age_M(i,1:nA)) > threshold & sum(Final_Wgt_len_age_F(i,1:nA)) > threshold & Final_Wgt_len_age_ALL(i,end) > threshold
            %% both male and female hake have samples at ith length bin for aged hake (bio-sample station 2) and unaged hake (bio-sample station 1)
            Uaged2aged_mat_M(i,:)=Final_Wgt_len_age_M(i,end)*Final_Wgt_len_age_M(i,1:nA)/sum(Final_Wgt_len_age_M(i,1:nA));
            Uaged2aged_mat_F(i,:)=Final_Wgt_len_age_F(i,end)*Final_Wgt_len_age_F(i,1:nA)/sum(Final_Wgt_len_age_F(i,1:nA));
        end
        if any(isnan(Uaged2aged_mat_M))
            disp(i)
        end
    end
    
    %% combined aged and age-distributed from unaged biomass
    Final_Wgt_len_age_M0=Final_Wgt_len_age_M(:,1:length(para.bio.hake_age_bin))+Uaged2aged_mat_M;
    Final_Wgt_len_age_F0=Final_Wgt_len_age_F(:,1:length(para.bio.hake_age_bin))+Uaged2aged_mat_F;
    Final_Wgt_len_age_ALL0=Final_Wgt_len_age_M0+Final_Wgt_len_age_F0;
    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM=Final_Wgt_len_age_M0;
    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF=Final_Wgt_len_age_F0;
    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL=Final_Wgt_len_age_ALL0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start of 2/9/2013 revision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     %% abundance tables
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustM=Final_out.Len_Age_Matrix_AcoustM;
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustF=Final_out.Len_Age_Matrix_AcoustF;
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL=Final_out.Len_Age_Matrix_AcoustALL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of 2/9/2013 revision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% construct 1D kriged NASC, abundance, and biomass table (grid cell)
% data.final.table.kriged_biomass0_description
%   Columns 1 through 4
%     'Lat'    'Lon'    'stratum'    'NASC'
%   Columns 5 through 10
%     'ntk_male'    'ntk_female'    'ntk_total'    'wgt_male'    'wgt_female'    'wgt_total'  
%   Columns 11 through 12
%     'sig_b'    'krig_CV'

    data.final.table.kriged_biomass0=[data.out.krig.lat data.out.krig.lon haul_stratum_id(:) NASC(:) ntk_male(:) ntk_female(:) ntk(:) wgt_int_M(:) wgt_int_F(:) wgt_int_ALL(:) sig_b(:)];
    data.final.table.kriged_biomass0_description={'Lat','Lon','stratum','NASC','ntk_male','ntk_female','ntk_total','wgt_male','wgt_female','wgt_total','sig_b','krig_CV','krig_SD'};
else
    %% kried NASC or number density  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if para.proc.kriging_input == 2
        %% kriged NASC
        data.out.krig.gv(data.out.krig.gv<0 | isnan(data.out.krig.gv) == 1)=0;
        NASC=data.out.krig.gv;
        nntk=zeros(size(NASC));
        for i=1:size(strata_lat,1)
            indx=find(haul_stratum_id == i);
            nntk(indx)=NASC(indx)/para.bio.strata(i).sig_b;   % number density
        end
    else
        %% Kriged fish number density
        nntk=data.out.krig.gv;
    end
    ntk=nntk.*data.in.US_CAN_Mesh_area;           % total # hake within mesh cells
    dS=0;
    for ii=1:length(data.out.krig.gv)             % loop through grid cells
        if data.bio.strata(haul_stratum_id(ii)).LaveALL < para.bio.Lmin  % remove YOY hake
            ntk(ii)=0;
        end
        NASC(ii)=ntk(ii)*data.bio.strata(haul_stratum_id(ii)).sig_b;
        if haul_stratum_id(ii) ~= 0 & ~isnan(haul_stratum_id(ii))
            sig_b(ii)=data.bio.strata(haul_stratum_id(ii)).sig_b;
        else
            sig_b(ii)=1e10;
        end
        if ~isnan(ntk(ii))
            %% proportion the numbers+
            %% Length station (#1)
            N_len_M=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_M_proportion;
            N_len_F=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_F_proportion;
            %% length and age station (#2)
            N_len_age_M=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion;
            N_len_age_F=ntk(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion;
            %% update Len-Age matrix for unaged hake
            Final_out.Len_Age_Matrix_AcoustM(:,end)=Final_out.Len_Age_Matrix_AcoustM(:,end)+(N_len_M*para.bio.strata(haul_stratum_id(ii)).Len_key_Mn(:));
            Final_out.Len_Age_Matrix_AcoustF(:,end)=Final_out.Len_Age_Matrix_AcoustF(:,end)+(N_len_F*para.bio.strata(haul_stratum_id(ii)).Len_key_Fn(:));
            Final_out.Len_Age_Matrix_AcoustALL(:,end)=Final_out.Len_Age_Matrix_AcoustALL(:,end) ...
                + (N_len_M*para.bio.strata(haul_stratum_id(ii)).Len_key_Mn(:) + N_len_F*para.bio.strata(haul_stratum_id(ii)).Len_key_Fn(:));
            %% update Len-Age matrix for aged hake
            if N_len_age_M ~=0
                Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)+(N_len_age_M*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Mn(:,1:end));
                Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1) ...
                    +(N_len_age_M*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Mn(:,1:end));
            end
            if N_len_age_F ~=0
                Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)+(N_len_age_F*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Fn(:,1:end));
                Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1) ...
                    +(N_len_age_F*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_Fn(:,1:end));
            end
            %% mean weight per hake in the specified stratum
            [len_wgt_M, len_wgt_F, len_wgt_ALL]=get_wgt_len_key(haul_stratum_id(ii));
            %%  sum of the proportions of station1 and station 2
            ntk_male(ii)=(ntk(ii)*(para.bio.strata(haul_stratum_id(ii)).Len_M_proportion+para.bio.strata(haul_stratum_id(ii)).Len_Age_M_proportion));
            ntk_female(ii)=(ntk(ii)*(para.bio.strata(haul_stratum_id(ii)).Len_F_proportion+para.bio.strata(haul_stratum_id(ii)).Len_Age_F_proportion));

            %% convert average number of fish to averaged weight in kgch
            Wgt_male_int(ii)=ntk_male(ii)*len_wgt_ALL;                                % in kg
            Wgt_female_int(ii)=ntk_female(ii)*len_wgt_ALL;                            % in kg
            Wgt_unsexed_int(ii)=(ntk(ii)-ntk_male(ii)-ntk_female(ii))*len_wgt_ALL;  % in kg for unsexed hake
            Wgt_ALL_int(ii)=Wgt_male_int(ii)+Wgt_female_int(ii)+Wgt_unsexed_int(ii);

            %% acumulative biomass
            Final_Wgt_len_age_M=Final_Wgt_len_age_M+Wgt_male_int(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
            Final_Wgt_len_age_F=Final_Wgt_len_age_F+Wgt_female_int(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;
            Final_Wgt_len_age_ALL=Final_Wgt_len_age_ALL+Wgt_ALL_int(ii)*para.bio.strata(haul_stratum_id(ii)).Len_Age_key_wgt_ALLn;            
        else  % NASC(ii) or ntk(ii) is NaN
            NASC(ii)=0;
            Wgt_male_int(ii)=0;                                % in kg
            Wgt_female_int(ii)=0;                              % in kg
            Wgt_unsexed_int(ii)=0;  % in kg for unsexed hake
            Wgt_ALL_int(ii)=0;
            ntk_male(ii)=0;
            ntk_female(ii)=0;
            ntk(ii)=0;
            stratum(ii)=haul_stratum_id(ii);
        end  % end of non-NaN loop
    end   % end grid cell loop
    
    %% total biomass
    data.final.grand_male_biomass=nansum(Wgt_male_int)*1e-6;            % kg --> kmt or 1000 tons
    data.final.grand_female_biomass=nansum(Wgt_female_int)*1e-6;        % kg --> kmt or 1000 tons
    data.final.grand_tot_biomass=nansum(Wgt_ALL_int)*1e-6;              % kg --> kmt or 1000 tons
    
    %% create final length-age-sex structured tables
    %% abundance tables
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustM=Final_out.Len_Age_Matrix_AcoustM;
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustF=Final_out.Len_Age_Matrix_AcoustF;
    data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL=Final_out.Len_Age_Matrix_AcoustALL;

    %% kriged biomass tables
    %% unsexed
    unsex_total_biomass=data.final.grand_tot_biomass2y-(data.final.grand_male_biomass+data.final.grand_female_biomass);
    data.final.grand.unsexed_biomass=nansum(Wgt_unsexed_int)*1e-6;

    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM=Final_Wgt_len_age_M*1e-3;     % kmt --> mmt
    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF=Final_Wgt_len_age_F*1e-3;     % kmt --> mmt
    data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL=Final_Wgt_len_age_ALL*1e-3; % kmt --> mmt

    var=Wgt_ALL_int(:)./data.in.US_CAN_Mesh_area;   % weight or biomass density
    var1=data.out.krig.ge;                          % this is a normalized and dimensionless quantity
    C0=data.out.vario.c0;
    cv_s=data.out.krig.ge_s;                        % sample varaince
    %% construct 1D kriged NASC, abundance, and biomass table (grid cell)
% data.final.table.kriged_biomass0_description
%   Columns 1 through 4
%     'Lat'    'Lon'    'stratum'    'NASC'
%   Columns 5 through 10
%     'ntk_male'    'ntk_female'    'ntk_total'    'wgt_male'    'wgt_female'    'wgt_total'  
%   Columns 11 through 12
%     'sig_b'    'krig_CV'

    data.final.table.kriged_biomass0=[data.out.krig.lat data.out.krig.lon haul_stratum_id(:) NASC(:) ntk_male(:) ntk_female(:) ntk(:) Wgt_male_int(:) Wgt_female_int(:) Wgt_ALL_int(:) sig_b(:)];
    data.final.table.kriged_biomass0_description={'Lat','Lon','stratum','NASC','ntk_male','ntk_female','ntk_total','wgt_male','wgt_female','wgt_total','sig_b','krig_CV' };
end
%% add length and age labels in the 2D arrays
%% biomass
if para.proc.exclude_age1 == 1
    data=modify_abundance_biomass_matrices(data,para);
end
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM*1e-3];   % kmt --> mmt
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF*1e-3];
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL=[[0 para.bio.hake_age_bin]; para.bio.hake_len_bin(:) data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL*1e-3];
%% abundance
data.final.table.kriged_Num_Len_Age_Matrix_AcoustM=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) round(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM)];
data.final.table.kriged_Num_Len_Age_Matrix_AcoustF=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) round(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF)];
data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL=[[0 para.bio.hake_age_bin nan]; para.bio.hake_len_bin(:) round(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL)];

return

