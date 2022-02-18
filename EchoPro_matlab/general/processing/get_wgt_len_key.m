function [len_wgt_M, len_wgt_F, len_wgt_ALL]=get_wgt_len_key(stratum_id)
% get length-weight regression keys for Male, Female, and ALL
% i.e. averaged weight per fish (M, F, ALL) based on the length distribution in the specified stratum 
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global data para bio_dat

if para.acoust.TS_station_num == 1
    switch para.proc.len_wgt_method  % using two FSCS' sample stations
        case 1  %% average length-weight-key per fish from whole survey region using length distribution from whole survey region
            if para.proc.survey_len_wgt_key_sex == 1  %% with gender difference
                len_wgt_M=data.bio.strata(stratum_id).Len_key_Mn*data.bio.len_wgt_M(:);
                len_wgt_F=data.bio.strata(stratum_id).Len_key_Fn*data.bio.len_wgt_F(:);
                len_wgt_ALL=data.bio.strata(stratum_id).Len_key_Nn*data.bio.len_wgt_ALL(:);
            else   %% with no gender difference, i.e. use both gender data in len-age regression analysis
                len_wgt_M=data.bio.strata(stratum_id).Len_key_Mn*data.bio.len_wgt_ALL(:);
                len_wgt_F=data.bio.strata(stratum_id).Len_key_Fn*data.bio.len_wgt_ALL(:);
                len_wgt_ALL=data.bio.strata(stratum_id).Len_key_Nn*data.bio.len_wgt_ALL(:);
            end
        case 2  %% average length-weight-key per fish from whole survey region using length distribution from whole survey region with sex difference
                %% this option is the same as para.proc.len_wgt_method=1 & para.proc.survey_len_wgt_key_sex = 1
            len_wgt_M=data.bio.strata(stratum_id).Len_key_Mn*data.bio.len_wgt_M(:);
            len_wgt_F=data.bio.strata(stratum_id).Len_key_Fn*data.bio.len_wgt_F(:);
            len_wgt_ALL=data.bio.strata(stratum_id).Len_key_Nn*data.bio.len_wgt_ALL(:);
        case 3  %% average length-weight-key per fish for every stratum using length distribution within specified stratum
            len_wgt_M= data.bio.strata(stratum_id).Len_key_Mn*data.bio.strata(stratum_id).Len_wgt_key_M(:);
            len_wgt_F= data.bio.strata(stratum_id).Len_key_Fn*data.bio.strata(stratum_id).Len_wgt_key_F(:);
            len_wgt_ALL= data.bio.strata(stratum_id).Len_key_Nn*data.bio.strata(stratum_id).Len_wgt_key_ALL(:);
        case 4  %% average length-weight-key per fish for every stratum using mean length within specified stratum
            LaveM=data.bio.strata(stratum_id).LaveM;   % mean length of hake in assigned stratum corresponding to trawl kk
            len_wgt_M=data.bio.len_wgt_all.reg_w0M*LaveM^data.bio.len_wgt_all.reg_pM;
            LaveF=data.bio.strata(stratum_id).LaveF;   % mean length of hake in assigned stratum corresponding to trawl kk
            len_wgt_F=data.bio.len_wgt_all.reg_w0F*LaveF^data.bio.len_wgt_all.reg_pF;
            len_wgt_ALL=data.bio.len_wgt_all.reg_w0*data.bio.strata(stratum_id).Lave^data.bio.len_wgt_all.reg_p;
    end
else  % using two FSCS' sample stations
    %% proportion of samples at station #2
    fac2_M=data.bio.strata(stratum_id).Len_Age_M_proportion;
    fac2_F=data.bio.strata(stratum_id).Len_Age_F_proportion;
    fac2_ALL=fac2_M+fac2_F;
    %% proportion of samples at station #1
    fac1_ALL=1-fac2_ALL;
    fac1_M=fac1_ALL*sum(data.bio.strata(stratum_id).Len_key_Mn);
    fac1_F=fac1_ALL*sum(data.bio.strata(stratum_id).Len_key_Fn);
    fac1_N=fac1_ALL*sum(data.bio.strata(stratum_id).Len_key_Nn);   % unsexed
    if fac1_M == 0 & fac1_F == 0 & fac2_M == 0 & fac2_F == 0
        % only unsexed data
%         fac1_M=0.5;
%         fac1_F=0.5;
    end
    fac1_M=fac1_M/(fac1_M+fac2_M);
    fac1_F=fac1_F/(fac1_F+fac2_F);
    fac1_ALL=fac1_ALL/(fac1_ALL+fac2_ALL);
    fac2_M=fac2_M/(fac1_M+fac2_M);
    fac2_F=fac2_F/(fac1_F+fac2_F);
    fac2_ALL=fac2_ALL/(fac1_ALL+fac2_ALL);
    switch para.proc.len_wgt_method  %
        case 1  %% average length-weight-key per fish from whole survey region using length distribution from whole survey region
            Len_Age_key_ALL=data.bio.strata(stratum_id).Len_Age_key_ALLn;
            if para.proc.survey_len_wgt_key_sex == 1  %% with gender difference
                len_wgt_M=(fac1_M*data.bio.strata(stratum_id).Len_key_Mn+fac2_M*nansum(data.bio.strata(stratum_id).Len_Age_key_Mn'))*data.bio.len_wgt_M(:);
                len_wgt_F=(fac1_F*data.bio.strata(stratum_id).Len_key_Fn+fac2_F*nansum(data.bio.strata(stratum_id).Len_Age_key_Fn'))*data.bio.len_wgt_F(:);
                len_wgt_ALL=(fac1_ALL*data.bio.strata(stratum_id).Len_key_Nn+fac2_ALL*nansum(Len_Age_key_ALL'))*data.bio.len_wgt_ALL(:);
            else   %% with no gender difference, i.e. use both gender data in len-age regression analysis
                len_wgt_M=(fac1_M*data.bio.strata(stratum_id).Len_key_Mn+fac2_M*nansum(data.bio.strata(stratum_id).Len_Age_key_Mn'))*data.bio.len_wgt_ALL(:);
                len_wgt_F=(fac1_F*data.bio.strata(stratum_id).Len_key_Fn+fac2_F*nansum(data.bio.strata(stratum_id).Len_Age_key_Fn'))*data.bio.len_wgt_ALL(:);
                len_wgt_ALL=(fac1_ALL*data.bio.strata(stratum_id).Len_key_Nn+fac2_ALL*nansum(Len_Age_key_ALL'))*data.bio.len_wgt_ALL(:);
            end
        case 2  %% average length-weight-key per fish from whole survey region using length distribution from whole survey region
                %% this option is the same as para.proc.len_wgt_method=1 & para.proc.survey_len_wgt_key_sex = 1
            len_wgt_M=data.bio.strata(stratum_id).Len_key_Mn*data.bio.len_wgt_M(:);
            len_wgt_F=data.bio.strata(stratum_id).Len_key_Fn*data.bio.len_wgt_F(:);
            len_wgt_ALL=data.bio.strata(stratum_id).Len_key_Nn*data.bio.len_wgt_ALL(:);
        case 3  %% average length-weight-key per fish for every stratum using length distribution within specified stratum
            len_wgt_M= data.bio.strata(stratum_id).Len_key_Mn*data.bio.strata(stratum_id).Len_wgt_key_M(:);
            len_wgt_F= data.bio.strata(stratum_id).Len_key_Fn*data.bio.strata(stratum_id).Len_wgt_key_F(:);
            len_wgt_ALL= data.bio.strata(stratum_id).Len_key_Nn*data.bio.strata(stratum_id).Len_wgt_key_ALL(:);
        case 4  %% average length-weight-key per fish for every stratum using mean length within specified stratum
            LaveM=data.bio.strata(stratum_id).LaveM;   % mean length of hake in assigned stratum corresponding to trawl kk
            len_wgt_M=data.bio.len_wgt_all.reg_w0M*LaveM^data.bio.len_wgt_all.reg_pM;
            LaveF=data.bio.strata(stratum_id).LaveF;   % mean length of hake in assigned stratum corresponding to trawl kk
            len_wgt_F=data.bio.len_wgt_all.reg_w0F*LaveF^data.bio.len_wgt_all.reg_pF;
            len_wgt_ALL=data.bio.len_wgt_all.reg_w0*data.bio.strata(stratum_id).Lave^data.bio.len_wgt_all.reg_p;
    end
end
% if length(len_wgt_M) > 1 | length(len_wgt_F) > 1 | length(len_wgt_ALL) > 1
%     disp([length(len_wgt_M) length(len_wgt_F) length(len_wgt_ALL)])
%     disp(' =======================================')
% end
return