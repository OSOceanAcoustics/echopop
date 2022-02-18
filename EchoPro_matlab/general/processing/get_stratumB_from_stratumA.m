function [L1, Gender1, L2, Age2, Wgt2, Gender2, data] = get_stratumB_from_stratumA(data, iA, iB)
% construct some matrices for straum iA from Stratum iB
% Dezhang Chu, NOAA Fisheries, NWFSC
% 10-27-2018
% 03-21-2020,  copy all length/age/wgt/sex of station 2, rather than only part of them

% data.bio.strata(iA) = data.bio.strata(iB);
% L1 = data.bio.strata(iB).length1;
% Gender1 = data.bio.strata(iB).gender1;
% L2 = data.bio.strata(iB).length2;
% Age2 = data.bio.strata(iB).age2;
% Wgt2 = data.bio.strata(iB).weight2;
% Gender2 = data.bio.strata(iB).gender2;



%% only part of the samples in station 2
% Lave_iA = data.bio.strata(iA).Lave_j;
% ind = find(data.bio.strata(iB).length2 >= 0.8*Lave_iA & data.bio.strata(iB).length2 <= 1.2*Lave_iA);
ind = 1:length(data.bio.strata(iB).length2);
L1 = data.bio.strata(iA).length1;
Gender1 = data.bio.strata(iA).gender1;
L2 = data.bio.strata(iB).length2(ind);
Age2 = data.bio.strata(iB).age2(ind);
Wgt2 = data.bio.strata(iB).weight2(ind);
Gender2 = data.bio.strata(iB).gender2(ind);
data.bio.strata(iA).length2 = data.bio.strata(iB).length2(ind);
data.bio.strata(iA).age2 = data.bio.strata(iB).age2(ind);
data.bio.strata(iA).weight2 = data.bio.strata(iB).weight2(ind);
data.bio.strata(iA).gender2 = data.bio.strata(iB).gender2(ind);
data.bio.strata(iA).Len_Age_key_M = data.bio.strata(iB).Len_Age_key_M;
data.bio.strata(iA).Len_Age_key_wgt_M = data.bio.strata(iB).Len_Age_key_wgt_M;
data.bio.strata(iA).Len_Age_key_F = data.bio.strata(iB).Len_Age_key_F;
data.bio.strata(iA).Len_Age_key_wgt_F = data.bio.strata(iB).Len_Age_key_wgt_F;
data.bio.strata(iA).Len_Age_key_ALL = data.bio.strata(iB).Len_Age_key_ALL;
data.bio.strata(iA).Len_Age_key_wgt_ALL = data.bio.strata(iB).Len_Age_key_wgt_ALL;
data.bio.strata(iA).Len_Age_key_Mn = data.bio.strata(iB).Len_Age_key_Mn;
data.bio.strata(iA).Len_Age_key_Fn = data.bio.strata(iB).Len_Age_key_Fn;
data.bio.strata(iA).Len_Age_key_ALLn = data.bio.strata(iB).Len_Age_key_ALLn;
data.bio.strata(iA).matrix_NtotM_wgt = data.bio.strata(iB).matrix_NtotM_wgt;
data.bio.strata(iA).matrix_NtotF_wgt = data.bio.strata(iB).matrix_NtotF_wgt;
data.bio.strata(iA).matrix_NtotALL_wgt = data.bio.strata(iB).matrix_NtotALL_wgt;
data.bio.strata(iA).Len_Age_key_wgt_Mn = data.bio.strata(iB).Len_Age_key_wgt_Mn;
data.bio.strata(iA).Len_Age_key_wgt_Fn = data.bio.strata(iB).Len_Age_key_wgt_Fn;
data.bio.strata(iA).Len_Age_key_wgt_ALLn = data.bio.strata(iB).Len_Age_key_wgt_ALLn;
data.bio.strata(iA).Len_Age_M_proportion = data.bio.strata(iB).Len_Age_M_proportion;
data.bio.strata(iA).Len_Age_F_proportion = data.bio.strata(iB).Len_Age_F_proportion;
data.bio.strata(iA).Len_M_proportion = data.bio.strata(iB).Len_M_proportion;
data.bio.strata(iA).Len_F_proportion = data.bio.strata(iB).Len_F_proportion;
data.bio.strata(iA).Len_Age_M_wgt_proportion = data.bio.strata(iB).Len_Age_M_wgt_proportion;
data.bio.strata(iA).Len_Age_F_wgt_proportion = data.bio.strata(iB).Len_Age_F_wgt_proportion;
data.bio.strata(iA).Len_M_wgt_proportion = data.bio.strata(iB).Len_M_wgt_proportion;
data.bio.strata(iA).Len_F_wgt_proportion = data.bio.strata(iB).Len_F_wgt_proportion;
data.bio.strata(iA).Len_wgt_key_M = data.bio.strata(iB).Len_wgt_key_M;
data.bio.strata(iA).Len_wgt_key_F = data.bio.strata(iB).Len_wgt_key_F;
data.bio.strata(iA).Len_wgt_key_F_ind = data.bio.strata(iB).Len_wgt_key_F_ind;
data.bio.strata(iA).Len_wgt_key_ALL = data.bio.strata(iB).Len_wgt_key_ALL;
data.bio.strata(iA).Len_wgt_key_ALL_ind = data.bio.strata(iB).Len_wgt_key_ALL_ind;

end