function data=modify_abundance_biomass_matrices(data,para);
%% if age-1 is excluded, the abundance and biomass matries based on trawls need to 
%% be revised. The leaked abundance and biomass to age-1 should be distributed to 
%% quantities of other ages
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/2/2013


if size(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM,2) > length(para.bio.hake_age_bin)+1 % 20 age-bins + unaged column
    %% with length as the first column
    ind0=2;
    row_offset=1;
else
    %% without length as the first column
    ind0=1;
    row_offset=0;
end
%% abundance at the i'th length bin
%% male
num_age1_sum=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(1+row_offset:end,ind0));
num_ind=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end-1);
num_age_prop=num_age1_sum*num_ind/sum(nansum(num_ind));
data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(row_offset+1+1:end,ind0)=0;
rev_num=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end-1)+num_age_prop;
data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end-1)=rev_num;
%% female
num_age1_sum=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0));
num_ind=data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end-1);
num_age_prop=num_age1_sum*num_ind/sum(nansum(num_ind));
data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0)=0;
rev_num=data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end-1)+num_age_prop;
data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end-1)=rev_num;
%% ALL
num_age1_sum=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0));
num_ind=data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end-1);
num_age_prop=num_age1_sum*num_ind/sum(nansum(num_ind));
data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0)=0;
rev_num=data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end-1)+num_age_prop;
data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end-1)=rev_num;

%% biomass at the i'th length bin
%% male
wgt_age1_sum=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0));
wgt_ind=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end);
wgt_age_prop=wgt_age1_sum*wgt_ind/sum(nansum(wgt_ind));
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0)=0;
rev_wgt=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end)+wgt_age_prop;
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(row_offset+1:end,ind0+1:end)=rev_wgt;
%% female
wgt_age1_sum=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0));
wgt_ind=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end);
wgt_age_prop=wgt_age1_sum*wgt_ind/sum(nansum(wgt_ind));
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0)=0;
rev_wgt=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end)+wgt_age_prop;
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(row_offset+1:end,ind0+1:end)=rev_wgt;
%% ALL
wgt_age1_sum=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0));
wgt_ind=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end);
wgt_age_prop=wgt_age1_sum*wgt_ind/sum(nansum(wgt_ind));
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0)=0;
rev_wgt=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end)+wgt_age_prop;
data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(row_offset+1:end,ind0+1:end)=rev_wgt;

return