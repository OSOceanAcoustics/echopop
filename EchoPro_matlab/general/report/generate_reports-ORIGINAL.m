function generate_reports(hdl)
% generate lenght-age-sex structured reports (excel spreadsheet files, .xlsx files)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global data para

tic
%% write acoustically weighted un-kriged lenght-age-gender abundance tables 
files=dir([para.proc.output_filepath '\*.xlsx']);
for i=1:length(files)
    delete([para.proc.output_filepath '\' files(i).name])
end
if get(hdl.radio_assessment,'value') == 1
    disp('write acoustically weighted un-kriged lenght-age-gender abundance tables ...')
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Age (Male)''},1,''L1:L1'');'];
    cmd1txtAx=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Un-aged''},1,''V2:V2'');'];
    %   cmd1txtLsum=['xlswrite(''' para.proc.output_filepath '\len_age_abundance_table.xlsx'',{''Age (Male)''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',data.final.table.Len_Age_Matrix_AcoustM,1,''A2:V42'');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Age (Female)''},2,''L1:L1'');'];
    cmd2txtAx=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Un-aged''},2,''V2:V2'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',data.final.table.Len_Age_Matrix_AcoustF,2,''A2:V42'');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Age (All)''},3,''L1:L1'');'];
    cmd3txtAx=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',{''Un-aged''},3,''V2:V2'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'',data.final.table.Len_Age_Matrix_AcoustALL,3,''A2:V42'');'];
    
    %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1txtAx)
    eval(cmd1)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-aged'},1,'V2');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},1,'W2');
    sum_over_lenM=sum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_ageM=sum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_ageM,1,'W3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},1,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},1,'A45');
    nM_aged=sum(sum_over_lenM(1:end-1));              % total number of aged fish
    p1=sum_over_lenM(1)/nM_aged;                      % age1 proportion
    tot_age2_M=(nM_aged+sum_over_lenM(end))*(1-p1);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],tot_age2_M,1,'B45');
%     xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM(2:end)),1,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Over Age:'},1,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_ageM),1,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-kriged Acoustically Weighted Abundance (Male)'},1,'I47');
    
    %% Female    
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2txtAx)
    eval(cmd2)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-aged'},2,'V2');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},2,'W2');
    sum_over_lenF=sum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_ageF=sum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_ageF,2,'W3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},2,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},2,'A45');
    nF_aged=sum(sum_over_lenF(1:end-1));              % total number of aged fish
    p1=sum_over_lenF(1)/nF_aged;                      % age1 proportion
    tot_age2_F=(nF_aged+sum_over_lenF(end))*(1-p1);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],tot_age2_F,2,'B45');
%     xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenF(2:end)),2,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Over Age:'},2,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_ageF),2,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-kriged Acoustically Weighted Abundance (Female)'},2,'I47');
    
    %% ALL    
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3txtAx)
    eval(cmd3)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-aged'},3,'V2');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Subtotal'},3,'W2');
    sum_over_len=sum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_len,3,'B43');
    sum_over_age=sum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum_over_age,3,'W3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},3,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},3,'A45');
    nALL_aged=sum(sum_over_len(1:end-1));              % total number of aged fish
    p1=sum_over_len(1)/nALL_aged;                      % age1 proportion
    tot_age2_ALL=(nALL_aged+sum_over_len(end))*(1-p1);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],tot_age2_ALL,3,'B45');
%     xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_len(2:end)),3,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Over Age:'},3,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_age),3,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Male+Female:'},3,'G44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'H44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_abundance_table.xlsx'],{'Un-kriged Acoustically Weighted Abundance (ALL)'},3,'I47');
    
 
    %% write acoustically weighted kriged lenght-age-gender abundance tables
    disp('write acoustically weighted kriged lenght-age-gender abundance tables ...')
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Age ''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',data.final.table.kriged_Num_Len_Age_Matrix_AcoustM,1,''A2:V42'');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Age ''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',data.final.table.kriged_Num_Len_Age_Matrix_AcoustF,2,''A2:V42'');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',{''Age ''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'',data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL,3,''A2:V42'');'];
    
    %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Un-aged'},1,'V2');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},1,'W2');
    sum_over_lenM=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_ageM=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_ageM,1,'W3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},1,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},1,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM(2:end)),1,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Over Age:'},1,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_ageM),1,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Kriged Acoustically Weighted Abundance (Male)'},1,'I47');
    
    %% Female    
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Un-aged'},2,'V2');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},2,'W2');
    sum_over_lenF=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_ageF=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_ageF,2,'W3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},2,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},2,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenF(2:end)),2,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Over Age:'},2,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_ageF),2,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Kriged Acoustically Weighted Abundance (Female)'},2,'I47');
    
    %% ALL
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'un-aged'},3,'V2');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Subtotal'},3,'W2');
    sum_over_len=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_len,3,'B43');
    sum_over_age=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum_over_age,3,'W3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age1+)'},3,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Total (age2+)'},3,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_len(2:end)),3,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Over Age:'},3,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_age),3,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Male+Female:'},3,'G44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'H44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_abundance_table.xlsx'],{'Kriged Acoustically Weighted Abundance (ALL)'},3,'I47');
    
    %% write acoustically weighted un-kriged lenght-age-gender biomass tables
    disp('write acoustically weighted un-kriged lenght-age-gender biomass tables ...')
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Age ''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',data.final.table.Wgt_Len_Age_Matrix_AcoustM,1,''A2:U42'');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Age ''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',data.final.table.Wgt_Len_Age_Matrix_AcoustF,2,''A2:U42'');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',{''Age ''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'',data.final.table.Wgt_Len_Age_Matrix_AcoustALL,3,''A2:U42'');'];
    
     %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},1,'V2');
    sum_over_lenM=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_ageM=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_ageM,1,'V3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},1,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},1,'A45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM(2:end)),1,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Over Age:'},1,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_ageM),1,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Un-kriged Acoustically Weighted Biomass (mmt) (Male)'},1,'I47');
    
    %% Female
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},2,'V2');
    sum_over_lenF=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_ageF=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_ageF,2,'V3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},2,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},2,'A45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenF(2:end)),2,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Over Age:'},2,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_ageF),2,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Un-kriged Acoustically Weighted Biomass (mmt) (Female)'},2,'I47');
    
    %% ALL
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Subtotal'},3,'V2');
    sum_over_len=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_len,3,'B43');
    sum_over_age=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum_over_age,3,'V3');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},3,'A44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},3,'A45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_len(2:end)),3,'B45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Over Age:'},3,'D44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_age),3,'E44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Male+Female:'},3,'G44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'H44');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Male+Female:'},3,'G45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM(2:end))+sum(sum_over_lenF(2:end)),3,'H45');
    xlswrite([para.proc.output_filepath '\un-kriged_len_age_biomass_table.xlsx'],{'Un-kriged Acoustically Weighted Biomass (mmt) (ALL)'},3,'I47');
    
    %% write acoustically weighted kriged lenght-age-gender biomass tables
    disp('write acoustically weighted kriged lenght-age-gender biomass tables ...')
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Age ''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM,1,''A2:U42'');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Age ''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF,2,''A2:U42'');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',{''Age ''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'',data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL,3,''A2:U42'');'];
    
    %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},1,'V2');
    sum_over_lenM=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_ageM=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_ageM,1,'V3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},1,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},1,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM(2:end)),1,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Over Age:'},1,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_ageM),1,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Kriged Acoustically Weighted Biomass (mmt) (Male)'},1,'I47');
    
    %% Female
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},2,'V2');
    sum_over_lenF=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_ageF=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_ageF,2,'V3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},2,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},2,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenF(2:end)),2,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Over Age:'},2,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_ageF),2,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Kriged Acoustically Weighted Biomass (mmt) (Female)'},2,'I47');
    
    %% ALL
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Subtotal'},3,'V2');
    sum_over_len=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_len,3,'B43');
    sum_over_age=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum_over_age,3,'V3');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age1+)'},3,'A44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Total (age2+)'},3,'A45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_len(2:end)),3,'B45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Over Age:'},3,'D44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_age),3,'E44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Male+Female:'},3,'G44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'H44');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Male+Female:'},3,'G45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],sum(sum_over_lenM(2:end))+sum(sum_over_lenF(2:end)),3,'H45');
    xlswrite([para.proc.output_filepath '\kriged_len_age_biomass_table.xlsx'],{'Kriged Acoustically Weighted Biomass (mmt) (ALL)'},3,'I47');
    
    if para.bio_data_type ~= 3 % Not A-Shop data
    %% write total length-haul counts tables
    disp('write total length-haul counts tables ...')
    Letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    [m,n]=size(data.final.table.len_haul_M);
    n0=floor(n/26);
    n1=floor((n-1)/26);
    if n0 == 0
        rangeL_text=[Letters(n+1) '2']
        rangeL_value=[Letters(n+1) '3'];
    else
        rangeL_text=[Letters(n0) Letters(n-n0*26+1) '2'];
        rangeL_value=[Letters(n0) Letters(n-n0*26+1) '3'];
    end
    switch n1
        case 0
            rangeL=['A2:' Letters(n) num2str(m+1)];
        case 1
            rangeL=['A2:A' Letters(n-n1*26) num2str(m+1)];
        case 2
            rangeL=['A2:B' Letters(n-n1*26) num2str(m+1)];
        case 3
            rangeL=['A2:C' Letters(n-n1*26) num2str(m+1)];
        case 4
            rangeL=['A2:D' Letters(n-n1*26) num2str(m+1)];
        case 5
            rangeL=['A2:E' Letters(n-n1*26) num2str(m+1)];
    end
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Haul Number (Male)''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',data.final.table.len_haul_M,1,''' rangeL ''');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Haul Number (Female)''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',data.final.table.len_haul_F,2,''' rangeL ''');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',{''Haul Number (ALL)''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\total_len_haul_counts_table.xlsx'',data.final.table.len_haul_ALL,3,''' rangeL ''');'];
    
    %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end

    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},1,rangeL_text);
    sum_over_lenM=sum(data.final.table.len_haul_M(2:end,2:end));
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_haulM=sum(data.final.table.len_haul_M(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_haulM,1,rangeL_value);    
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Total'},1,'A44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Un-Aged Length-Haul Counts (Male)'},1,'I47');
    
    %% Female
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
%     [m,n]=size(data.final.table.len_haul_F);
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end

    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},2,rangeL_text);
    sum_over_lenF=sum(data.final.table.len_haul_F(2:end,2:end));
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_haulF=sum(data.final.table.len_haul_F(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_haulF,2,rangeL_value);    
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Total'},2,'A44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Un-Aged Length-Haul Counts (Female)'},2,'I47');
    
    %% ALL
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
%     [m,n]=size(data.final.table.len_haul_ALL);
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Subtotal'},3,rangeL_text);
    sum_over_len=sum(data.final.table.len_haul_ALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_len,3,'B43');
    sum_over_haul=sum(data.final.table.len_haul_ALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum_over_haul,3,rangeL_value);    
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Total'},3,'A44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Male+Female'},3,'D44');
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'E44');    
    xlswrite([para.proc.output_filepath '\total_len_haul_counts_table.xlsx'],{'Un-Aged Length-Haul Counts (ALL)'},3,'I47');
    
    %% write aged length-haul counts tables
    disp('write aged length-haul counts tables ...')
    Letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    [m,n]=size(data.final.table.aged_len_haul_M);
    n0=floor(n/26);
    n1=floor((n-1)/26);
    if n0 == 0
        rangeL_text=[Letters(n+1) '2']
        rangeL_value=[Letters(n+1) '3'];
    else
        rangeL_text=[Letters(n0) Letters(n-n0*26+1) '2'];
        rangeL_value=[Letters(n0) Letters(n-n0*26+1) '3'];
    end
    switch n1
        case 0
            rangeL=['A2:' Letters(n) num2str(m+1)];
        case 1
            rangeL=['A2:A' Letters(n-n1*26) num2str(m+1)];
        case 2
            rangeL=['A2:B' Letters(n-n1*26) num2str(m+1)];
        case 3
            rangeL=['A2:C' Letters(n-n1*26) num2str(m+1)];
        case 4
            rangeL=['A2:D' Letters(n-n1*26) num2str(m+1)];
        case 5
            rangeL=['A2:E' Letters(n-n1*26) num2str(m+1)];
    end
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Haul Number (Male)''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',data.final.table.aged_len_haul_M,1,''' rangeL ''');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Haul Number (Female)''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',data.final.table.aged_len_haul_F,2,''' rangeL ''');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',{''Haul Number (ALL)''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'',data.final.table.aged_len_haul_ALL,3,''' rangeL ''');'];
    
    %% Male
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},1,'A43');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},1,rangeL_text);
    sum_over_lenM=sum(data.final.table.aged_len_haul_M(2:end,2:end));
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_lenM,1,'B43');
    sum_over_haulM=sum(data.final.table.aged_len_haul_M(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_haulM,1,rangeL_value);    
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Total'},1,'A44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum(sum_over_lenM),1,'B44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Aged Length-Haul Counts (Male)'},1,'I47');
    
    %% Female
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
    [m,n]=size(data.final.table.aged_len_haul_F);
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},2,'A43');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},2,rangeL_text);
    sum_over_lenF=sum(data.final.table.aged_len_haul_F(2:end,2:end));
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_lenF,2,'B43');
    sum_over_haulF=sum(data.final.table.aged_len_haul_F(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_haulF,2,rangeL_value);    
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Total'},2,'A44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum(sum_over_lenF),2,'B44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Aged Length-Haul Counts (Female)'},2,'I47');
    
    %% ALL
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
%     [m,n]=size(data.final.table.aged_len_haul_ALL);
%     n2=floor((n+1)/26);
%     rangeL_text=[Letters(n+1-n2*26) '2'];
%     if n2 > 0
%         rangeL_text=[Letters(n2) Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n2) Letters(n+1-n2*26) '3'];
%     else
%         rangeL_text=[Letters(n+1-n2*26) '2'];
%         rangeL_value=[Letters(n+1-n2*26) '3'];
%     end
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},3,'A43');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Subtotal'},3,rangeL_text);
    sum_over_len=sum(data.final.table.aged_len_haul_ALL(2:end,2:end));
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_len,3,'B43');
    sum_over_haul=sum(data.final.table.aged_len_haul_ALL(2:end,2:end),2);
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum_over_haul,3,rangeL_value);    
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Total'},3,'A44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum(sum_over_len),3,'B44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Male+Female'},3,'D44');
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],sum(sum_over_lenM)+sum(sum_over_lenF),3,'E44');    
    xlswrite([para.proc.output_filepath '\aged_len_haul_counts_table.xlsx'],{'Aged Length-Haul Counts (ALL)'},3,'I47');
    
    %% write compact length-haul tables
    disp('write compact length-haul tables ...')
    [m,n]=size(data.final.table.compact_len_haul_M);
    n1=floor((n-1)/26);
    switch n1
        case 0
            rangeL=['A2:' Letters(n) num2str(m+1)];
        case 1
            rangeL=['A2:A' Letters(n-n1*26) num2str(m+1)];
        case 2
            rangeL=['A2:B' Letters(n-n1*26) num2str(m+1)];
        case 3
            rangeL=['A2:C' Letters(n-n1*26) num2str(m+1)];
        case 4
            rangeL=['A2:D' Letters(n-n1*26) num2str(m+1)];
        case 5
            rangeL=['A2:E' Letters(n-n1*26) num2str(m+1)];
    end
    cmd1txtL=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Length (cm)''},1,''A1:A1'');'];
    cmd1txtA=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Haul Number (Male)''},1,''L1:L1'');'];
    cmd1=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',data.final.table.compact_len_haul_M,1,''' rangeL ''');'];
    cmd2txtL=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Length (cm)''},2,''A1:A1'');'];
    cmd2txtA=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Haul Number (Female)''},2,''L1:L1'');'];
    cmd2=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',data.final.table.compact_len_haul_F,2,''' rangeL ''');'];
    cmd3txtL=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Length (cm)''},3,''A1:A1'');'];
    cmd3txtA=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',{''Haul Number (ALL)''},3,''L1:L1'');'];
    cmd3=['xlswrite(''' para.proc.output_filepath '\compact_len_haul_table.xlsx'',data.final.table.compact_len_haul_ALL,3,''' rangeL ''');'];
    eval(cmd1txtL)
    eval(cmd1txtA)
    eval(cmd1)
    xlswrite([para.proc.output_filepath '\compact_len_haul_table.xlsx'],{'Compact Length-Haul Counts (Male)'},1,'I44');
    eval(cmd2txtL)
    eval(cmd2txtA)
    eval(cmd2)
    xlswrite([para.proc.output_filepath '\compact_len_haul_table.xlsx'],{'Compact Length-Haul Counts (Female)'},2,'I44');
    eval(cmd3txtL)
    eval(cmd3txtA)
    eval(cmd3)
    xlswrite([para.proc.output_filepath '\compact_len_haul_table.xlsx'],{'Compact Length-Haul Counts (ALL)'},3,'I44');
    end % not write A-Shop data  
    %% write biomass density, NASC, number density table: input  for kriging 
    disp('write biomass output table for kriging ...')
    var=data.final.table.biomass(:,[5 6 21 9 18]);
    krig_input_header={'Lat','Lon','Biomass density','NASC','Number density'};
    xlswrite([para.proc.output_filepath '\kriging_input.xlsx'],krig_input_header,1,'A1');
    xlswrite([para.proc.output_filepath '\kriging_input.xlsx'],var,1,'A2');
    
    %% write un-kriged biomass table output for record with zeros
    disp('write un-kriged biomass table output for record with zeros ...')
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_output-' date '_0.xlsx'],data.final.table.biomass_description,1);
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_output-' date '_0.xlsx'],data.final.table.biomass,1,'A2');
    
    %% write un-kriged biomass table output for record with zeros removed
    disp('write un-kriged biomass table output for record with zeros removed ...')
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_output-' date '_1.xlsx'],data.final.table.biomass_description,1);
    ind=find(data.final.table.biomass(:,15) > eps);
    out_dat=data.final.table.biomass(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_output-' date '_1.xlsx'],out_dat,1,'A2');
    
    %% write kriged biomass table output for record with zeros
    disp('write kriged biomass table output for record with zeros ...')
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_output-' date '_0.xlsx'],data.final.table.kriged_biomass0_description,1);
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_output-' date '_0.xlsx'],data.final.table.kriged_biomass0,1,'A2');
    
    %% write kriged biomass table output for record with zeros removed
    disp('write kriged biomass table output for record with zeros removed ...')
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_output-' date '_1.xlsx'],data.final.table.kriged_biomass0_description,1);
    ind=find(data.final.table.kriged_biomass0(:,7) > eps);
    out_dat=data.final.table.kriged_biomass0(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_output-' date '_1.xlsx'],out_dat,1,'A2');
    
end
generate_aged_reports(1)

if get(hdl.radio_oracle,'value') == 1
    disp('oracle!')
end
toc
return