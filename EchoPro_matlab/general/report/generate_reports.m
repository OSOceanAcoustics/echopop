function generate_reports(hdl)
    % generate lenght-age-sex structured reports (excel spreadsheet files, .xlsx files)
    %
    %% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
    %% Last Modification:       3/24/2013
    
    global data para
    
    tic
    %% write acoustically weighted un-kriged lenght-age-gender abundance tables 
    files=dir([para.proc.output_filepath '/*.xlsx']);
    for i=1:length(files)
        delete([para.proc.output_filepath '/' files(i).name])
    end
    if get(hdl.radio_assessment,'value') == 1
        disp('write acoustically weighted un-kriged lenght-age-gender abundance tables ...')
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Age (Male)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1txtAx=['writematrix([''Un-aged''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''V2:V2'');'];
        %   cmd1txtLsum=['writematrix([''Age (Male)''],''' para.proc.output_filepath '/len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.Len_Age_Matrix_AcoustM, ''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''A2:V42'');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Age (Female)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2txtAx=['writematrix([''Un-aged''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''V2:V2'');'];
        cmd2=['writematrix(data.final.table.Len_Age_Matrix_AcoustF, ''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''A2:V42'');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Age (All)''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3txtAx=['writematrix([''Un-aged''],''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''V2:V2'');'];
        cmd3=['writematrix(data.final.table.Len_Age_Matrix_AcoustALL, ''' para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''A2:V42'');'];
        
        %% Male
        eval(cmd1txtL)
        eval(cmd1txtA)
        eval(cmd1txtAx)
        eval(cmd1)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Un-aged'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','W2');
        sum_over_lenM=sum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_ageM=sum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end),2);
        writematrix(sum_over_ageM,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A45');
        nM_aged=sum(sum_over_lenM(1:end-1));              % total number of aged fish
        p1=sum_over_lenM(1)/nM_aged;                      % age1 proportion
        tot_age2_M=(nM_aged+sum_over_lenM(end))*(1-p1);
        writematrix(tot_age2_M,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B45');
    %     writematrix(sum(sum_over_lenM(2:end)),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','D44');
        writematrix(sum(sum_over_ageM),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','E44');
        writematrix(['Un-kriged Acoustically Weighted Abundance (Male)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','I47');
        
        %% Female    
        eval(cmd2txtL)
        eval(cmd2txtA)
        eval(cmd2txtAx)
        eval(cmd2)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Un-aged'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','W2');
        sum_over_lenF=sum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_ageF=sum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end),2);
        writematrix(sum_over_ageF,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A45');
        nF_aged=sum(sum_over_lenF(1:end-1));              % total number of aged fish
        p1=sum_over_lenF(1)/nF_aged;                      % age1 proportion
        tot_age2_F=(nF_aged+sum_over_lenF(end))*(1-p1);
        writematrix(tot_age2_F,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B45');
    %     writematrix(sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','D44');
        writematrix(sum(sum_over_ageF),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','E44');
        writematrix(['Un-kriged Acoustically Weighted Abundance (Female)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','I47');
        
        %% ALL    
        eval(cmd3txtL)
        eval(cmd3txtA)
        eval(cmd3txtAx)
        eval(cmd3)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['Un-aged'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','W2');
        sum_over_len=sum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_age=sum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end),2);
        writematrix(sum_over_age,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A45');
        nALL_aged=sum(sum_over_len(1:end-1));              % total number of aged fish
        p1=sum_over_len(1)/nALL_aged;                      % age1 proportion
        tot_age2_ALL=(nALL_aged+sum_over_len(end))*(1-p1);
        writematrix(tot_age2_ALL,[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B45');
    %     writematrix(sum(sum_over_len(2:end)),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_age),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','E44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','G44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','H44');
        writematrix(['Un-kriged Acoustically Weighted Abundance (ALL)'],[para.proc.output_filepath '/un-kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','I47');
        
     
        %% write acoustically weighted kriged lenght-age-gender abundance tables
        disp('write acoustically weighted kriged lenght-age-gender abundance tables ...')
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM, ''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',1,''Range'',''A2:V42'');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF, ''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',2,''Range'',''A2:V42'');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL, ''' para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'',''Sheet'',3,''Range'',''A2:V42'');'];
        
        %% Male
        eval(cmd1txtL)
        eval(cmd1txtA)
        eval(cmd1)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Un-aged'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','W2');
        sum_over_lenM=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_ageM=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end),2);
        writematrix(sum_over_ageM,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','A45');
        writematrix(sum(sum_over_lenM(2:end)),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','D44');
        writematrix(sum(sum_over_ageM),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','E44');
        writematrix(['Kriged Acoustically Weighted Abundance (Male)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',1,'Range','I47');
        
        %% Female    
        eval(cmd2txtL)
        eval(cmd2txtA)
        eval(cmd2)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Un-aged'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','W2');
        sum_over_lenF=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_ageF=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end),2);
        writematrix(sum_over_ageF,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','A45');
        writematrix(sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','D44');
        writematrix(sum(sum_over_ageF),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','E44');
        writematrix(['Kriged Acoustically Weighted Abundance (Female)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',2,'Range','I47');
        
        %% ALL
        eval(cmd3txtL)
        eval(cmd3txtA)
        eval(cmd3)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['un-aged'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','V2');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','W2');
        sum_over_len=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_age=sum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
        writematrix(sum_over_age,[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','W3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','A45');
        writematrix(sum(sum_over_len(2:end)),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_age),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','E44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','G44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','H44');
        writematrix(['Kriged Acoustically Weighted Abundance (ALL)'],[para.proc.output_filepath '/kriged_len_age_abundance_table.xlsx'],'Sheet',3,'Range','I47');
        
        %% write acoustically weighted un-kriged lenght-age-gender biomass tables
        disp('write acoustically weighted un-kriged lenght-age-gender biomass tables ...')
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.Wgt_Len_Age_Matrix_AcoustM, ''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''A2:U42'');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.Wgt_Len_Age_Matrix_AcoustF, ''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''A2:U42'');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.Wgt_Len_Age_Matrix_AcoustALL, ''' para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''A2:U42'');'];
        
         %% Male
        eval(cmd1txtL)
        eval(cmd1txtA)
        eval(cmd1)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','V2');
        sum_over_lenM=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_ageM=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2);
        writematrix(sum_over_ageM,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A45');
        writematrix(sum(sum_over_lenM(2:end)),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','D44');
        writematrix(sum(sum_over_ageM),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','E44');
        writematrix(['Un-kriged Acoustically Weighted Biomass (mmt) (Male)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','I47');
        
        %% Female
        eval(cmd2txtL)
        eval(cmd2txtA)
        eval(cmd2)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','V2');
        sum_over_lenF=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_ageF=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2);
        writematrix(sum_over_ageF,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A45');
        writematrix(sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','D44');
        writematrix(sum(sum_over_ageF),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','E44');
        writematrix(['Un-kriged Acoustically Weighted Biomass (mmt) (Female)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','I47');
        
        %% ALL
        eval(cmd3txtL)
        eval(cmd3txtA)
        eval(cmd3)
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','V2');
        sum_over_len=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_age=sum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
        writematrix(sum_over_age,[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A45');
        writematrix(sum(sum_over_len(2:end)),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_age),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','E44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','G44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','H44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','G45');
        writematrix(sum(sum_over_lenM(2:end))+sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','H45');
        writematrix(['Un-kriged Acoustically Weighted Biomass (mmt) (ALL)'],[para.proc.output_filepath '/un-kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','I47');
        
        %% write acoustically weighted kriged lenght-age-gender biomass tables
        disp('write acoustically weighted kriged lenght-age-gender biomass tables ...')
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM, ''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',1,''Range'',''A2:U42'');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF, ''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',2,''Range'',''A2:U42'');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Age ''],''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL, ''' para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'',''Sheet'',3,''Range'',''A2:U42'');'];
        
        %% Male
        eval(cmd1txtL)
        eval(cmd1txtA)
        eval(cmd1)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','V2');
        sum_over_lenM=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_ageM=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2);
        writematrix(sum_over_ageM,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','A45');
        writematrix(sum(sum_over_lenM(2:end)),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','D44');
        writematrix(sum(sum_over_ageM),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','E44');
        writematrix(['Kriged Acoustically Weighted Biomass (mmt) (Male)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',1,'Range','I47');
        
        %% Female
        eval(cmd2txtL)
        eval(cmd2txtA)
        eval(cmd2)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','V2');
        sum_over_lenF=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_ageF=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2);
        writematrix(sum_over_ageF,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','A45');
        writematrix(sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','D44');
        writematrix(sum(sum_over_ageF),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','E44');
        writematrix(['Kriged Acoustically Weighted Biomass (mmt) (Female)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',2,'Range','I47');
        
        %% ALL
        eval(cmd3txtL)
        eval(cmd3txtA)
        eval(cmd3)
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','V2');
        sum_over_len=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_age=sum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2);
        writematrix(sum_over_age,[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','V3');
        writematrix(['Total (age1+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Total (age2+)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','A45');
        writematrix(sum(sum_over_len(2:end)),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','B45');
        writematrix(['Over Age:'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_age),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','E44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','G44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','H44');
        writematrix(['Male+Female:'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','G45');
        writematrix(sum(sum_over_lenM(2:end))+sum(sum_over_lenF(2:end)),[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','H45');
        writematrix(['Kriged Acoustically Weighted Biomass (mmt) (ALL)'],[para.proc.output_filepath '/kriged_len_age_biomass_table.xlsx'],'Sheet',3,'Range','I47');
        
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
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Haul Number (Male)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.len_haul_M, ''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''' rangeL ''');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Haul Number (Female)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.len_haul_F, ''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''' rangeL ''');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Haul Number (ALL)''],''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.len_haul_ALL, ''' para.proc.output_filepath '/total_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''' rangeL ''');'];
        
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
    
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range',rangeL_text);
        sum_over_lenM=sum(data.final.table.len_haul_M(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_haulM=sum(data.final.table.len_haul_M(2:end,2:end),2);
        writematrix(sum_over_haulM,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Un-Aged Length-Haul Counts (Male)'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',1,'Range','I47');
        
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
    
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range',rangeL_text);
        sum_over_lenF=sum(data.final.table.len_haul_F(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_haulF=sum(data.final.table.len_haul_F(2:end,2:end),2);
        writematrix(sum_over_haulF,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Un-Aged Length-Haul Counts (Female)'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',2,'Range','I47');
        
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
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range',rangeL_text);
        sum_over_len=sum(data.final.table.len_haul_ALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_haul=sum(data.final.table.len_haul_ALL(2:end,2:end),2);
        writematrix(sum_over_haul,[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Male+Female'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','E44');    
        writematrix(['Un-Aged Length-Haul Counts (ALL)'],[para.proc.output_filepath '/total_len_haul_counts_table.xlsx'],'Sheet',3,'Range','I47');
        
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
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Haul Number (Male)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.aged_len_haul_M, ''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',1,''Range'',''' rangeL ''');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Haul Number (Female)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.aged_len_haul_F, ''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',2,''Range'',''' rangeL ''');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Haul Number (ALL)''],''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.aged_len_haul_ALL, ''' para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'',''Sheet'',3,''Range'',''' rangeL ''');'];
        
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
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range',rangeL_text);
        sum_over_lenM=sum(data.final.table.aged_len_haul_M(2:end,2:end));
        writematrix(sum_over_lenM,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range','B43');
        sum_over_haulM=sum(data.final.table.aged_len_haul_M(2:end,2:end),2);
        writematrix(sum_over_haulM,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range','A44');
        writematrix(sum(sum_over_lenM),[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range','B44');
        writematrix(['Aged Length-Haul Counts (Male)'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',1,'Range','I47');
        
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
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range',rangeL_text);
        sum_over_lenF=sum(data.final.table.aged_len_haul_F(2:end,2:end));
        writematrix(sum_over_lenF,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range','B43');
        sum_over_haulF=sum(data.final.table.aged_len_haul_F(2:end,2:end),2);
        writematrix(sum_over_haulF,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range','A44');
        writematrix(sum(sum_over_lenF),[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range','B44');
        writematrix(['Aged Length-Haul Counts (Female)'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',2,'Range','I47');
        
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
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','A43');
        writematrix(['Subtotal'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range',rangeL_text);
        sum_over_len=sum(data.final.table.aged_len_haul_ALL(2:end,2:end));
        writematrix(sum_over_len,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','B43');
        sum_over_haul=sum(data.final.table.aged_len_haul_ALL(2:end,2:end),2);
        writematrix(sum_over_haul,[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range',rangeL_value);    
        writematrix(['Total'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','A44');
        writematrix(sum(sum_over_len),[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','B44');
        writematrix(['Male+Female'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','D44');
        writematrix(sum(sum_over_lenM)+sum(sum_over_lenF),[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','E44');    
        writematrix(['Aged Length-Haul Counts (ALL)'],[para.proc.output_filepath '/aged_len_haul_counts_table.xlsx'],'Sheet',3,'Range','I47');
        
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
        cmd1txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',1,''Range'',''A1:A1'');'];
        cmd1txtA=['writematrix([''Haul Number (Male)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',1,''Range'',''L1:L1'');'];
        cmd1=['writematrix(data.final.table.compact_len_haul_M, ''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',1,''Range'',''' rangeL ''');'];
        cmd2txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',2,''Range'',''A1:A1'');'];
        cmd2txtA=['writematrix([''Haul Number (Female)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',2,''Range'',''L1:L1'');'];
        cmd2=['writematrix(data.final.table.compact_len_haul_F, ''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',2,''Range'',''' rangeL ''');'];
        cmd3txtL=['writematrix([''Length (cm)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',3,''Range'',''A1:A1'');'];
        cmd3txtA=['writematrix([''Haul Number (ALL)''],''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',3,''Range'',''L1:L1'');'];
        cmd3=['writematrix(data.final.table.compact_len_haul_ALL, ''' para.proc.output_filepath '/compact_len_haul_table.xlsx'',''Sheet'',3,''Range'',''' rangeL ''');'];
        eval(cmd1txtL)
        eval(cmd1txtA)
        eval(cmd1)
        writematrix(['Compact Length-Haul Counts (Male)'],[para.proc.output_filepath '/compact_len_haul_table.xlsx'],'Sheet',1,'Range','I44');
        eval(cmd2txtL)
        eval(cmd2txtA)
        eval(cmd2)
        writematrix(['Compact Length-Haul Counts (Female)'],[para.proc.output_filepath '/compact_len_haul_table.xlsx'],'Sheet',2,'Range','I44');
        eval(cmd3txtL)
        eval(cmd3txtA)
        eval(cmd3)
        writematrix(['Compact Length-Haul Counts (ALL)'],[para.proc.output_filepath '/compact_len_haul_table.xlsx'],'Sheet',3,'Range','I44');
        end % not write A-Shop data  
        %% write biomass density, NASC, number density table: input  for kriging 
        disp('write biomass output table for kriging ...')
        var=data.final.table.biomass(:,[5 6 21 9 18]);
        % krig_input_header={'Lat','Lon','Biomass density','NASC','Number density'};
        krig_input_header=['Lat','Lon','Biomass density','NASC','Number density'];
        writematrix(krig_input_header,[para.proc.output_filepath '/kriging_input.xlsx'],'Sheet',1,'Range','A1');
        writematrix(var,[para.proc.output_filepath '/kriging_input.xlsx'],'Sheet',1,'Range','A2');
        
        %% write un-kriged biomass table output for record with zeros
        disp('write un-kriged biomass table output for record with zeros ...')
        writecell(data.final.table.biomass_description,[para.proc.output_filepath '/EchoPro_un-kriged_output-' date '_0.xlsx'],'Sheet',1);
        writematrix(data.final.table.biomass,[para.proc.output_filepath '/EchoPro_un-kriged_output-' date '_0.xlsx'],'Sheet',1,'Range','A2');
        
        %% write un-kriged biomass table output for record with zeros removed
        disp('write un-kriged biomass table output for record with zeros removed ...')
        writecell(data.final.table.biomass_description,[para.proc.output_filepath '/EchoPro_un-kriged_output-' date '_1.xlsx'],'Sheet',1);
        ind=find(data.final.table.biomass(:,15) > eps);
        out_dat=data.final.table.biomass(ind,:);
        writematrix(out_dat,[para.proc.output_filepath '/EchoPro_un-kriged_output-' date '_1.xlsx'],'Sheet',1,'Range','A2');
        
        %% write kriged biomass table output for record with zeros
        disp('write kriged biomass table output for record with zeros ...')
        writecell(data.final.table.kriged_biomass0_description,[para.proc.output_filepath '/EchoPro_kriged_output-' date '_0.xlsx'],'Sheet',1);
        writematrix(data.final.table.kriged_biomass0,[para.proc.output_filepath '/EchoPro_kriged_output-' date '_0.xlsx'],'Sheet',1,'Range','A2');
        
        %% write kriged biomass table output for record with zeros removed
        disp('write kriged biomass table output for record with zeros removed ...')
        writecell(data.final.table.kriged_biomass0_description,[para.proc.output_filepath '/EchoPro_kriged_output-' date '_1.xlsx'],'Sheet',1);
        ind=find(data.final.table.kriged_biomass0(:,7) > eps);
        out_dat=data.final.table.kriged_biomass0(ind,:);
        writematrix(out_dat,[para.proc.output_filepath '/EchoPro_kriged_output-' date '_1.xlsx'],'Sheet',1,'Range','A2');
        
    end
    generate_aged_reports(1)
    
    if get(hdl.radio_oracle,'value') == 1
        disp('oracle!')
    end
    toc
    return