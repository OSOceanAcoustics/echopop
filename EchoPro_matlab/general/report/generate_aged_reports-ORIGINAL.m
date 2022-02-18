function generate_aged_reports(hdl)
% generate aged_abundance/biomass structured reports (excel spreadsheet files, .xlsx files)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       5/7/2016

global data para

tic
%% write acoustically weighted un-kriged lenght-age-gender abundance tables 
% files=dir([para.proc.output_filepath '\*aged_*.xlsx']);
% for i=1:length(files)
%     delete([para.proc.output_filepath '\' files(i).name])
% end
if para.proc.exclude_age1 == 1
    age_ind=2:length(para.bio.hake_age_bin);
else
    age_ind=1:length(para.bio.hake_age_bin);
end
% if get(hdl.radio_assessment,'value') == 1
if hdl == 1
    n=length(data.bio.strata);
    age_str={mat2cell(1:length(para.bio.hake_age_bin),1,length(para.bio.hake_age_bin))};
    age_str={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};
    for i=1:n
        if isempty(data.bio.strata) ~= 1 & isempty(data.bio.strata(i).Len_Age_key_wgt_Mn) ~= 1
            if sum(sum(data.bio.strata(i).Len_Age_key_wgt_Mn(:,age_ind))) == 0
                data.bio.strata(i).num_age_key_wgt_male=sum(data.bio.strata(i).Len_Age_key_wgt_Mn(:,age_ind));
            else
                data.bio.strata(i).num_age_key_wgt_male=sum(data.bio.strata(i).Len_Age_key_wgt_Mn(:,age_ind))/sum(sum(data.bio.strata(i).Len_Age_key_wgt_Mn(:,age_ind)));
            end
            if sum(sum(data.bio.strata(i).Len_Age_key_wgt_Fn(:,age_ind))) == 0
                data.bio.strata(i).num_age_key_wgt_female=sum(data.bio.strata(i).Len_Age_key_wgt_Fn(:,age_ind));
            else
                data.bio.strata(i).num_age_key_wgt_female=sum(data.bio.strata(i).Len_Age_key_wgt_Fn(:,age_ind))/sum(sum(data.bio.strata(i).Len_Age_key_wgt_Fn(:,age_ind)));
            end
            if sum(sum(data.bio.strata(i).Len_Age_key_wgt_ALLn(:,age_ind))) == 0
                data.bio.strata(i).num_age_key_wgt_all=sum(data.bio.strata(i).Len_Age_key_wgt_ALLn(:,age_ind));
            else
                data.bio.strata(i).num_age_key_wgt_all=sum(data.bio.strata(i).Len_Age_key_wgt_ALLn(:,age_ind))/sum(sum(data.bio.strata(i).Len_Age_key_wgt_ALLn(:,age_ind)));
            end
        end
    end
    
    %% write un-kriged aged biomass table output for record with zeros
    disp('write un-kriged aged biomass output for record with zeros ...')
    unkriged_index_M=[1 5:7 13];
    unkriged_index_F=[1 5:7 14];
    unkriged_index_ALL=[1 5:7 15];
    for i=1:size(data.final.table.biomass,1);
       unkriged_aged_wgt_M(i,age_ind)=data.bio.strata(data.final.table.biomass(i,7)).num_age_key_wgt_male*data.final.table.biomass(i,unkriged_index_M(end));
       unkriged_aged_wgt_F(i,age_ind)=data.bio.strata(data.final.table.biomass(i,7)).num_age_key_wgt_female*data.final.table.biomass(i,unkriged_index_F(end));
       unkriged_aged_wgt_ALL(i,age_ind)=data.bio.strata(data.final.table.biomass(i,7)).num_age_key_wgt_all*data.final.table.biomass(i,unkriged_index_ALL(end));
    end
    %% ALL (Male + Female)
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],{'ALL (Male + Female)'},1,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass_description(unkriged_index_ALL),1,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,1,'F2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass(:,unkriged_index_ALL),1,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],unkriged_aged_wgt_ALL,1,'F3');
    %% Male
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],{'Male'},2,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass_description(unkriged_index_M),2,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,2,'F2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass(:,unkriged_index_M),2,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],unkriged_aged_wgt_M,2,'F3');
    %% Female
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],{'Female'},3,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass_description(unkriged_index_F),3,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,3,'F2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.biomass(:,unkriged_index_F),3,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_0.xlsx'],unkriged_aged_wgt_F,3,'F3');
    
    %% write un-kriged aged biomass table output for record with zeros removed
    disp('write un-kriged aged biomass output for record with zeros removed ...')
    %% ALL (Male + Female)
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],{'ALL (Male + Female)'},1,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.biomass_description(unkriged_index_ALL),1,'A2');
    ind=find(data.final.table.biomass(:,unkriged_index_ALL(end)) > eps);              % biomass column for ALL (male+female)
    unkriged_out_dat_ALL=data.final.table.biomass(ind,:);      
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_out_dat_ALL(:,unkriged_index_ALL),1,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,1,'F2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_aged_wgt_ALL(ind,:),1,'F3');
    %% Male
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],{'Male'},2,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.biomass_description(unkriged_index_M),2,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,2,'F2');
    ind=find(data.final.table.biomass(:,unkriged_index_M(end)) > eps);             % biomass column for male only
    unkriged_out_dat_M=data.final.table.biomass(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_out_dat_M(:,unkriged_index_M),2,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_aged_wgt_M(ind,:),2,'F3');
    %% Female
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],{'Female'},3,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.biomass_description(unkriged_index_F),3,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,3,'F2');
    ind=find(data.final.table.biomass(:,unkriged_index_F(end)) > eps);             % biomass column for female only
    unkriged_out_dat_F=data.final.table.biomass(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_out_dat_F(:,unkriged_index_F),3,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_un-kriged_aged_output-' para.survey_year '_1.xlsx'],unkriged_aged_wgt_F(ind,:),3,'F3');
    
    
    
    
    %% write kriged aged biomass table output for record with zeros
    kriged_index_M=[1:3 8];
    kriged_index_F=[1:3 9];
    kriged_index_ALL=[1:3 10];
    for i=1:size(data.final.table.kriged_biomass0,1);
        kriged_aged_wgt_M(i,age_ind)=data.bio.strata(data.final.table.kriged_biomass0(i,3)).num_age_key_wgt_male*data.final.table.kriged_biomass0(i,kriged_index_M(end));
        kriged_aged_wgt_F(i,age_ind)=data.bio.strata(data.final.table.kriged_biomass0(i,3)).num_age_key_wgt_female*data.final.table.kriged_biomass0(i,kriged_index_F(end));
        kriged_aged_wgt_ALL(i,age_ind)=data.bio.strata(data.final.table.kriged_biomass0(i,3)).num_age_key_wgt_all*data.final.table.kriged_biomass0(i,kriged_index_ALL(end));
    end
    disp('write kriged aged biomass table output for record with zeros ...')
    %% ALL (Male + Female)
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],{'ALL (Male + Female)'},1,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_ALL),1,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,1,'E2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0(:,kriged_index_ALL),1,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],kriged_aged_wgt_ALL,1,'E3');
    %% Male
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],{'Male'},2,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_M),2,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,2,'E2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0(:,kriged_index_M),2,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],kriged_aged_wgt_M,2,'E3');
    %% Female
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],{'Female'},3,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_F),3,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],age_str,3,'E2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],data.final.table.kriged_biomass0(:,kriged_index_F),3,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_0.xlsx'],kriged_aged_wgt_F,3,'E3');
    
    %% write kriged aged biomass table output for record with zeros removed
    disp('write kriged aged biomass table output for record with zeros removed ...')
    %% ALL (Male + Female)
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],{'ALL (Male + Female)'},1,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_ALL),1,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,1,'E2');
    ind=find(data.final.table.kriged_biomass0(:,kriged_index_ALL(end)) > eps);     % biomass column for ALL (male+female)
    kriged_out_dat_ALL=data.final.table.kriged_biomass0(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_out_dat_ALL(:,kriged_index_ALL),1,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_aged_wgt_ALL(ind,:),1,'E3');
    %%  Male
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],{'Male'},2,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_M),2,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,2,'E2');
    ind=find(data.final.table.kriged_biomass0(:,kriged_index_M(end)) > eps);     % biomass column for male only
    kriged_out_dat_M=data.final.table.kriged_biomass0(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_out_dat_M(:,kriged_index_M),2,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_aged_wgt_M(ind,:),2,'E3');
    %%  Female
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],{'Female'},3,'A1');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],data.final.table.kriged_biomass0_description(kriged_index_F),3,'A2');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],age_str,3,'E2');
    ind=find(data.final.table.kriged_biomass0(:,kriged_index_F(end)) > eps);     % biomass column for female only
    kriged_out_dat_F=data.final.table.kriged_biomass0(ind,:);
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_out_dat_F(:,kriged_index_F),3,'A3');
    xlswrite([para.proc.output_filepath '\EchoPro_kriged_aged_output-' para.survey_year '_1.xlsx'],kriged_aged_wgt_F(ind,:),3,'E3');
    
end


toc
return