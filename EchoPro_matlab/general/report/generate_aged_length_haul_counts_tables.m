function generate_aged_length_haul_counts_tables
% generate lenght-age-sex structured reports (excel spreadsheet files, .xlsx files)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       12/22/2016

global data para

files=dir([para.proc.output_filepath '\*aged_*.xlsx']);
for i=1:length(files)
    delete([para.proc.output_filepath '\' files(i).name])
end

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