function save_visualization_results
% save visualization results (xlsx file)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Modification:       4/4/2013
%% Modification        4/18/2018  - add "length vs age" pahes 2 &3 outputs


global data para hdl

if get(hdl.vis.radio_visual_un_kriged,'value') == 1
    data_type_str='unkriged_';
elseif get(hdl.vis.radio_visual_kriged,'value') == 1
    data_type_str='kriged_';
elseif get(hdl.vis.radio_visual_biological,'value') == 1
    data_type_str='biological_';
end
filename0=[data_type_str data.visual.func_str '_vs_' data.visual.var_str '.xlsx'];

ind_beg=find(filename0 =='(');
ind_end=find(filename0 ==')');
for i=length(ind_beg):-1:1
    filename0(ind_beg(i)-1:ind_end(i))=[];
end

ind=find(filename0 == ' ');
filename0(ind)='_';

ind=find(filename0 == '<' | filename0 == '>');
filename0(ind)='';

if isfield(para,'output_dir')
    cd(para.output_dir)
else
    cd(para.home_dir)
end
    
[filename, pathname, filterindex] = uiputfile( ...
{'*.xlsx','Excel Workbook (*.xlsx)';
 '*.xls', 'Excel 97-2003 Workbook (*.xls)';...
 '*.csv','CSV (Comma Delimited (*.csv)'},...
 'Save Output File as',filename0);


if ~isstr(filename)
    return
else
    para.output_dir=pathname;
    cd(para.home_dir)
end

[pathstr, name, ext] = fileparts(filename);

full_fname=fullfile(pathname,filename);



if get(hdl.vis.radio_visual_var_len_age,'value') == 1
    xlswrite(full_fname,{data.visual.func1_str},1,'J1');
    if size(data.visual.func1,1) > length(data.visual.var1)  % un-kriged or kriged data (sex specific)
        %% male
        xlswrite(full_fname,{data.visual.var1_str},1,'A2');
        xlswrite(full_fname,{[data.visual.var2_str  ' (Male) ']},1,'J2');
        xlswrite(full_fname,' ',1,'A3');
        xlswrite(full_fname,data.visual.var1,1,'B3:U3');
        xlswrite(full_fname,data.visual.var2,1,'A4');
        func1_aged_male=data.visual.func1(1:40,:);
        xlswrite(full_fname,func1_aged_male,1,'B4');
        %% female (row offset = 43)
        xlswrite(full_fname,{[data.visual.var2_str  ' (Female) ']},1,'J45');
        xlswrite(full_fname,{data.visual.var1_str},1,'A45');
        xlswrite(full_fname,0,1,'A46');
        xlswrite(full_fname,data.visual.var1,1,'B46:U46');
        xlswrite(full_fname,data.visual.var2,1,'A47');
        func1_aged_male=data.visual.func1(41:80,:);
        xlswrite(full_fname,func1_aged_male,1,'B47');
        %% ALL (row offset = 86)
        xlswrite(full_fname,{data.visual.var1_str},1,'A88');
        xlswrite(full_fname,{[data.visual.var2_str  ' (ALL) ']},1,'J88');
        xlswrite(full_fname,' ',1,'A89');
        xlswrite(full_fname,data.visual.var1,1,'B89:U89');
        xlswrite(full_fname,data.visual.var2,1,'A90');
        func1_aged_male=data.visual.func1(81:120,:);
        xlswrite(full_fname,func1_aged_male,1,'B90');
        if get(hdl.vis.radio_visual_func_number,'value') == 1
            %% Un-aged Abundance
            %% Male
            xlswrite(full_fname,{'Un-Aged (x1000)'},1,'V3');
            func1_unaged_male=data.visual.func2(:,1);
            xlswrite(full_fname,func1_unaged_male,1,'V4');
            %% Female
            xlswrite(full_fname,{'Un-Aged (x1000)'},1,'V46');
            func1_unaged_female=data.visual.func2(:,2);
            xlswrite(full_fname,func1_unaged_female,1,'V47');
            %% ALL
            xlswrite(full_fname,{'Un-Aged (x1000)'},1,'V89');
            func1_unaged_ALL=data.visual.func2(:,3);
            xlswrite(full_fname,func1_unaged_ALL,1,'V90');
        end
    else  % biological (non-sex-specific)
        %% ALL
        xlswrite(full_fname,{data.visual.var1_str},1,'A2');
        xlswrite(full_fname,{[data.visual.var2_str  ' (ALL) ']},1,'J2');
        xlswrite(full_fname,' ',1,'A3');
        xlswrite(full_fname,data.visual.var1',1,'B3:U3');
        xlswrite(full_fname,data.visual.var2,1,'A4');
        xlswrite(full_fname,data.visual.func1',1,'B4');
    end
elseif get(hdl.vis.radio_visual_biological,'value') ~= 1   ...
        |  (get(hdl.vis.radio_visual_biological,'value') == 1 & get(hdl.vis.radio_visual_func_number,'value') == 1) ...
        |  (get(hdl.vis.radio_visual_biological,'value') == 1 & get(hdl.vis.radio_visual_func_biomass,'value') == 1)
    if size(data.visual.func1,2) > 1
        if size(data.visual.func1_str,2) ~= 2
            data.visual.func1_str={{data.visual.func1_str},{'(Male, Female, ALL)'}};
        end
    else
        if ~iscell(data.visual.func1_str)
            data.visual.func1_str={data.visual.func1_str};
        end
    end
    if get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
        xlswrite(full_fname,{data.visual.var1_str},1,'A1');
        xlswrite(full_fname,{data.visual.var2_str},1,'B1');
        if size(data.visual.func1_str,2) > 1
            xlswrite(full_fname,data.visual.func1_str{1},1,'C1');
            xlswrite(full_fname,data.visual.func1_str{2},1,'D1');
        else
            xlswrite(full_fname,data.visual.func1_str,1,'C1');
        end
        ind=find(isnan(data.visual.func1) == 1);
        data.visual.func1(ind)=0;
        xlswrite(full_fname,data.visual.var1,1,'A2');
        xlswrite(full_fname,data.visual.var2,1,'B2');
        xlswrite(full_fname,data.visual.func1,1,'C2');
    else
        xlswrite(full_fname,{data.visual.var1_str},1,'A1');
        if size(data.visual.func1_str,2) > 1
            xlswrite(full_fname,data.visual.func1_str{1},1,'B1');
            xlswrite(full_fname,data.visual.func1_str{2},1,'C1');
        else
            xlswrite(full_fname,data.visual.func1_str,1,'B1');
        end
        ind=find(isnan(data.visual.func1) == 1);
        data.visual.func1(ind)=0;
        xlswrite(full_fname,data.visual.var1,1,'A2');
        xlswrite(full_fname,data.visual.func1,1,'B2');
        if get(hdl.vis.radio_visual_var_age,'value') == 1 & get(hdl.vis.radio_visual_func_number,'value') == 1
            xlswrite(full_fname,{'Un-Aged'},1,'A22');
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            xlswrite(full_fname,{'Male'},1,'A2');
            xlswrite(full_fname,{'Female'},1,'A3');
            xlswrite(full_fname,{'Un-Aged'},1,'A4');
        end
    end
elseif get(hdl.vis.radio_visual_biological,'value') == 1    % biological
    if get(hdl.vis.radio_visual_func_length,'value') == 1 | get(hdl.vis.radio_visual_func_age,'value') == 1
        if get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            xlswrite(full_fname,{data.visual.var1_str},1,'A1');
            xlswrite(full_fname,{data.visual.var2_str},1,'B1');
            xlswrite(full_fname,{data.visual.func1_str},1,'C1');
            xlswrite(full_fname,data.visual.var1,1,'A2');
            xlswrite(full_fname,data.visual.var2,1,'B2');
            xlswrite(full_fname,data.visual.func1,1,'C2');
        else
            if get(hdl.vis.radio_visual_var_weight,'value') ~= 1 & get(hdl.vis.radio_visual_var_survey_region,'value') ~= 1
                xlswrite(full_fname,{'Counts'},1,'A1');
                xlswrite(full_fname,{data.visual.func1_str},1,'A2');
                xlswrite(full_fname,{data.visual.var1_str},1,'B2');
                xlswrite(full_fname,' ',1,'A3');
                xlswrite(full_fname,data.visual.var2,1,'B3');
                xlswrite(full_fname,data.visual.var1,1,'A4');
                xlswrite(full_fname,data.visual.var1,1,'A4');
                xlswrite(full_fname,data.visual.func1,1,'B4');
                if get(hdl.vis.radio_visual_func_age,'value') == 1 & get(hdl.vis.radio_visual_var_length,'value') == 1
                    %% age vs length
                    xlswrite(full_fname,{'Age'},2,'A1');
                    xlswrite(full_fname,data.visual.len_wgt.age_bin(:),2,'A2');
                    xlswrite(full_fname,{'Length (cm)'},2,'B1');
                    xlswrite(full_fname,data.visual.len_wgt.yfit(:),2,'B2');
                    xlswrite(full_fname,{'Regression P'},3,'A1');
                    xlswrite(full_fname,data.visual.len_wgt.p,3,'A2');
                elseif get(hdl.vis.radio_visual_func_length,'value') == 1 & get(hdl.vis.radio_visual_var_age,'value') == 1
                    %% length vs age
                    xlswrite(full_fname,{'Age'},2,'A1');
                    xlswrite(full_fname,data.visual.len_wgt.age_bin(:),2,'A2');
                    xlswrite(full_fname,{'Length (cm)'},2,'B1');
                    xlswrite(full_fname,data.visual.len_wgt.yfit(:),2,'B2');
                    xlswrite(full_fname,{'Regression P'},3,'A1');
                    xlswrite(full_fname,data.visual.len_wgt.p,3,'A2');
                    
                end
            else
                % weight vs length
                xlswrite(full_fname,{data.visual.var1_str},1,'A1');
                xlswrite(full_fname,{data.visual.func1_str},1,'B1');
                xlswrite(full_fname,data.visual.var1,1,'A2');
                xlswrite(full_fname,data.visual.func1,1,'B2');
                xlswrite(full_fname,{data.visual.var1_str},2,'A1');
                xlswrite(full_fname,{data.visual.func1_str},2,'B1');
                xlswrite(full_fname,data.bio.len_wgt_all.len(:),2,'A2');
                xlswrite(full_fname,data.bio.len_wgt_all.wgt(:),2,'B2');
                xlswrite(full_fname,{'w0 and p'},3,'A1');
                xlswrite(full_fname,data.bio.len_wgt_all.reg_w0,3,'A2');
                xlswrite(full_fname,data.bio.len_wgt_all.reg_p,3,'A3');                
            end
        end
    else
        xlswrite(full_fname,{data.visual.var1_str},1,'A1');
        xlswrite(full_fname,{data.visual.var2_str},1,'B1');
    end
end