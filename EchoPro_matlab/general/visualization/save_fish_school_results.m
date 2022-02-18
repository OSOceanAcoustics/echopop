function save_visualization_results
% save visualization results (xlsx file)
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global data para hdl

filename0=['fish_school_' data.visual.school.func1_str '_vs_' data.visual.school.var1_str '.xlsx'];

ind_beg=find(filename0 =='(');
ind_end=find(filename0 ==')');
filename0(ind_beg-1:ind_end)=[];

ind=find(filename0 == ' ');
filename0(ind)='_';


[filename, pathname, filterindex] = uiputfile( ...
{'*.xlsx','Excel Workbook (*.xlsx)';
 '*.xls', 'Excel 97-2003 Workbook (*.xls)';...
 '*.csv','CSV (Comma Delimited (*.csv)';...
 '*.txt','Text (*.txt)';...
 '*.*',  'All Files (*.*)'},...
 'Save Output File as',filename0);

if ~isstr(filename)
    return
end

[pathstr, name, ext] = fileparts(filename);

full_fname=fullfile(pathname,filename);

switch ext
    case {'.xlsx','.xls'}
        xlswrite(full_fname,{data.visual.school.var1_str},1,'A1');
        xlswrite(full_fname,{data.visual.school.func1_str},1,'C1');
        xlswrite(full_fname,data.visual.school.var1,1,'A2');
        xlswrite(full_fname,data.visual.school.func1,1,'C2');        
    case '.csv'
        if size(data.visual.var1,1) == size(data.visual.school.func1,1)
           csvwrite(full_fname,[data.visual.var1 data.visual.school.func1]);
        else
           fprintf('Array dimensions are not consistent!!\n')
        end
    case '.txt'
        if size(data.visual.var1,1) == size(data.visual.school.func1,1)
           data_out=[data.visual.var1 data.visual.school.func1];
           save(full_fname,'data_out','-ascii');
        else
           fprintf('Array dimensions are not consistent!!\n')
        end
    otherwise
end