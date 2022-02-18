function		file_browser3d(window,index,type)
%% function		file_browser3d(window,index) handles all brower pushbutton selection
% window_index = index of process task 
%            1 - Data Preparation
%            2 - Variogram
%            3 - Krig
%            2 - Visualization
% index = index for specific options
% type  = Load/Save parameter data type, only for Krig option
%            1 - Varigram parameters only
%            2 - Kriging parameters only
%            3 - Both of 1 and 2
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

HDIR=para.krig_home_dir;

switch  window
    case 1							% data preparation window
        if index == 1		% data input conversion filename
            eval(['cd ''' para.file_dir.data_conversion ''''])
            if para.platform == 1
                [filename, filepath]=uigetfile('*.m; *.dll','Select A Data-Conversion File');
            elseif para.platform == 2
                [filename, filepath]=uigetfile('*','Select A Data-Conversion File');
            elseif para.platform == 3
                [filename, filepath]=uigetfile('*.m','Select A Data-Conversion File');
            end
            eval(['cd ''' HDIR ''''])
            if isstr(filename)				% file exists
                [pathi, truncated_filename, ext]=fileparts(filename);
                para.dataprep.dat_conv_fname=truncated_filename;
                set(hdl.dataprep.dat_conv_fname,'string',filename);
                para.file_dir.data_conversion=filepath;
                cmd=['addpath ' filepath];
                eval(cmd);
            end
        elseif index == 2  % input data filename
            eval(['cd ''' para.file_dir.datafile ''''])
            if para.platform == 1
                [filename, filepath]=uigetfile('*.dat; *.txt','Select an Input File');
            elseif para.platform == 2
                [filename, filepath]=uigetfile('*','Select an Input File');
            elseif para.platform == 3
                [filename, filepath]=uigetfile('*.dat','Select an Input File');
            end
            eval(['cd '''  HDIR ''''])
            if isstr(filename)				% file exists
                para.dataprep.filename=[filepath filename];
                set(hdl.dataprep.file,'string',para.dataprep.filename);
                [pathi, truncated_filename, ext]=fileparts(para.dataprep.filename);
                para.dataprep.fileID=truncated_filename;
                set(hdl.dataprep.fileID,'string',para.dataprep.fileID);
                dataprep3d(1)
                para.file_dir.datafile=filepath;
            end
        elseif index == 3  % save data format filename
            if isfield(data,'in')
                save_data_format_info(1)
            else
                krig_message(1,'Data have not been loaded yet !!');
            end
        end
    case 2									% fit variogram/correlogram window
        eval(['cd ''' para.file_dir.parafile ''''])
        switch index
            case 1                             % load a parameter file
                [filename, filepath]=uigetfile('*.mat','Select a Model Parameter File');
                if isstr(filename)
                    fname=[filepath filename];
                    para.vario.para_file=fname;
                    load_para_file(fname,1);					% new variogram/correlogram window parameters
                    set3dvariopara(2);							% set both edit field and slider positions
                    para.file_dir.parafile=filepath;
                end
            case 2                                % save a parameter file
                [pathi, truncated_filename, ext]=fileparts(para.dataprep.filename);
                para.dataprep.para_file=[truncated_filename '_para'];
                [filename, filepath]=uiputfile('*.mat','Save the Model Parameter File as',para.dataprep.para_file);
                if isstr(filename)
                    fname=[filepath filename];
                    para.vario.para_file=fname;
                    save_para_file(fname,1);            % save parameters only
                    para.file_dir.parafile=filepath;
                end
        end
        eval(['cd ''' HDIR ''''])
    case 3									% kriging window
        switch index
            case 1			 % parameter file
                if get(hdl.krig.vario_para,'value') == 1		% variogram parameter file
                    para_type=1;
                elseif get(hdl.krig.krig_para,'value') == 1	% kriging parameter file
                    para_type=2;
                elseif get(hdl.krig.both_para,'value') == 1	% both variogram and kriging parameter file
                    para_type=3;
                else
                    para_type=3;
                    krig_message(1,'You need to select at least one option !!!');
                    return
                end
                eval(['cd ''' para.file_dir.parafile ''''])
                if get(hdl.krig.load_para,'value') == 1				 % load parameter file
                    [filename, filepath]=uigetfile('*.mat','Select a Parameter File');
                    if isstr(filename)
                        para.krig.para_file_in=[filepath filename];
                        load_para_file(para.krig.para_file_in,para_type);
                        para.krig.save_para=0;
                    else
                        set(hdl.krig.load_para,'value',0);
                        set(hdl.krig.para_file_browser,'enable','off');
                        set(hdl.krig.vario_para,'enable','off');
                        set(hdl.krig.krig_para,'enable','off');
                        set(hdl.krig.both_para,'enable','off');
                        return
                    end
                    if para_type >= 2
                        set3dkrigpara(1);
                    end			% set kriging parameter from variable struct
                elseif get(hdl.krig.save_para,'value')	== 1			 % save parameter file
                    if isfield(para.krig,'batch_file_proc')
                        para.krig.bat_proc_cnt=1;
                    end
                    if ~isfield(para.dataprep,'para_file')
                        if isfield(para.dataprep,'fileID')
                            para.dataprep.para_file=[para.dataprep.fileID '_para'];
                        else
                            para.dataprep.para_file=pwd;
                        end
                    end
                    [filename, filepath]=uiputfile('*.mat','Save to a Parameter File',para.dataprep.para_file);
                    if isstr(filename)
                        para.krig.para_file_out=[filepath filename];
                        save_para_file(para.krig.para_file_out,1);		% save parameters only
                    else
                        set(hdl.krig.para_file_browser,'enable','off');
                        set(hdl.krig.save_para,'value',0);
                        return
                    end
                end
                eval(['cd ''' HDIR ''''])
                para.file_dir.parafile=filepath;
            case 2										% load data file
                pp=get(hdl.krig.load_data_file,'value');
                eval(['cd ''' para.file_dir.datafile ''''])
                [filename, filepath]=uigetfile('*.dat; *.txt','Select an Input Data File');
                eval(['cd ''' HDIR ''''])
                if ~isstr(filename)
                    set(hdl.krig.data_file_browser,'enable','off');
                    set(hdl.krig.load_data_file,'value',0);
                    set(hdl.krig.data_file,'visible','off');
                    return;
                end
                para.krig.load_data_file=1;
                para.file_dir.datafile=filepath;
                para.krig.data_file=[filepath filename];
                set(hdl.krig.data_file,'string',para.krig.data_file);
                [pathi, truncated_filename, ext]=fileparts(para.krig.data_file);
                set(hdl.krig.fileID,'string',truncated_filename);
                para.dataprep.fileID=truncated_filename;
                if para.krig.load_data_format_file == 0
                    loaddatfile(2);			% load data file from kriging window
                    datachk(2);					% calling datachk.m from kriging window using the existing data format
                else
                    load_data_format_info
                    loaddatfile(2);			% load data file from kriging window
                    datachk(3);					% calling datachk.m from kriging window using a specified data format
                end
                if ~isfield(para.dataprep,'transform_index') para.dataprep.transform_index=1;end
                data.in.tv=datatransform(1,data.in.v,para.dataprep.transform_index);	% Forward Data Transformation
                para.krig.load_data_file=1;
                para.dispkrig.Qcheck=0;  % initilize cross validation status
            case 3						% load parameter file from menu bar
                set(hdl.krig.load_para,'value',1);
                radio_action(3,1);
                switch type
                    case 1		% load variogram parameters only
                        set(hdl.krig.vario_para,'value',1);
                        radio_action(3,3);
                    case 2		% load kriging parameters only
                        set(hdl.krig.krig_para,'value',1);
                        radio_action(3,4);
                    case 3		% load both variogram and kriging parameters
                        set(hdl.krig.both_para,'value',1);
                        radio_action(3,5);
                end
                file_browser3d(3,1,1);
            case 4						% save parameter file from menu bar
                set(hdl.krig.save_para,'value',1);
                radio_action(3,2);
                file_browser3d(3,1,2);
            case 5						% load data file from menu bar
                set(hdl.krig.load_data_file,'value',1);
                radio_action(3,6);
                file_browser3d(3,2);
            case 6                      % batch process filename-list file
                eval(['cd ''' para.file_dir.batch_filename '''']);
                [filename, filepath]=uigetfile('*.dat','Select a File');
                eval(['cd ''' HDIR ''''])
                para.file_dir.batch_filename=filepath;
                if ~isstr(filename)
                    return;
                end
                para.krig.batch_data_file=[filepath filename];
            case 7                                     % batch process log file
                eval(['cd ''' para.file_dir.batch_log ''''])
                if isfield(para,'dataprep') & isfield(para.dataprep,'fileID')
                    batch_log_file_path=para.dataprep.fileID;
                else
                    batch_log_file_path='';
                end
                [filename, filepath]=uiputfile('*.log','Define a Log Filename',batch_log_file_path);
                eval(['cd ''' HDIR ''''])
                if ~isstr(filename)
                    return;
                end
                if length(filename) > 4 & strcmp(filename(end-3:end),'.log')
                    para.krig.batch_log_file=[filepath filename];
                elseif isempty(find(filename == '.'))
                    para.krig.batch_log_file=[filepath filename '.log'];
                end
                para.file_dir.batch_log=filepath;
            case 8                                      % load customized grid file
                eval(['cd ''' para.file_dir.gridfile ''''])
                [filename, filepath]=uigetfile('*.dat; *.txt','Define a Customized Grid Filename');
                eval(['cd ''' HDIR ''''])
                if ~isstr(filename)
                    return;
                end
                para.krig.load_griddata_file=1;
                para.file_dir.gridfile=filepath;
                para.krig.grid_file=[filepath filename];
                get_set_gridfile_para
                set3dkrigpara(2)
            case 9                                     % load data format file
                eval(['cd ''' para.file_dir.datafile ''''])
                [filename, filepath]=uigetfile('*.mat','Select an Input Data Format File');
                eval(['cd ''' HDIR ''''])
                para.krig.load_data_format_file=1;
                if ~isstr(filename)
                    set(hdl.krig.data_format_file_browser,'enable','off');
                    set(hdl.krig.load_data_format_file,'value',0);
                    set(hdl.krig.data_file,'visible','off');
                    return;
                end
                para.file_dir.data_format_file=filepath;
                para.krig.data_format_file=[filepath filename];
                set(hdl.krig.data_file,'string',para.krig.data_format_file);
                load_data_format_info
        end
    case 4									% visualization window
        switch index
            case 1% parameter file
                eval(['cd ''' para.file_dir.mat_file_in ''''])
                [filename, filepath]=uigetfile( '*.mat','Load a .mat File');
                eval(['cd ''' HDIR ''''])
                fname=[filepath filename];
                if isstr(filename)
                    loaddatfile(3,fname);			% load data file from visualization window
                    para.file_dir.mat_file_in=filepath;
                end
            case 2								% save krig output to a file
                default_fname=get(hdl.dispkrig3d.fileID,'string');
                indx=find(default_fname == '.')-1;
                if isempty(indx)
                    indx=length(default_fname);
                end
                eval(['cd ''' para.file_dir.mat_file_out ''''])
                [filename, filepath]=uiputfile([default_fname(1:indx) '.mat'],'Save output to a File',para.dataprep.fileID);
                eval(['cd ''' HDIR ''''])
                if isstr(filename)
                    para.save_output.filename=[filepath filename];
                    save_para_file(para.save_output.filename,2);		% save everything
                    para.file_dir.mat_file_out=filepath;
                end
            case 3								% save krig figure to a file
                filemenufcn(gcbf,'FileSaveAs')
        end
        
end										% end window switch
eval(['chdir ' '''' para.home_dir '''']);

