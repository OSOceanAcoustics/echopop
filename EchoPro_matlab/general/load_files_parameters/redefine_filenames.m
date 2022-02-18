function varargout = redefine_filenames(varargin)
% LOAD_DATA M-file for redefine_filenames.fig
%      LOAD_DATA, by itself, creates a new LOAD_DATA or raises the existing
%      singleton*.
%
%      H = LOAD_DATA returns the handle to a new LOAD_DATA or the handle to
%      the existing singleton*.
%
%      LOAD_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_DATA.M with the given input arguments.
%
%      LOAD_DATA('Property','Value',...) creates a new LOAD_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before redefine_filenames_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to redefine_filenames_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help redefine_filenames

% Last Modified by GUIDE v2.5 04-Jul-2011 01:02:50

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @redefine_filenames_OpeningFcn, ...
                   'gui_OutputFcn',  @redefine_filenames_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before redefine_filenames is made visible.
function redefine_filenames_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to redefine_filenames (see VARARGIN)

global para hdl

% Choose default command line output for redefine_filenames
handles.output = hObject;

% Update handles structure
%% Synchrinize the main panel with the current data update panel
%% set survey year 
set(findobj(hObject,'tag','popup_survey_year'),'value',get(findobj(hdl.main,'tag','popup_survey_year'),'value'));
%% default source of dataset (US, CAN, or US&CAN)
switch para.proc.source
    case 1
         set(findobj(hObject,'Tag','radio_US'),'value',1)
    case 2
         set(findobj(hObject,'Tag','radio_CAN'),'value',1)
    case 3
         set(findobj(hObject,'Tag','radio_ALL'),'value',1)
end
% Update handles structure
guidata(hObject, handles);
return

% --- Outputs from this function are returned to the command line.
function varargout = redefine_filenames_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return

% --- Executes on button press in pb_bio_ength.
function pb_load_update_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_ength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl
hdl.redefine_filenames=handles;
update_proc_parameters(hObject)
return

% --- Executes on button press in pb_load_advance.
function pb_load_advance_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load_advance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl
hdl.load=handles;
return


% --- Executes on button press in pbQuit.
function pbQuit_Callback(hObject, eventdata, handles)
% hObject    handle to pbQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close
return

% --- Executes on button press in pb_bio_gear.
function pb_bio_gear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_gear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Trawl'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end

hdl.load=handles;


% --- Executes on button press in pb_bio_transect_haul.
function pb_bio_trawl_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_transect_haul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Trawl'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end

hdl.load=handles;


% --- Executes on button press in pb_bio_catch.
function pb_bio_catch_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_catch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Catch'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;

return

% --- Executes on button press in pb_bio_ength.
function pb_bio_length_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_ength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Length'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;

return

% --- Executes on button press in pb_bio_pecimen.
function pb_bio_specimen_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_pecimen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl 
[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Specimen'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;
return


% --- Executes on button press in pb_bio_transect_haul.
function pb_bio_transect_haul_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_transect_haul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Transect-Hauls'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end

hdl.load=handles;


% --- Executes on button press in pb_acoust_datatype.
function pb_acoust_datatype_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

if get(handles.radio_acoust_EK60_dir,'value') == 1 | get(handles.radio_acoust_EK500_dir,'value') == 1
    filepath=uigetdir('C:\','Pick an EK60 or EK500 raw data Root Directory');
    if filepath ~= 0
        set(hObject,'UserData',filepath);
    end
else
    [filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)'; ...
    '*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a processed NASC file']);
    if filename ~= 0
        set(hObject,'UserData',[filepath filename]);
    end
end

hdl.load=handles;

return


% --- Executes on button press in pb_acoust_cal.
function pb_acoust_cal_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.txt','Text files (*.txt)';'*.dat','ASCII Files (*.dat)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Calibration'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;


% --- Executes on button press in pb_acoust_evr.
function pb_acoust_evr_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_evr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

filepath=uigetdir('C:\','Pick an EchoView region files Root Directory');
if filepath ~= 0
   set(hObject,'UserData',filepath);
end

hdl.load=handles;


% --- Executes on button press in pb_acoust_transect_region.
function pb_acoust_transect_region_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_transect_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.txt','Text files (*.txt)';'*.dat','ASCII Files (*.dat)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a '' Transect-Region '' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;


% --- Executes on button press in pb_acoust_VL_spacing.
function pb_acoust_VL_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_VL_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.csv','csv files (*.csv)';'*.txt','Text files (*.txt)';'*.dat','ASCII Files (*.dat)';'*.*',' ALL (*.*)'},['Pick a Vessel Log Spacing file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end
hdl.load=handles;


% --- Executes on button press in pb_VL_lat_lon.
function pb_VL_lat_lon_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Vl_lat_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

filepath=uigetdir('C:\','Pick Vessel Log Directory');
if filepath ~= 0
   set(hObject,'UserData',filepath);
end
hdl.load=handles;

% --- Executes on button press in pb_acoust_Ev_bot.
function pb_acoust_Ev_bot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_Ev_bot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl
filepath=uigetdir('C:\','Pick an EchoView Detected Bottom Directory');
if filepath ~= 0
   set(hObject,'UserData',filepath);
end

hdl.load=handles;


% --- Executes on button press in pb_bio_transect_haul.
function pb_bio_strata_Callback(hObject, eventdata, handles)
% hObject    handle to pb_bio_strata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.xls','Excel files (*.xls)';'*.xlsx','Excel Files (*.xlsx)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Strata'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end

hdl.load=handles;

% --- Executes on button press in pb_bio_transect_haul.
function radio_acoust_species_mix_Callback(hObject, eventdata, handles)
% hObject    handle to radio_acoust_species_mix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para

if get(hObject,'value')
    set(handles.pb_acoust_mix_trawl,'Enable','on');
    set(handles.text_acoust_mix_trawl,'Enable','on');
    para.acoust.mix_region=1;
else
    set(handles.pb_acoust_mix_trawl,'Enable','off');    
    set(handles.text_acoust_mix_trawl,'Enable','off');
    para.acoust.mix_region=0;
end

hdl.load=handles;

% --- Executes on button press in radio_acoust_2sta_TS.
function radio_acoust_2sta_TS_Callback(hObject, eventdata, handles)
% hObject    handle to radio_acoust_2sta_TS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para

if get(hObject,'value')
    para.acoust.TS_station_num=2;
else
    para.acoust.TS_station_num=1;
end
hdl.load=handles;


% --- Executes on button press in pb_bio_transect_haul.
function pb_acoust_mix_trawl_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoust_mix_trawl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl

[filename,filepath]=uigetfile({'*.txt','Text Files (*.txt)';'*.xls','Excel files (*.xls)';'*.csv','csv files (*.csv)';'*.*',' ALL (*.*)'},['Pick a ''Strata'' file']);
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
end

hdl.load=handles;

% --- Executes on button press in radio_US.
function radio_US_Callback(hObject, eventdata, handles)
% hObject    handle to radio_US (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proc_mix_region
global hdl para
para.proc.source=1;
hdl.load=handles;
% construct parameter filename
para_setting_filename=['proc_parameters_' para.survey_year];
eval(para_setting_filename);

return

% --- Executes on button press in radio_CAN.
function radio_CAN_Callback(hObject, eventdata, handles)
% hObject    handle to radio_CAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proc_mix_region
global hdl para 
para.proc.source=2;
hdl.load=handles;
% construct parameter filename
para_setting_filename=['proc_parameters_' para.survey_year];
eval(para_setting_filename);

return

% --- Executes on button press in radio_ALL.
function radio_ALL_Callback(hObject, eventdata, handles)
% hObject    handle to radio_ALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proc_mix_region
global hdl para

para.proc.source=3;
hdl.load=handles;
% construct parameter filename
para_setting_filename=['proc_parameters_' para.survey_year];
eval(para_setting_filename);


return



