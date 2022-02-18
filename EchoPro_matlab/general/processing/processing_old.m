function varargout = Processing(varargin)
% PROCESSING M-file for Processing.fig
%      PROCESSING, by itself, creates a new PROCESSING or raises the existing
%      singleton*.
%
%      H = PROCESSING returns the handle to a new PROCESSING or the handle to
%      the existing singleton*.
%
%      PROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSING.M with the given input arguments.
%
%      PROCESSING('Property','Value',...) creates a new PROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Processing_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Processing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Processing

% Last Modified by GUIDE v2.5 04-April-2013 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Processing_OpeningFcn, ...
                   'gui_OutputFcn',  @Processing_OutputFcn, ...
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


% --- Executes just before Processing is made visible.
function Processing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Processing (see VARARGIN)

global para
% Choose default command line output for Processing
handles.output = hObject;

if para.proc.stratification_index == 0
    % INPFC
    set(handles.popup_proc_stratification,'value', 2);
elseif para.proc.stratification_index == 1
    % K-S
    set(handles.popup_proc_stratification,'value', 1);
end
   
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Processing wait for user response (see UIRESUME)
% uiwait(handles.Processing);
% hdl.proc=handles;
return

% --- Outputs from this function are returned to the command line.
function varargout = Processing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return


% --- Executes on selection change in popup_proc_scale.
function popup_proc_scale_Callback(hObject, eventdata, handles)
% hObject    handle to popup_proc_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_proc_scale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_proc_scale
% disp('popup_proc_scale')
global hdl
hdl.proc=handles;
return

% --- Executes during object creation, after setting all properties.
function popup_proc_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_proc_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% global hdl 

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% hdl.proc=handles;

return

% --- Executes on selection change in popup_trawl_data.
function popup_trawl_data_Callback(hObject, eventdata, handles)
% hObject    handle to popup_trawl_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_trawl_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_trawl_data
% disp('popup_trawl_data')
global hdl

hdl.proc=handles;
%% trawl data:  1 = Acoustic Survey; 2 = Bottom Trawl Survey; 3 = Observer Data
para.bio_data_type = get(handles.popup_trawl_data,'value');         

return

% --- Executes during object creation, after setting all properties.
function popup_trawl_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_trawl_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global hdl 

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% hdl.proc=handles;

return

function edit_start_transect_num_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start_transect_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start_transect_num as text
%        str2double(get(hObject,'String')) returns contents of edit_start_transect_num as a double
% global hdl
% hdl.proc=handles;
para.proc.start_transect=str2num(get(hObject,'string'));
return

% --- Executes during object creation, after setting all properties.
function edit_start_transect_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start_transect_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% global hdl 

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% hdl.proc=handles;

return

function edit_end_transect_num_Callback(hObject, eventdata, handles)
% hObject    handle to edit_end_transect_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_end_transect_num as text
%        str2double(get(hObject,'String')) returns contents of edit_end_transect_num as a double
global hdl para
% hdl.proc=handles;
para.proc.end_transect=str2num(get(hObject,'string'));

return

% --- Executes during object creation, after setting all properties.
function edit_end_transect_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_end_transect_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% global hdl

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% hdl.proc=handles;
return

function edit_transect_reduction_fraction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_transect_reduction_fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_transect_reduction_fraction as text
%        str2double(get(hObject,'String')) returns contents of edit_transect_reduction_fraction as a double
global hdl para
hdl.proc=handles;
para.proc.transect_reduction_fraction=str2num(get(hObject,'string'));

return

function popup_transect_reduction_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start_transect_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global hdl para
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
para.proc.transect_reduction_mode=get(hObject,'value');
% hdl=handles;
return

% --- Executes on selection change in popup_transect_reduction_mode.
function popup_transect_reduction_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popup_transect_reduction_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_transect_reduction_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_transect_reduction_mode
global hdl para
hdl.proc=handles;
para.proc.transect_reduction_mode=get(hObject,'value');

return

% --- Executes on selection change in edit_bootstrap_limit.
function edit_bootstrap_limit_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bootstrap_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edit_bootstrap_limit contents as cell array

global hdl para
hdl.proc=handles;
para.proc.bootstrap.limit=str2num(get(hObject,'string'));

return

% --- Executes on button press in pb_start.
function pb_start_Callback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para
hdl.proc=handles;
para.proc.bootstrap.cnt=0;
proc_acoustic_data
return

% --- Executes on button press in pb_quit.
function pb_quit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hdl.proc=handles;
close
return


% --- Executes on button press in pb_proc_advanced.
function pb_proc_advanced_Callback(hObject, eventdata, handles)
% hObject    handle to pb_proc_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl para data
hdl.proc=handles;
hdl.proc.adv_proc=advanced_processing;   % popup window for advanced process
%% set correct  processing parameters
hdl_radio_excluding_age1=findobj(hdl.proc.adv_proc,'Tag','radio_excluding_age1');
set(hdl_radio_excluding_age1,'value',para.proc.exclude_age1);
hdl_popupmenu_len_wgt_method=findobj(hdl.proc.adv_proc,'Tag','popupmenu_len_wgt_method');
set(hdl_popupmenu_len_wgt_method,'value',para.proc.len_wgt_method);

return


% --- Executes on selection change in popup_processing_stratification.
function popup_proc_stratification_Callback(hObject, eventdata, handles)
% hObject    handle to popup_processing_stratification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_processing_stratification contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_processing_stratification
global hdl para 
hdl.proc=handles;
if get(hObject,'value') == 1
    %% post-stratification: KS analysis, or trawl-based stratification
    para.proc.KS_stratification=1;
    para.proc.stratification_index=1;
else
    %% geographically defined stratification
    para.proc.KS_stratification=0;
    para.proc.stratification_index=0;
end
return

% --- Executes during object creation, after setting all properties.
function popup_proc_stratification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_processing_stratification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global hdl para 

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% if para.proc.stratification_index > 1
if para.proc.stratification_index ~= 1
   para.proc.KS_stratification=0;
   set(hObject,'value',2);
end
% hdl.proc=handles;
return

% --- Executes on selection change in popup_proc_variable.
function popup_proc_variable_Callback(hObject, eventdata, handles)
% hObject    handle to popup_proc_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_proc_variable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_proc_variable
global hdl para
hdl.proc=handles;
para.proc.kriging_input=get(hObject,'value');   % 1 = biomass density; 2 = NASC;  3 = number density
return

% --- Executes on button press in radio_proc_mix_region.
function radio_proc_mix_region_Callback(hObject, eventdata, handles)
% hObject    handle to radio_proc_mix_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proc_mix_region
global hdl
hdl.proc=handles;
return

% --- Executes on button press in radio_proc_include_age1.
function radio_proc_include_age1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_proc_include_age1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_proc_include_age1
global hdl para
para.proc.exclude_age1 = 1 - get(hObject, 'value');          % 0 = include age 1 hake,   1 = exclude age 1 hake
if para.proc.exclude_age1 == 1
    para.acoust.filename.processed_data=para.acoust.filename.processed_data_age2;
    para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age2;
else
    para.acoust.filename.processed_data=para.acoust.filename.processed_data_age1;
    para.bio_acoust.filename.Transect_region_haul=para.bio_acoust.filename.Transect_region_haul_age1;
end

hdl.proc=handles;
return

% --- Executes on button press in radio_extrapolation.
function radio_extrapolation_Callback(hObject, eventdata, handles)
% hObject    handle to radio_extrapolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_extrapolation
global hdl para
para.proc.extrapolation=1;   % With extrapolation for kriging
hdl.proc=handles;
para_setting_filename=['proc_parameters_' para.survey_year];
eval(para_setting_filename);

return

% --- Executes on button press in radio_no_extrapolation.
function radio_no_extrapolation_Callback(hObject, eventdata, handles)
% hObject    handle to radio_no_extrapolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_no_extrapolation
global hdl para
para.proc.extrapolation=0;   % No extrapolation for kriging
hdl.proc=handles;
para_setting_filename=['proc_parameters_' para.survey_year];
eval(para_setting_filename);
return

% % --- Executes on button press in radio_extrapolation.
% function radio_extrapolation_Callback(hObject, eventdata, handles)
% % hObject    handle to radio_extrapolation (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of radio_extrapolation
% global hdl para
% para.proc.extrapolation=1;   % With extrapolation for kriging
% hdl.proc=handles;
% para_setting_filename=['proc_parameters_' para.survey_year];
% eval(para_setting_filename);
% 
% return

% --- Executes on button press in radio_agedata.
function radio_agedata_Callback(hObject, eventdata, handles)
% hObject    handle to radio_agedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_agedata

global hdl para

hdl.proc=handles;
if   para.proc.age_data_status == 2
    set(handles.radio_agedata,'Value',1);
elseif para.proc.age_data_status == 1
    set(handles.radio_agedata,'Value',0);
end
para.proc.age_data_status=2-get(handles.radio_agedata,'Value');
return

% --- Executes on button press in radio_kriging_proc.
function radio_kriging_proc_Callback(hObject, eventdata, handles)
% hObject    handle to radio_kriging_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_kriging_proc

global hdl para

para.proc.kriging=get(handles.radio_kriging_proc,'value');
hdl.proc=handles;
return

% --- Executes on button press in radio_default_para
function radio_default_para_Callback(hObject, eventdata, handles)
% hObject    handle to radio_default_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_kriging_proc

global hdl para

para.proc.default_para=get(handles.radio_default_para,'value');
hdl.proc=handles;
return



