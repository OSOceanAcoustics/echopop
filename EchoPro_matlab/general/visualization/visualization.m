function varargout = visualization(varargin)
% visualization M-file for visualization.fig
%      VISUALIZATION, by itself, creates a new VISUALIZATION or raises the existing
%      singleton*.
%
%      H = VISUALIZATION returns the handle to a new VISUALIZATION or the handle to
%      the existing singleton*.
%
%      VISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZATION.M with the given input arguments.
%
%      VISUALIZATION('Property','Value',...) creates a new VISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualization_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visualization

% Last Modified by GUIDE v2.5 01-Jul-2011 19:51:27

% Begin initialization code - DO NOT EDIT
%global hdl para

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualization_OpeningFcn, ...
                   'gui_OutputFcn',  @visualization_OutputFcn, ...
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


% --- Executes just before visualization is made visible.
function visualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualization (see VARARGIN)

global hdl

% Choose default command line output for visualization
handles.output = hObject;
hdl.vis=handles;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes visualization wait for user response (see UIRESUME)
% uiwait(handles.visualization);
hObject=hdl.vis.radio_visual_kriged;
set(handles.radio_visual_func_length,'value',0);
set(handles.radio_visual_func_biomass,'value',1);
set(handles.radio_visual_var_lat_lon,'value',1);
data_type_selection(hObject, eventdata, handles)
return


% --- Outputs from this function are returned to the command line.
function varargout = visualization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global hdl

set(hdl.vis.radio_visual_un_kriged,'value',1)
varargout{1} = handles.output;
hdl.vis=handles;
return

% --- Executes on button press in pbStart.
function pbStart_Callback(hObject, eventdata, handles)
% hObject    handle to pbStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.vis=handles;
proc_visualization;
return

% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.vis.school=handles;
save_visualization_results;
return

% --- Executes on button press in pbQuit.
function pbQuit_Callback(hObject, eventdata, handles)
% hObject    handle to pbQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global hdl
%hdl.vis=handles;
close
return

% --- Executes on button press in pb_fish_school.
function pb_fish_school_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fish_school (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl para
if str2num(para.survey_year) >= 2003
    hdl.vis=handles;
    hdl.fish_school=fish_school;
else
    errordlg('No adequate information for survey years pripor to 2003!!', 'Error Message')
end
return

% --- Executes on button press in radio_visual_un_kriged.
function radio_visual_un_kriged_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_un_kriged (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para 
hdl.vis=handles;
data_type_selection(hObject, eventdata, handles);
return


% --- Executes on button press in radio_visual_kriged.
function radio_visual_kriged_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_kriged (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para 
hdl.vis=handles;
data_type_selection(hObject, eventdata, handles);
return

% --- Executes on button press in radio_visual_biological.
function radio_visual_biological_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_biological (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para 
hdl.vis=handles;
data_type_selection(hObject, eventdata, handles);
return

% --- Executes on button press in radio_visual_func_length.
function radio_visual_func_length_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_func_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
global hdl para
hdl.vis=handles;
set_enabled_selection(hObject, eventdata, handles );
return

% --- Executes on button press in radio_visual_func_age.
function radio_visual_func_age_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_func_age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para 
hdl.vis=handles;
set_enabled_selection(hObject, eventdata, handles );
return


% --- Executes on button press in radio_visual_func_NASC.
function radio_visual_func_NASC_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_func_NASC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para
hdl.vis=handles;
set_enabled_selection(hObject, eventdata, handles );
return

% --- Executes on button press in radio_visual_func_number.
function radio_visual_func_number_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_func_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para
hdl.vis=handles;
set_enabled_selection(hObject, eventdata, handles );
return


% --- Executes on button press in radio_visual_func_biomass.
function radio_visual_func_biomass_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_func_biomass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para
hdl.vis=handles;
set_enabled_selection(hObject, eventdata, handles );
return


% --- Executes on button press in radio_visual_var_transect.
function radio_visual_var_transect_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_transect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl para data
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_lat.
function radio_visual_var_lat_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_strata.
function radio_visual_var_strata_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_strata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_trawl.
function radio_visual_var_trawl_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_trawl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_length.
function radio_visual_var_length_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_age.
function radio_visual_var_age_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return


% --- Executes on button press in radio_visual_var_gender.
function radio_visual_var_gender_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_gender (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_weight.
function radio_visual_var_weight_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return


% --- Executes on button press in radio_visual_layer_bot_depth.
function radio_visual_var_layer_depth_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_layer_bot_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_bot_depth.
function radio_visual_var_bot_depth_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_bot_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

% --- Executes on button press in radio_visual_var_survey_region.
function radio_visual_var_survey_region_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_len_survey_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return


% --- Executes on button press in radio_visual_var_lat_lon.
function radio_visual_var_lat_lon_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_var_lat_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return


% --- Executes on button press in radio_visual_var_len_age.
function radio_visual_var_len_age_Callback(hObject, eventdata, handles)
% hObject    handle to radio_visual_hist_len_age (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl 
hdl.vis=handles;
return

