function varargout = advanced_processing(varargin)
% ADVANCED_PROCESSING M-file for advanced_processing.fig
%      ADVANCED_PROCESSING, by itself, creates a new ADVANCED_PROCESSING or raises the existing
%      singleton*.
%
%      H = ADVANCED_PROCESSING returns the handle to a new ADVANCED_PROCESSING or the handle to
%      the existing singleton*.
%
%      ADVANCED_PROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCED_PROCESSING.M with the given input arguments.
%
%      ADVANCED_PROCESSING('Property','Value',...) creates a new ADVANCED_PROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before advanced_processing_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to advanced_processing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help advanced_processing

% Last Modified by GUIDE v2.5 26-Oct-2011 08:53:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @advanced_processing_OpeningFcn, ...
                   'gui_OutputFcn',  @advanced_processing_OutputFcn, ...
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


% --- Executes just before advanced_processing is made visible.
function advanced_processing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to advanced_processing (see VARARGIN)

% Choose default command line output for advanced_processing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes advanced_processing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = advanced_processing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return


% --- Executes on button press in pb_browse_input_filepath.
function pb_browse_input_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to pb_browse_input_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para
hdl.proc.adv_proc=handles;

filepath=uigetdir(para.home_dir,'Pick a Directory');
if filepath ~= 0
    set(hObject,'UserData',filepath);
    para.acoust.filename.EV_export_input_filepath=filepath;
    set(hdl.proc.adv_proc.text_input_file_path,'string',para.acoust.filename.EV_export_input_filepath);
end
return

% --- Executes on button press in pb_browse_mixed_filepath.
function pb_browse_mixed_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to pb_browse_mixed_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para
hdl.proc.adv_proc=handles;

filepath=uigetdir(para.home_dir,'Pick a Directory');
if filepath ~= 0
    set(hObject,'UserData',filepath);
    para.acoust.filename.EV_export_mixed_filepath=filepath;
    set(hdl.proc.adv_proc.text_mixed_file_path,'string',para.acoust.filename.EV_export_mixed_filepath);
end
return

% --- Executes on button press in pb_browse_output.
function pb_browse_output_Callback(hObject, eventdata, handles)
% hObject    handle to pb_browse_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl para
hdl.proc.adv_proc=handles;

[filename,filepath]=uiputfile({'*.xlsx','xlsx files (*.xlsx)';'*.*',' ALL (*.*)'},'Specify an output NASC file');
if filename ~= 0
    set(hObject,'UserData',[filepath filename]);
    para.acoust.filename.proc_EV_NASC=[filepath filename];
    set(hdl.proc.adv_proc.text_output_filename,'string',para.acoust.filename.proc_EV_NASC);
end
return

% --- Executes on button press in pb_apply.
function pb_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pb_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.proc.adv_proc=handles;
generate_EV_EI_NASC_table
return


% --- Executes on button press in radio_output_filename.
function radio_output_filename_Callback(hObject, eventdata, handles)
% hObject    handle to radio_output_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_output_filename

global hdl
hdl.proc.adv_proc=handles;
return

% --- Executes on button press in radio_1yp.
function radio_1yp_Callback(hObject, eventdata, handles)
% hObject    handle to radio_1yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_1yp

global hdl
hdl.proc.adv_proc=handles;
return


% --- Executes on button press in radio_2yp.
function radio_2yp_Callback(hObject, eventdata, handles)
% hObject    handle to radio_2yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_2yp

global hdl
hdl.proc.adv_proc=handles;
return

% --- Executes on button press in radio_1y_only.
function radio_1y_only_Callback(hObject, eventdata, handles)
% hObject    handle to radio_1y_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_1y_only

global hdl
hdl.proc.adv_proc=handles;
return

% --- Executes on button press in radio_mixed.
function radio_mixed_Callback(hObject, eventdata, handles)
% hObject    handle to radio_mixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_mixed

global hdl
hdl.proc.adv_proc=handles;
return


% --- Executes on button press in radio_excluding_age1.
function radio_excluding_age1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_excluding_Age1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_excluding_Age1

global hdl para
hdl.proc.adv_proc=handles;
para.proc.exclude_age1=1-get(hObject,'Value');
%% changed on 11-04-2021 to correct the radio handle
set(hdl.proc.adv_proc.radio_exclude_age1,'Value',para.proc.exclude_age1);
return

% --- Executes on button press in radio_ordered_by_transects.
function radio_ordered_by_transects_Callback(hObject, eventdata, handles)
% hObject    handle to radio_ordered_by_transects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_ordered_by_transects

global hdl para
hdl.proc.adv_proc=handles;
para.proc.ordered_by_transects=get(hObject,'Value');
return


% --- Executes on button press in popupmenu_len_wgt_method.
function popupmenu_len_wgt_method_Callback(hObject, eventdata, handles)
% hObject    handle to radio_excluding_Age1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of popupmenu_len_wgt_method

global hdl para
hdl.proc.adv_proc=handles;
para.proc.len_wgt_method=get(hObject,'Value');
return

