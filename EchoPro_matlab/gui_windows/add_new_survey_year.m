function varargout = add_new_survey_year(varargin)
% ADD_NEW_SURVEY_YEAR MATLAB code for add_new_survey_year.fig
%      ADD_NEW_SURVEY_YEAR, by itself, creates a new ADD_NEW_SURVEY_YEAR or raises the existing
%      singleton*.
%
%      H = ADD_NEW_SURVEY_YEAR returns the handle to a new ADD_NEW_SURVEY_YEAR or the handle to
%      the existing singleton*.
%
%      ADD_NEW_SURVEY_YEAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADD_NEW_SURVEY_YEAR.M with the given input arguments.
%
%      ADD_NEW_SURVEY_YEAR('Property','Value',...) creates a new ADD_NEW_SURVEY_YEAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before add_new_survey_year_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to add_new_survey_year_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help add_new_survey_year

% Last Modified by GUIDE v2.5 29-Nov-2013 16:08:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @add_new_survey_year_OpeningFcn, ...
                   'gui_OutputFcn',  @add_new_survey_year_OutputFcn, ...
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


% --- Executes just before add_new_survey_year is made visible.
function add_new_survey_year_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to add_new_survey_year (see VARARGIN)

% Choose default command line output for add_new_survey_year
global hdl
handles.output = hObject;
hdl.add_new_survey_year=hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes add_new_survey_year wait for user response (see UIRESUME)
% uiwait(handles.add_new_survey_name);


% --- Outputs from this function are returned to the command line.
function varargout = add_new_survey_year_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edt_new_survey_year_Callback(hObject, eventdata, handles)
% hObject    handle to edt_new_survey_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_new_survey_year as text
%        str2double(get(hObject,'String')) returns contents of edt_new_survey_year as a double
global para hdl
tmp=char(get(handles.edt_new_survey_year,'string'));
available_survey_years=get(hdl.survey_year,'string');
available_survey_years(length(available_survey_years)+1)={tmp};
set(hdl.survey_year,'string',available_survey_years);
set(hdl.survey_year,'value',length(available_survey_years));
para.survey_year=tmp(1:4);
delete(handles.add_new_survey_name)
initialization(para.survey_year);

% --- Executes during object creation, after setting all properties.
function edt_new_survey_year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_new_survey_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global para hdl
    
tmp=char(get(handles.edt_new_survey_year,'string'));
available_survey_years=get(hdl.survey_year,'string');
available_survey_years(length(available_survey_years)+1)={tmp};
set(hdl.survey_year,'string',available_survey_years);
para.survey_year=tmp(1:4);
set(hdl.survey_year,'value',length(available_survey_years));
delete(handles.add_new_survey_name)
initialization(para.survey_year);
return
