function varargout = MainWindow(varargin)
% MAINWINDOW M-file for MainWindow.fig
%      MAINWINDOW, by itself, creates a new MAINWINDOW or raises the existing
%      singleton*.
%
%      H = MAINWINDOW returns the handle to a new MAINWINDOW or the handle to
%      the existing singleton*.
%
%      MAINWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINWINDOW.M with the given input arguments.
%
%      MAINWINDOW('Property','Value',...) creates a new MAINWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainWindow_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainWindow

% Last Modified by GUIDE v2.5 04-Jul-2011 14:13:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @MainWindow_OutputFcn, ...
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


% --- Executes just before MainWindow is made visible.
function MainWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainWindow (see VARARGIN)

% Choose default command line output for MainWindow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainWindow wait for user response (see UIRESUME)
% uiwait(handles.MainWindow);
global hdl
hdl.main=handles;


% --- Outputs from this function are returned to the command line.
function varargout = MainWindow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global hdl
varargout{1} = handles.output;
hdl.main=handles;

% --- Executes on popupmenu in survey_year.
function popup_survey_year_Callback(hObject, eventdata, handles)
% hObject    handle to popup_survey_year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para hdl
if  get(handles.popup_survey_year,'value') == 1
    hdl.survey_year=hObject;
    add_new_survey_year
else
    available_survey_years=get(hObject,'string');
    selected_year_ind=get(hObject,'value');
    tmp=char(available_survey_years(selected_year_ind));
    para.survey_year=tmp(1:4);
%     initialization(para.survey_year);
    switch para.survey_year
        case {'1995','2005','2007'}
            para.proc.source=1;                 % 1 = US only
        case '2003'
            para.proc.source=2;                 % 2 = CAN only            
        otherwise
            para.proc.source=3;                 % 3 = Both US and CAN            
    end    
end
%% Last Modification:  3/25/2020  --> 2008 and 2010 IVC not available
if strcmp(para.survey_year, '2008') == 1 | strcmp(para.survey_year, '2010') == 1
    f = warndlg(sprintf('Selection of ''%s'' is currently not avialable!',tmp), 'Warning Message');
else
    load_proc_parameters(para.survey_year);         % popup datafile names
end

return

% --- Executes on popupmenu in survey_region.
function popup_survey_region_Callback(hObject, eventdata, handles)
% hObject    handle to popup_survey_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para hdl

% Survey Region: 1 = US; 2 CAN; 3 = all
survey_region_index = get(handles.popup_survey_region,'value');
if survey_region_index <= 3
    para.proc.source = get(handles.popup_survey_region,'value');
    para.platform_name = 'FSV';
else
    para.proc.source = get(handles.popup_survey_region,'value') - 3;
    para.platform_name = 'SD';
    para.bio_data_type = 3;
end
load_proc_parameters(para.survey_year);

return

% --- Executes on button press in pbQuit.
function pbQuit_Callback(hObject, eventdata, handles)
% hObject    handle to pbQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para

%% close all windows 
if isfield(hdl,'redefine_filenames') & ishandle(hdl.redefine_filenames)
    delete(hdl.redefine_filenames)
end
if isfield(hdl,'show_parameter') & ishandle(hdl.show_parameter)
    delete(hdl.show_parameter)
end
if isfield(hdl,'processing') & ishandle(hdl.processing)
    delete(hdl.processing)
end
if isfield(hdl,'visualization') & ishandle(hdl.visualization)
    delete(hdl.visualization)
end
if isfield(hdl,'report') & ishandle(hdl.report)
    delete(hdl.report)
end
if isfield(hdl,'fish_school') & ishandle(hdl.fish_school)
    delete(hdl.fish_school)
end

close all
return

% --- Executes on button press in pb_refresh.
function pb_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to pb_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global para

survey_year_ind=get(handles.popup_survey_year,'value');
survey_years=get(handles.popup_survey_year,'string');
initialization
para.survey_year=char(survey_years(survey_year_ind));
% load_proc_parameters(para.survey_year); % popup datafile names
return


% --- Executes on button press in pb_reload_para.
function pb_redefine_filenames_Callback(hObject, eventdata, handles)
% hObject    handle to pb_reload_para (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para
para.tasks.opr_indx=1;
get_operation_window(hObject,handles);

% --- Executes on button press in pb_show_parameter.
function pb_show_parameter_Callback(hObject, eventdata, handles)
% hObject    handle to pb_show_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para

para.tasks.opr_indx=2;
get_operation_window(hObject,handles);


% --- Executes on button press in pb_processing.
function pb_processing_Callback(hObject, eventdata, handles)
% hObject    handle to pb_processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para
para.tasks.opr_indx=3;
get_operation_window(hObject,handles);


% --- Executes on button press in pb_visualization.
function pb_visualization_Callback(hObject, eventdata, handles)
% hObject    handle to pb_visualization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para
para.tasks.opr_indx==4;
get_operation_window(hObject,handles);

% --- Executes on button press in pb_report.
function pb_report_Callback(hObject, eventdata, handles)
% hObject    handle to pb_report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global para
para.tasks.opr_indx=5;
get_operation_window(hObject,handles);




