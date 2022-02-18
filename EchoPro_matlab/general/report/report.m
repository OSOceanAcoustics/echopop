function varargout = report(varargin)
% REPORT M-file for report.fig
%      REPORT, by itself, creates a new REPORT or raises the existing
%      singleton*.
%
%      H = REPORT returns the handle to a new REPORT or the handle to
%      the existing singleton*.
%
%      REPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REPORT.M with the given input arguments.
%
%      REPORT('Property','Value',...) creates a new REPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before report_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to report_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help report

% Last Modified by GUIDE v2.5 04-Jul-2011 03:06:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @report_OpeningFcn, ...
                   'gui_OutputFcn',  @report_OutputFcn, ...
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


% --- Executes just before report is made visible.
function report_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to report (see VARARGIN)

% Choose default command line output for report
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes report wait for user response (see UIRESUME)
% uiwait(handles.Reports);


% --- Outputs from this function are returned to the command line.
function varargout = report_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radio_assessment.
function radio_assessment_Callback(hObject, eventdata, handles)
% hObject    handle to radio_assessment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_assessment
global hdl
hdl.proc=handles;
return



% --- Executes on button press in radio_oracle.
function radio_oracle_Callback(hObject, eventdata, handles)
% hObject    handle to radio_oracle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_oracle
global hdl
hdl.proc=handles;
return

% --- Executes on button press in pb_start.
function pb_output_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl para
hdl.proc.report=handles;
para.proc.output_filepath=uigetdir('Specify an Output Filepath!');
set(handles.text_table_output_filepath,'string',para.proc.output_filepath);
return


% --- Executes on button press in pb_start.
function pb_start_Callback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.proc=handles;
generate_reports(handles);
return
