function varargout = show_parameters(varargin)
% SHOW_PARAMETERS M-file for show_parameters.fig
%      SHOW_PARAMETERS, by itself, creates a new SHOW_PARAMETERS or raises the existing
%      singleton*.
%
%      H = SHOW_PARAMETERS returns the handle to a new SHOW_PARAMETERS or the handle to
%      the existing singleton*.
%
%      SHOW_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_PARAMETERS.M with the given input arguments.
%
%      SHOW_PARAMETERS('Property','Value',...) creates a new SHOW_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_parameters

% Last Modified by GUIDE v2.5 04-Jul-2011 03:36:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @show_parameters_OutputFcn, ...
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


% --- Executes just before show_parameters is made visible.
function show_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_parameters (see VARARGIN)

% Choose default command line output for show_parameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_acoustical.
function pb_acoustical_Callback(hObject, eventdata, handles)
% hObject    handle to pb_acoustical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para
hdl.show_para=handles;
disp_str=get_disp_str(hdl.show_para.text_parameters,para,3);
return

% --- Executes on button press in pb_general.
function pb_general_Callback(hObject, eventdata, handles)
% hObject    handle to pb_general (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hdl data para
hdl.show_para=handles;
disp_str=get_disp_str(hdl.show_para.text_parameters,para,1);

return

% --- Executes on button press in pb_biological.
function pb_biological_Callback(hObject, eventdata, handles)
% hObject    handle to pb_biological (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl data para
hdl.show_para=handles;
disp_str=get_disp_str(hdl.show_para.text_parameters,para,2);
return
