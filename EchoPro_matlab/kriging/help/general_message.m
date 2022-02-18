function varargout = general_message(varargin)
% GENERAL_MESSAGE M-file for general_message.fig
%      GENERAL_MESSAGE, by itself, creates a new GENERAL_MESSAGE or raises the existing
%      singleton*.
%
%      H = GENERAL_MESSAGE returns the handle to a new GENERAL_MESSAGE or the handle to
%      the existing singleton*.
%
%      GENERAL_MESSAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENERAL_MESSAGE.M with the given input arguments.
%
%      GENERAL_MESSAGE('Property','Value',...) creates a new GENERAL_MESSAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before general_message_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to general_message_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help general_message

% Last Modified by GUIDE v2.5 04-May-2004 17:12:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @general_message_OpeningFcn, ...
                   'gui_OutputFcn',  @general_message_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before general_message is made visible.
function general_message_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to general_message (see VARARGIN)

% Choose default command line output for general_message
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

htxt=findobj(hObject,'Style','text');
set(htxt,'string',varargin{2});
hlistbox=findobj(hObject,'Style','listbox');
hpushbutton=findobj(hObject,'Style','pushbutton');
set(hlistbox,'string',varargin{3});
%% determine the sizes of the text and listbox objects, as well as the
%% window size based on the input contents
  v0=get(hObject,'position');
  v1=get(hpushbutton,'position');
  v2=get(hlistbox,'position');
  v3=get(htxt,'position');
  if iscell(varargin{2})
     n3=length(varargin{2});
  else
     n3=size(varargin{2},1);
  end
  if iscell(varargin{3})
     n2=length(varargin{3});
  else
     n2=size(varargin{3},1);
  end
  n2=min(n2,15);
  v2(4)=n2+5;
  v3(4)=n3+0.5;
  v2(2)=v1(2)+v1(4)+2.0;
  v3(2)=v2(2)+v2(4)+0.5;
  n0=v3(2)+v3(4);
  v0(4)=n0+3;
%  [v1;v2;v3;v0]
  set(hObject,'Position',v0,'units','characters');
  set(htxt,'Position',v3,'units','characters');
  set(hlistbox,'Position',v2,'units','characters');

% UIWAIT makes general_message wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = general_message_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


