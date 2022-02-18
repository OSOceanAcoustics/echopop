function varargout = fish_school(varargin)
% FISH_SCHOOL M-file for fish_school.fig
%      FISH_SCHOOL, by itself, creates a new FISH_SCHOOL or raises the existing
%      singleton*.
%
%      H = FISH_SCHOOL returns the handle to a new FISH_SCHOOL or the handle to
%      the existing singleton*.
%
%      FISH_SCHOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FISH_SCHOOL.M with the given input arguments.
%
%      FISH_SCHOOL('Property','Value',...) creates a new FISH_SCHOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fish_school_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fish_school_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fish_school

% Last Modified by GUIDE v2.5 27-Aug-2011 19:44:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fish_school_OpeningFcn, ...
                   'gui_OutputFcn',  @fish_school_OutputFcn, ...
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


% --- Executes just before fish_school is made visible.
function fish_school_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fish_school (see VARARGIN)

% Choose default command line output for fish_school
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fish_school wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fish_school_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbStart.
function pbStart_Callback(hObject, eventdata, handles)
% hObject    handle to pbStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.vis.school=handles;
proc_fish_school_visualization;
return

% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdl
hdl.vis.school=handles;
save_fish_school_results;
return

% --- Executes on button press in pbQuit.
function pbQuit_Callback(hObject, eventdata, handles)
% hObject    handle to pbQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data
if isfield(data.visual,'school')
    if isfield(data.visual.school,'func')
        data.visual.school=rmfield(data.visual.school,'func');
        data.visual.school=rmfield(data.visual.school,'func_str');
    end
    if isfield(data.visual.school,'var')
        data.visual.school=rmfield(data.visual.school,'var');
        data.visual.school=rmfield(data.visual.school,'var_str');
    end
    data.visual=rmfield(data.visual,'school');
end
close
return

% --- Executes on button press in radio_school_length.
function radio_school_length_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_length
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,6);
data.visual.school.func_str='School Length (nm)';
return

% --- Executes on button press in radio_school_height.
function radio_school_height_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_height
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,7);
data.visual.school.func_str='School Height (m)';
return

% --- Executes on button press in radio_school_area.
function radio_school_area_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_area
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,8)*1e-6;
data.visual.school.func_str='School Area (km^2)';
return

% --- Executes on button press in radio_school_NASC.
function radio_school_NASC_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_NASC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_NASC
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,9);
data.visual.school.func_str='NASC (m^2/nmi^2)';
return

% --- Executes on button press in radio_school_density.
function radio_school_density_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_density
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,10);
data.visual.school.func_str='Fish Density (no./m^3)';
return

% --- Executes on button press in radio_biomass.
function radio_school_biomass_Callback(hObject, eventdata, handles)
% hObject    handle to radio_biomass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_biomass
global hdl para data
hdl.vis.school=handles;
data.visual.school.func=data.final.table.school(:,11)*1e-6;
data.visual.school.func_str='Biomass (kmt)';
return


% --- Executes on button press in radio_transect.
function radio_transect_Callback(hObject, eventdata, handles)
% hObject    handle to radio_transect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_length
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=data.final.table.school(:,1);
data.visual.school.var_str='Transect';
return

% --- Executes on button press in radio_lat.
function radio_lat_Callback(hObject, eventdata, handles)
% hObject    handle to radio_lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_length
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=data.final.table.school(:,3);
data.visual.school.var_str='Latitude';
return

% --- Executes on button press in radio_region_seq_no.
function radio_region_seq_no_Callback(hObject, eventdata, handles)
% hObject    handle to radio_region_seq_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_length
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=data.final.table.school(:,2);
data.visual.school.var_str='Region Sequence No.';
return

% --- Executes on button press in radio_school_depth.
function radio_school_depth_Callback(hObject, eventdata, handles)
% hObject    handle to radio_school_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_length
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=data.final.table.school(:,5);
data.visual.school.var_str='Depth (m)';
return


% --- Executes on button press in radio_hist.
function radio_school_hist_Callback(hObject, eventdata, handles)
% hObject    handle to radio_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_height
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=[];
data.visual.school.var_str='Counts';
data.visual.school.histogram=1;
return


% --- Executes on button press in radio_PDF.
function radio_school_PDF_Callback(hObject, eventdata, handles)
% hObject    handle to radio_PDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_school_height
global hdl para data
hdl.vis.school=handles;
data.visual.school.var=[];
data.visual.school.var_str='PDF';
return



