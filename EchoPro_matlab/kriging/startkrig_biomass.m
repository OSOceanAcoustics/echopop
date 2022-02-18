%% The EchoPro Kriging Software Package is developed by
%% Dezhang Chu
%% Last modification: 3/20/2013


function    [data, para, hdl]=startkrig_biomass;


global para	 		% all setting and processing parameters
global data  		% input and output data
global hdl			% handles for all windows
global color		% color RGB

%% Matlab Version
Ver=ver;
Version=[];
for i = 1:length(Ver)
    p=strcmp(Ver(i).Name,'MATLAB');
    q=strcmp(Ver(i).Name,'MATLAB Toolbox');
    if p == 1 | q == 1
        Version = floor(str2num(Ver(i).Version(1:3)));
        break
    end
end
if isempty(Version)
    disp(['Can''t determine what version of your Matlab software'])
    return
end

para.Matlab_Version=Version;
krig_initialization;

para.proc.default_para=get(hdl.proc.radio_default_para,'value');
if para.proc.default_para ~= 1
   dataprep3dfig;
else
   non_GUI_kriging_process
end

return