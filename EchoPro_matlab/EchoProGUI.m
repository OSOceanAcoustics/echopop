% Program EchoProGUI processes the acoustic (Exported from EchoView or raw EK60) and biological
% trawl sampling (FSCS) data to obtain length-age-sex structured hake (other species) biomass and
% variance estimates, and generate a set of output files (in MS Office spreadsheet format, i.e. .xlsx files)
%
% Written by Dezhang Chu, NOAA Fisheries, NWFSC
% last revision 3/23/2013
% last revision 3/23/2020 --> add SD option

% clear all; 
clear
close all; 
clear global;
clear persistent;

global data para hdl        %% global variables

home_dir=pwd;               %% set home directory
ind=find(home_dir == ':');
para.drive=home_dir(1:ind);
machine_bits=64;            %% bits of the computer


para.home_dir=home_dir;
addpath(genpath(home_dir))
% para.data_root_dir='C:\Projects\EchoPro\Historical Summary (for Kriging)\';
% para.data_root_dir='F:\Historical Summary (for Kriging)\';
% para.data_root_dir='F:\Backup Laptop 8283 -2016-09-23\Projects\EchoPro\Historical Summary (for Kriging)\';
% para.data_root_dir='N:\Survey.Acoustics\Survey Time Series Analysis\Historical Summary (for Kriging)\';
para.data_root_dir='/usr/mayorgadat/workmain/acoustics/2021-NWFSC-EchoPro-HakeIGP/EchoPro/EchoProGUI_reorg/inputs/';

if machine_bits == 64
    cmd0=['addpath(genpath(''' home_dir '\gui_windows64''))'];
else
    cmd0=['addpath(genpath(''' home_dir '\gui_windows32''))'];
end
eval(cmd0)

hdl.main=MainWindow;           %% Main GUI window (top or the 1st level GUI window)

survey_year_indx=get(findobj(hdl.main,'Tag','popup_survey_year'),'string');
para.survey_year=char(survey_year_indx(get(findobj(hdl.main,'Tag','popup_survey_year'),'value')));
initialization(para.survey_year)          %% initialization

return


