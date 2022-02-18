function	initialization2011x();
% function initialization perform initialization

global para color data hdl

    
pp=path;
%% Optimization Toolbox
OptimOption=0;
if strfind(pp,'optim')
    OptimOption=1;
end
HDIR=para.home_dir;
PlatForm=1;                 % Window OS
markersize=6;               % for Window OS

% set parameters
para.krig_home_dir=HDIR;
para.platform=PlatForm;
para.optim=OptimOption;
%para.status=0;
para.file_dir.data_conversion=HDIR;
para.file_dir.datafile=HDIR;
para.file_dir.data_format_file=HDIR;
para.file_dir.gridfile=HDIR;
para.file_dir.parafile=HDIR;
para.file_dir.batch_filename=HDIR;
para.file_dir.batch_log=HDIR;
para.file_dir.mat_file_in=HDIR;
para.file_dir.mat_file_out=HDIR;

para.dataprep.filename='';
para.dataprep.ext_prog=0;
para.dataprep.dat_conv_fname='';
para.dataprep.xy_switch=0;
para.krig.data_format_file=[];
para.status.dataprepfig=0;
para.status.dataprep=0;
para.status.variogramfig=0;
para.status.variogram=0;
para.status.krigingfig=0;
para.status.kriging=0;
para.status.dispkrigfig=0;
para.status.dispkrig=0;

hdl.status.dataprepfig=0;
hdl.status.variogramfig=0;
hdl.status.krigingfig=0;
hdl.status.dispkrigfig=0;

%% set default parameters
para.dataprep.ytox=1;
para.dataprep.ztox=1;
para.dataprep.ext_prog=0;
para.dataprep.filter_type=2;				% default filter = mean
para.dataprep.reduct_fac=1;
para.dataprep.filter_supt=1;
para.krig.load_data_format_file=0;
para.dataprep.data_disptype=1;
para.dataprep.data_disptype0=1;

para.vario.max_nugt=1;
para.vario.max_sill=1.5;
para.vario.max_powr=4.0;
para.vario.max_range=sqrt(2);			% normalized range
para.vario.max_hole=4*pi/para.vario.max_range;
para.vario.max_lscl=para.vario.max_range;
para.vario.res0=0.002;
para.vario.range0=0.06;
para.vario.disp_range0=0.06;
para.vario.load_para=0;
para.vario.para_file='';
para.krig.load_para=0;
para.krig.vario_para=0;
para.krig.krig_para=0;
para.krig.both_para=1;
para.krig.para_file_in='';
para.krig.load_data_file=0;
para.krig.batch_file_proc=0;
para.krig.bat_proc_cnt=0;
para.krig.data_file='';

para.dispkrig.markersize=markersize;
para.dispkrig.customized_grid_data_markersize=4;
para.krig.load_griddata_file=0;
para.dispkrig.trackline.dispflag=1;

% color
color.background=[0.8 0.8 0.8];
color.grey=[0.75 0.75 0.75];
color.dark_grey=[0.65 0.65 0.65];
color.blue=[0.753 0.753 0.753];

% data structure
data.in.dim=2;          % default 2D data

% graphic handle
hdl.object.edit_w=0.04;
hdl.object.pushbtn_w=0.05;
hdl.object.pushbtn_l=0.14;
hdl.object.popmenu=0.04;
hdl.object.radio_w=0.03;
hdl.object.txt_w8=0.03;
hdl.object.txt_w10=0.05;

if para.platform == 2   % If Unix OS
   getXvisual;
end
hdl.msg.h0=[];

warning off

