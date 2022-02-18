function		load_para_file(filename,para_type)
%  Load parameters from a file
%% para_type = 1   		load variogram parameter file
%%		     = 2			load kriging parameter file
%% 			 = 3			load a file with both variogram and kriging parameters
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data

%% do not alter other parameters
para1=para;
data1=data;					

%% some load file parameters need to be preserved 

%% save particular parameters
HDIR=para.home_dir;
current_matlab_version=para.Matlab_Version;
load_vario_para1=para.vario.load_para;
vario_parafile=para.vario.para_file;

if para_type > 1
load_krig_loadpara=para.krig.load_para;                 % load parameter file flag
load_krig_variopara=para.krig.vario_para;               % load variogram parameter file flag
load_krig_krigpara=para.krig.krig_para;                 % load krig parameter file flag
load_krig_bothpara=para.krig.both_para;                 % load both variogram and krig file flag
krig_parafile_in=para.krig.para_file_in;                % parameter filename
krig_loaddata=para.krig.load_data_file;                 % load data file flag
krig_datafile=para.krig.data_file;                      % datafile name

end
%para1.dataprep.latlonfac=para.dataprep.latlonfac;

eval(['load '''  filename ''''])
switch para_type
	case 1									% variogram/correlogram parameters only
		para1.vario=para.vario;
	case 2									% kriging parameters only
		para1.krig=para.krig;
	case 3									% dataprep, variogram/correlogram, and kriging parameters
		para1.vario=para.vario;
		para1.krig=para.krig;
%		para1.dataprep=para.dataprep;
end
para=para1;
data=data1;

if para.krig.load_griddata_file == 1
    set(hdl.krig.customized_gridfile,'value',1);
    set(hdl.krig.gridfile_browser,'enable','on');
    krig_message(1,['Load Customer-specified Grid File:  ' para.krig.grid_file]);
end
%% recover particular parameters
para.Matlab_Version=current_matlab_version;
para.vario.load_para=load_vario_para1;
para.vario.para_file=vario_parafile;
para.home_dir=HDIR;
if para_type > 1
para.krig.load_para=load_krig_loadpara;
para.krig.vario_para=load_krig_variopara;
para.krig.krig_para=load_krig_krigpara;
para.krig.both_para=load_krig_bothpara;
para.krig.para_file_in=krig_parafile_in;
para.krig.load_data_file=krig_loaddata;
para.krig.data_file=krig_datafile;
end
return