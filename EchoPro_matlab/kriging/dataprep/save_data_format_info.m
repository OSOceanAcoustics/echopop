function   save_data_format_info(opt)
% function save_data_format_info.m saves the data file information and the
% associated preliminary processing parameters
% opt = 1   need to give a file name
%       2   use existing filename
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

para1.dataprep.reduct_fac=str2num(get(hdl.dataprep.reduct_fac,'string'));
para1.dataprep.filter_supt=str2num(get(hdl.dataprep.filter_supt,'string'));
para1.dataprep.filter_type=get(hdl.dataprep.filter,'value');				% filter type: 1-simple reduction,
para1.dataprep.var1_indx=get(hdl.dataprep.var1,'value');
para1.dataprep.var2_indx=get(hdl.dataprep.var2,'value');
para1.dataprep.x_axis_indx=get(hdl.dataprep.x_axis,'value');
para1.dataprep.y_axis_indx=get(hdl.dataprep.y_axis,'value');
if isfield(para.dataprep,'ext_prog');
   para1.dataprep.ext_prog=para.dataprep.ext_prog;
   para1.dataprep.dat_conv_fname=para.dataprep.dat_conv_fname;
   para1.file_dir.data_conversion=para.file_dir.data_conversion;
end
para1.dataprep.transform_index=para.dataprep.transform_index;
para1.dataprep.ytox=para.dataprep.ytox;
para1.dataprep.ztox=para.dataprep.ztox;
para1.dataprep.data_dim=2;
para1.dataprep.xlabel=para.dataprep.xlabel;
para1.dataprep.ylabel=para.dataprep.ylabel;
if data.in.dim == 3
   para1.dataprep.z_axis_indx=get(hdl.dataprep.z_axis,'value');
   para1.dataprep.var3_indx=get(hdl.dataprep.var3,'value');
   para1.dataprep.data_dim=3;
   para1.dataprep.zlabel=para.dataprep.zlabel;
end

HDIR=pwd;
eval(['cd ''' para.file_dir.data_format_file ''''])
if opt == 1
  [pathi truncated_filename ext ver]=fileparts(para.dataprep.filename);
  fname_ini=[truncated_filename '_data_format'];
  [filename, filepath]=uiputfile('*.mat','Define a Data Format Filename',fname_ini);
  if ~isstr(filename) 
    eval(['cd ''' HDIR ''''])
    return;
  end
  para.dataprep.data_format_filename=filename;
  para.file_dir.data_format_file=filepath;
end
eval(['cd ''' para.file_dir.data_format_file ''''])
if para.Matlab_Version == 7
   cmd=['save ''' para.dataprep.data_format_filename '''' ' para1   -nounicode'];
else
   cmd=['save ''' para.dataprep.data_format_filename '''' ' para1'];
end
eval(cmd)
eval(['cd ''' HDIR ''''])

return