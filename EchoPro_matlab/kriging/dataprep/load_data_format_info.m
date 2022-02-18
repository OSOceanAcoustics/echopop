function   load_data_format_info(opt)
% function save_data_format_info.m load the data file information and the
% associated preliminary processing parameters
%   opt = not used
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

cmd=['load ''' para.krig.data_format_file ''''];
eval(cmd)

var=whos;
i=1;
while i <= length(var)
  if strcmp(var(i).name,'para1') 
    break
  end
  if i == length(var)
      krig_message(1,'Not a valid file format file, try to load a different file !!');
      return
  end
  i=i+1;
end

para.dataprep.reduct_fac=para1.dataprep.reduct_fac;
para.dataprep.filter_supt=para1.dataprep.filter_supt;
para.dataprep.filter_type=para1.dataprep.filter_type;				% filter type: 1-simple reduction,
para.dataprep.var1_indx=para1.dataprep.var1_indx;
para.dataprep.var2_indx=para1.dataprep.var2_indx;
para.dataprep.x_axis_indx=para1.dataprep.x_axis_indx;
para.dataprep.y_axis_indx=para1.dataprep.y_axis_indx;
if isfield(para1.dataprep,'ext_prog') 
   para.dataprep.ext_prog=para1.dataprep.ext_prog;
   if  para1.dataprep.ext_prog == 1
     para.dataprep.dat_conv_fname=para1.dataprep.dat_conv_fname;
     para.file_dir.data_conversion=para1.file_dir.data_conversion;
     cmd=['addpath ' para.file_dir.data_conversion];
     eval(cmd);
   end
end

para.dataprep.transform_index=para1.dataprep.transform_index;
para.dataprep.ytox=para1.dataprep.ytox;
para.dataprep.ztox=para1.dataprep.ztox;
para.dataprep.data_dim=para1.dataprep.data_dim;
para.dataprep.xlabel=para1.dataprep.xlabel;
para.dataprep.ylabel=para1.dataprep.ylabel;
if para.dataprep.data_dim == 3
   para.var3_indx=para1.var3_indx;
   para.dataprep.z_axis_indx=para1.dataprep.z_axis_indx;
   para.dataprep.zlabel=para1.dataprep.zlabel;
end



return