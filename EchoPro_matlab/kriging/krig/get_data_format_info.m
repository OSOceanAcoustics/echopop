function   get_data_format_info
% function get_data_format_info.m obtains the data file information and the
% associated preliminary processing parameters
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global para hdl data

cmd=['load ' para.krig.data_format_file];
eval(cmd)
para.dataprep.reduct_fac=para1.dataprep.reduct_fac;
para.dataprep.filter_supt=para1.dataprep.filter_supt;
para.dataprep.filter_type=para1.dataprep.filter_type;				% filter type: 1-simple reduction
para.dataprep.transform_index=para1.dataprep.transform_index;
para.dataprep.var1_indx=para1.dataprep.var1_indx;
para.dataprep.var2_indx=para1.dataprep.var2_indx;
para.dataprep.x_axis_indx=para1.dataprep.x_axis_indx;
para.dataprep.y_axis_indx=para1.dataprep.y_axis_indx;
para.dataprep.data_dim =para1.dataprep.data_dim ;
para.dataprep.dat_conv_fname=para1.dataprep.dat_conv_fname;
para.file_dir.data_conversion=para1.file_dir.data_conversion;
para.dataprep.ytox=para1.dataprep.ytox;
para.dataprep.ztox=para1.dataprep.ztox;

if para.dataprep.data_dim == 3
   para.dataprep.z_axis_indx=para1.dataprep.z_axis_indx;
end

return