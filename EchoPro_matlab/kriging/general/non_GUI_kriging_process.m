function non_GUI_kriging_process

global para	 		% all setting and processing parameters
global data  		% input and output data

auto_load_vario_krig_para(para.krig.vario_krig_para_filename);

non_GUI_load_data
non_GUI_compute_vario
non_GUI_kriging

return

