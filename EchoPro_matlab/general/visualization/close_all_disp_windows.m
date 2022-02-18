function  close_all_disp_windows
%% close all popup visualization windows 
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

fig_num=findobj('type','figure');
[fig_hdl,n]=sort(fig_num);
fig_hdl=fig_hdl(fig_hdl<= 14);
close(fig_hdl)

return