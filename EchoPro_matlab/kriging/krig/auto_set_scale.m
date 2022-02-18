function auto_set_scale(hdl,para)
  xmin=str2num(get(hdl.krig.xmin,'string'));
  xmin_str=num2str(floor(xmin/para.krig.x_res)*para.krig.x_res);
  set(hdl.krig.xmin,'string',xmin_str)
  xmax=str2num(get(hdl.krig.xmax,'string'));
  xmax_str=num2str(ceil(xmax/para.krig.x_res)*para.krig.x_res);
  set(hdl.krig.xmax,'string',xmax_str)

  ymin=str2num(get(hdl.krig.ymin,'string'));
  ymin_str=num2str(floor(ymin/para.krig.y_res)*para.krig.y_res);
  set(hdl.krig.ymin,'string',ymin_str)
  ymax=str2num(get(hdl.krig.ymax,'string'));
  ymax_str=num2str(ceil(ymax/para.krig.y_res)*para.krig.y_res);
  set(hdl.krig.ymax,'string',ymax_str)

  dx_str=para.krig.x_res;
  set(hdl.krig.dx,'string',para.krig.x_res)

  dy_str=para.krig.y_res;
  set(hdl.krig.dy,'string',para.krig.y_res)
return