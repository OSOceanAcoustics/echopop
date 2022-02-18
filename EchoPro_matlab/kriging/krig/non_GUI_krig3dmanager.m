function		non_GUI_krig3dmanager()
%% get kriging parameters from the Kriging GUI panel and perfrom kriging
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

DEG2RAD=pi/180;

if para.proc.default_para ~= 1
  set3dkrigpara(2);	
end

%% variogram/correlogram model parameters
if ~isfield(para.vario,'model')
  krig_message(1,'Need to load parameter file first !!');
  return
end
vmod=para.vario.model;
if para.vario.corr == 1
  vmod = -vmod;
end
model_para=[para.vario.nugt para.vario.sill para.vario.lscl para.vario.powr para.vario.hole];

if para.proc.default_para ~= 1
    %% check whether the selected kriging region is reasonable
    xmin_set=str2num(get(hdl.krig.xmin,'string'));
    xmax_set=str2num(get(hdl.krig.xmax,'string'));
    if min(data.in.var1) > xmax_set | max(data.in.var1) < xmin_set
        krig_message(1,['X-coordinate settings are inconsistent with the input data, click on ''Refresh'' button to ' ...
            'get automatic settings !!']);
        return
    end
    ymin_set=str2num(get(hdl.krig.ymin,'string'));
    ymax_set=str2num(get(hdl.krig.ymax,'string'));
    if min(data.in.var2) > ymax_set | max(data.in.var2) < ymin_set
        krig_message(1,['Y-coordinate settings are inconsistent with the input data, click on ''Refresh'' button to ' ...
            'get automatic settings !!']);
        return
    end
    if data.in.dim == 3
        zmin_set=str2num(get(hdl.krig.zmin,'string'));
        zmax_set=str2num(get(hdl.krig.zmax,'string'));
        if min(data.in.var3) > zmax_set | max(data.in.var3) < zmin_set
            krig_message(1,['Z-coordinate settings are inconsistent with the input data, click on ''Refresh'' button to ' ...
                'get automatic settings !!']);
            return
        end
    end
end
%% kriging parameters
mod_opt=para.krig.model;
srad=para.krig.srad;
kmin=para.krig.kmin;
kmax=para.krig.kmax;
blk_opt=para.krig.scheme;
blknx=para.krig.blk_nx;
blkny=para.krig.blk_ny;
blknz=para.krig.blk_nz;
elim=para.krig.elim;
dx=para.krig.dx;
dy=para.krig.dy;
x=data.in.x;
y=data.in.y;
var=data.in.tv;
nx=round((para.krig.xmax-para.krig.xmin)/dx);
ny=round((para.krig.ymax-para.krig.ymin)/dy);
xg=linspace(para.krig.xmin,para.krig.xmax,nx);
yg=linspace(para.krig.ymin,para.krig.ymax,ny);
para.krig.nx=nx;
para.krig.ny=ny;

data.out.krig.nx=nx;
data.out.krig.ny=ny;
data.out.krig.xg=xg;
data.out.krig.yg=yg;

if data.in.dim == 3 
	dz=para.krig.dz;
	z=data.in.z;
	nz=round((para.krig.zmax-para.krig.zmin)/dz);
	data.out.krig.nz=nz;
 	para.krig.nz=nz;
    zg=linspace(para.krig.zmin,para.krig.zmax,nz);
	data.out.krig.zg=zg;
end


msg1='Kriging in progress, ...';
msg2=' ';
hmsg=krig_message(2,msg1,msg2);

t0=clock;
%get_set_gridfile_para;

%% find the west bound boundary from SCB to NWVI  
[xb,yb]=find_west_bound(data);

if data.in.dim == 2
    para.dataprep.latlonfac_g=cos(DEG2RAD*data.out.krig.gy);
    xp=(data.out.krig.gx-para.dataprep.x_offset).*para.dataprep.latlonfac_g./para.dataprep.x_norm;
    yp=(data.out.krig.gy-para.dataprep.y_offset)/para.dataprep.y_norm;
    [vp,ep,ep_s]=krig2d(xp,yp,x,y,var,model_para,hmsg,xb,yb);
else
   [X,Y,Z]=meshgrid(xg,yg,zg);
   xp=reshape(X,nx*ny*nz,1);
   yp=reshape(Y,nx*ny*nz,1);
   zp=reshape(Z,nx*ny*nz,1);
	if para.krig.model == 1			% remove the mean for simple kriging
		var_mean=mean(var);
		var=var-var_mean;
   end
   [vp,ep,ep_s]=krig3d(xp,yp,zp,x,y,z,var,model_para,hmsg);
	if para.krig.model == 1			% add in the mean for simple kriging
		vp=vp+var_mean;
	end
   Vg=reshape(vp,ny,nx,nz);
   Eg=reshape(ep,ny,nx,nz);
   Eg_s=reshape(ep_s,ny,nx,nz);
end
%% customized grid data
if para.krig.load_griddata_file == 1
   data.out.krig.gv=vp;
   data.out.krig.ge=ep;
   data.out.krig.ge_s=ep_s;
end

close(hmsg);

%% ensure the kriged values have the same sign as the input values
t1=etime(clock,t0);

if para.krig.batch_file_proc == 1
 %%% 1/5/02   
  if para.krig.bat_proc_cnt == 1
     para.krig.bat_proc_time=round(t1);
     msg1=[sprintf('Kriging on file # %d is done in %g (sec), \n\nTotal time = %g (sec)',para.krig.bat_proc_cnt,round(t1),para.krig.bat_proc_time)];
     hdl.krig.hmsg=krig_message(1,msg1);
  else
     msg1=[sprintf('Kriging on file # %d is done in %g (sec), \n\nTotal time = %g (sec)',para.krig.bat_proc_cnt,round(t1),para.krig.bat_proc_time)];
     hdl.krig.hmsg=krig_message(3,msg1,' ',hdl.krig.hmsg);
     para.krig.bat_proc_time=para.krig.bat_proc_time+round(t1);
  end   
  
else
   msg1=['Kriging is done in ' sprintf('%g (sec)',round(t1))];
   hdl.krig.hmsg=krig_message(1,msg1);
end

para.krig.hmsg=hdl.krig.hmsg;
para.status.kriging = 1;

para.dispkrig.Qcheck=0;
para.dispkrig.JKcheck=0;
para.dispkrig.Dblkrig=0;

%% save outputs to data.out structure
data.out.krig.xp=xp(:);
data.out.krig.yp=yp(:);
data.out.krig.var=vp(:);
data.out.krig.err=ep(:);        % kriging varaince
data.out.krig.err_s=ep_s(:);   % sample variance 
data.out.krig.xn=x(:);
data.out.krig.yn=y(:);
data.out.krig.vn=var(:);
if data.in.dim==3
	data.out.krig.zp=zp(:);
	data.out.krig.zn=z(:);
end
para.krig.proc_opt=2;
if  para.krig.loop ~= 1
    kriging_proc
else
    return
end

