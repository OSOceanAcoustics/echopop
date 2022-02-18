function			variogram_proc()
%% compute data-based  normalizedsemi-variogram/correlogram 
%%
%
%% INPUT
%    data.in.tv=f(x,y)		 (x,y) coordinates of transformed variable value
%    para.vario.res   = lag resolution
%    para.vario.range = range limit for semivariogram computation
%    para.vario
%            .angle = orientation angle from horizontal axis in degree
%            .ratio = aspect ratio of the scale in vertical (Y) to that of horizontal (X)
%% OUTPUT
%    data.out.vario.grd =  lag grid
%    data.out.vario.gammah = gamma(h) 	or gamma(grd)   semi-variogram or correlogram
%    data.out.vario.c0  =  variance at each lag
%    data.out.vario.cnt =  no. of pairs at each lag
%%
%%  Kriging Software Package  version 3.0,   December 29, 2001
%%  Copyright (c) 1998, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global para hdl data color


if para.status.variogram >= 1 & hdl.status.variogramfig == 1
   hdl_vario=findobj(hdl.vario.axes1,'type','line');
   delete(hdl_vario);
end
drawnow

%% data 
x=data.in.x;
y=data.in.y;
var=data.in.tv;
ni=length(var);				% number of data points of observations
if data.in.dim == 3
   z=data.in.z;
else
   z=zeros(ni,1);
end
indx_nan=find(isnan(var)|isnan(x)|isnan(y)|isnan(z));
var(indx_nan)=[];
x(indx_nan)=[];
y(indx_nan)=[];
z(indx_nan)=[];

if isfield(hdl.vario,'fileID')
    set(hdl.vario.fileID,'String',para.dataprep.fileID);
else
    hdl.vario.fileID = uicontrol('Parent',hdl.vario.h0, ...
        'Units','normalized', ...
        'BackgroundColor',color.dark_grey, ...
        'FontWeight','bold', ...
        'ForegroundColor',[0 0 0.501960784313725], ...
        'ListboxTop',0, ...
        'Position',[0.27 0.95 0.48 0.03], ...
        'String',para.dataprep.fileID, ...
        'Style','text', ...
        'Tag','VarioFileID');
end
% get all parameters needed to compute data-based variogram
para.vario.res=str2num(get(hdl.vario.res,'String'));
para.vario.azm_beg=str2num(get(hdl.vario.ang_beg_azm,'string'));
para.vario.azm_end=str2num(get(hdl.vario.ang_end_azm,'string'));
if para.vario.azm_beg < 0
    para.vario.azm_beg=0;
end
if para.vario.azm_end > 180
    para.vario.azm_end = 180;
end

para.vario.azm_res=str2num(get(hdl.vario.ang_res_azm,'string'));
para.vario.ang_wd_azm=str2num(get(hdl.vario.ang_wd_azm,'string'));
para.vario.azm_rot=str2num(get(hdl.vario.ang_rot_azm,'string'));

para.vario.dip_beg=str2num(get(hdl.vario.ang_beg_dip,'string'));
para.vario.dip_end=str2num(get(hdl.vario.ang_end_dip,'string'));
if para.vario.dip_beg < -90
    para.vario.dip_beg=-90;
end
if para.vario.dip_end > 90
    para.vario.dip_end = 90;
end

para.vario.dip_res=str2num(get(hdl.vario.ang_res_dip,'string'));
para.vario.ang_wd_dip=str2num(get(hdl.vario.ang_wd_dip,'string'));
para.vario.dip_rot=str2num(get(hdl.vario.ang_rot_dip,'string'));

para.vario.ytox_ratio=str2num(get(hdl.vario.ytox_ratio,'string'));
para.vario.ztox_ratio=str2num(get(hdl.vario.ztox_ratio,'string'));
para.vario.range=str2num(get(hdl.vario.range_edit,'string'));


res=para.vario.res;
range=para.vario.range;
nlag=round(range/res);

%% direction arrays
%azm0=para.vario.azm_beg:para.vario.azm_res:para.vario.azm_end;
%dip0=para.vario.dip_beg:para.vario.dip_res:para.vario.dip_end;
if get(hdl.vario.enabled,'value') == 1		% omni-directional
    azm0=0;
    dip0=0;
    azm=azm0;
    dip=dip0;
elseif get(hdl.vario.enable2d3d,'value') == 1   % anisotropic variogram
    azm0=para.vario.azm_beg:para.vario.azm_res:para.vario.azm_end;
    if data.in.dim == 2
        dip0=0;
    else
        dip0=para.vario.dip_beg:para.vario.dip_res:para.vario.dip_end;
    end
end
ndir_azm=length(azm0);
ndir_dip=length(dip0);
if ndir_azm == 1
    dang_azm=360;
    %	set(hdl.vario.ang_wd_azm,'string',dang_azm)
else
    dang_azm=str2num(get(hdl.vario.ang_wd_azm,'string'));
end
if ndir_dip == 1
    dang_dip=180;
    %	set(hdl.vario.ang_wd_dip,'string',dang_dip/2)
else
    dang_dip=str2num(get(hdl.vario.ang_wd_dip,'string'));
end

para.vario.nazm=length(azm0);
para.vario.ndip=length(dip0);

if data.in.dim == 3
    azm=reshape(azm0(ones(para.vario.ndip,1),:)',para.vario.nazm*para.vario.ndip,1);
    dip=reshape(dip0(ones(para.vario.nazm,1),:),para.vario.nazm*para.vario.ndip,1);
elseif get(hdl.vario.enabled,'value') ~= 1
    if para.vario.ndip > 1 & para.vario.nazm > 1
        err_msg={'Inconsistent parameter settings for anisotropic variogram/correlogram computation. ' ...
            '   ', ...
            'For a 2-D data set, an anisotropic variogram/correlogram can be computed with either a varying azimuth angle or a varying dip angle, but not both.'};
        h=krig_message(1,err_msg);
        return
    elseif para.vario.nazm > 1 & para.vario.ndip == 1
        azm=azm0(:);
        dip=dip0(1)*ones(para.vario.nazm,1);
    elseif para.vario.nazm == 1 & para.vario.ndip > 1
        dip=dip0(:);
        azm=azm0(1)*ones(para.vario.ndip,1);
    else
        krig_message(2,'Inconsistent parameter settings for anisotropic variogram/correlogram computation!')
    end
end


ndir=length(azm);
if para.vario.vario == 1
    ivtype=5;
elseif para.vario.corr == 1
    ivtype=4;
end

%% set default parameters
xltol=-1;											% distance tolerance -1 -> ras/2
atol=0.5*dang_azm*ones(ndir,1);
bandwh=0.5*range*ones(ndir,1);				% no banwidth constraint in azimuth
dtol=0.5*dang_dip*ones(ndir,1);
bandwd=0.5*range*ones(ndir,1);				% no banwidth constraint in dip
nvarg=para.vario.nvarg;
ivtail=para.vario.ivtail;
ivhead=para.vario.ivhead;
nvar=para.vario.nvar;
isill=para.vario.isill;

%% process anisotrophic data
if abs(para.vario.azm_rot) > eps |  abs(para.vario.dip_rot) > eps ...
        | abs(1-para.vario.ytox_ratio) > 1e-5 | abs(1-para.vario.ztox_ratio)
    %% observations
    cz_in=[x(:) y(:) z(:)];
    cz_out=coordtransform3d(cz_in);
    %	 figure;plot(x,y,'dr',cz(:,1),cz(:,2),'og');pause
    x=cz_out(:,1);y=cz_out(:,2);z=cz_out(:,3);
end

nnan=length(indx_nan);
nd=ni-nnan;				% usable number of data points

if nd > 50
    msg1='Computing Semivariogram/Correlogram, please be patient,... ';
    msg2=sprintf('usable data points = %g',nd);
    hdl0=krig_message(2,msg1,msg2);
    tic
else
    hdl0=[];
end
if nd < 1
    krig_message(1,'No usable data !!');
    return
end

%% call variogram function
[np, lag, gammah, hm, tm, hv, tv, xx, yy, zz] = variogram3d(nd, x, y, z, var, nlag, res, xltol,ndir, azm,...
    atol, bandwh, dip, dtol, bandwd, nvarg,ivtail, ivhead, ivtype, nvar,isill);

%% save variogram computation results
if get(hdl.vario.enabled,'value') == 1			% omni-directional
    data.out.vario.lag=lag;	% grid for plotting
    data.out.vario.gammah=gammah;
    data.out.vario.cnt=np;
    data.out.vario.hm=hm;
    data.out.vario.tm=tm;
    data.out.vario.hv=hv;
    data.out.vario.tv=tv;
    data.out.vario.xx=xx;
    data.out.vario.yy=yy;
    data.out.vario.zz=zz;
    data.out.vario.nlag=nlag;
    para.vario.max_sill=1.2*max(gammah);
    para.vario.sill=0.7*max(gammah);
    para.vario.max_nugt=max(gammah);
    set3dvariopara(2);
else
    data.out.vario.lag2d3d=lag;	% grid for plotting
    data.out.vario.gammah2d3d=gammah;
    data.out.vario.cnt2d3d=np;
    data.out.vario.hm2d3d=hm;
    data.out.vario.tm2d3d=tm;
    data.out.vario.hv2d3d=hv;
    data.out.vario.tv2d3d=tv;
    data.out.vario.xx=xx;
    data.out.vario.yy=yy;
    data.out.vario.zz=zz;
    data.out.vario.nlag2d3d=nlag;
    para.vario.max_sill=1.2*max(gammah);
    para.vario.sill=0.7*max(gammah);
    para.vario.max_nugt=max(gammah);
end
indx1=find(hv > 0 & tv > 0);
data.out.vario.c0_ht=mean(hv(indx1).*tv(indx1));
data.out.vario.c0=std(var)^2;
para.vario.c0=data.out.vario.c0;

close(hdl0);

if para.vario.nazm > 1 | para.vario.ndip > 1
    if para.vario.nazm> 1 & para.vario.ndip > 1
        para.vario.dim = 3;
    else
        para.vario.dim = 2;
    end
else
    para.vario.dim = 1;
end

%% plotting variogram
if get(hdl.vario.enabled,'value') == 1			% omni-directional
    plotvariogram1d(1);			%	plot 1-D data-based semi-variogram or correlogram
else
    plotvariogram2d3d(1);
end

% find the default length scale
if isfield(data.out.vario,'gammah')
    [tmp indx1]=max(data.out.vario.gammah(1:round(0.5*length(data.out.vario.gammah))));
    [tmp indxl]=min(abs(data.out.vario.gammah(1:indx1)-0.3));			% 6dB
    para.vario.lscl=data.out.vario.lag(indxl);
    if para.vario.lscl <= 0.01*max(data.out.vario.lag)
        para.vario.lscl=0.1*max(data.out.vario.lag);
    end
elseif isfield(data.out.vario,'gammah2d3d')
    [tmp indx1]=max(data.out.vario.gammah2d3d(1:round(0.5*length(data.out.vario.gammah2d3d))));
    [tmp indxl]=min(abs(data.out.vario.gammah2d3d(1:indx1)-0.3));			% 6dB
    para.vario.lscl=data.out.vario.lag2d3d(indxl);
    if para.vario.lscl <= 0.01*max(data.out.vario.lag2d3d)
        para.vario.lscl=0.1*max(data.out.vario.lag2d3d);
    end
end
set(hdl.vario.lscl_edit,'String',num2str(para.vario.lscl));
set(hdl.vario.lscl_slider,'Value',para.vario.lscl/para.vario.max_lscl);
if isfield(data.out.vario,'lag_theo')
    para.status.variogram=2;
else
    para.status.variogram=1;
end
