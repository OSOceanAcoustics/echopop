function getdefault3dkrigpara
% Get DeFault Kriging Parameters

%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global hdl para data

% default grid sizes
para.krig.nx=20;
para.krig.ny=20;
para.krig.nz=20;
if para.status.dataprep == 1
	para.krig.xmin0=min(data.in.x0);
	para.krig.xmax0=max(data.in.x0);
	para.krig.ymin0=min(data.in.y0);
	para.krig.ymax0=max(data.in.y0);
	para.krig.zmin0=min(data.in.z0);
	para.krig.zmax0=max(data.in.z0);
    para.krig.dx0=(para.krig.xmax0-para.krig.xmin0)/(para.krig.nx-1);		
    para.krig.dy0=(para.krig.ymax0-para.krig.ymin0)/(para.krig.ny-1);
    para.krig.dz0=(para.krig.zmax0-para.krig.zmin0)/(para.krig.nz-1);
else
	para.krig.xmin0=-0.5;
	para.krig.xmax0=0.5;
	para.krig.ymin0=-0.5;
	para.krig.ymax0=0.5;
	para.krig.zmin0=-0.5;
	para.krig.zmax0=0.5;
    para.krig.dx0=0.05;
    para.krig.dy0=0.05;
    para.krig.dz0=0.05;
end

para.krig.dx=0.05;		
para.krig.dy=0.05;
para.krig.dz=0.05;

%% kriging parameters
para.krig.model=2;
para.krig.scheme=1;
para.krig.blk_nx=3;
para.krig.blk_ny=3;
para.krig.blk_nz=3;
para.krig.srad=0.3;
para.krig.kmin=10;
para.krig.kmax=30;
para.krig.elim=2.5;
para.krig.eps=1e-14;
para.krig.ratio=1e-3;

%% loading/saving parameters
para.krig.load_para=0;
para.krig.save_para=0;
para.krig.vario_para=0;
para.krig.krig_para=0;
para.krig.both_para=1;
para.krig.para_file_in='';
para.krig.para_file_out='';
para.krig.load_data_file=0;
para.krig.data_file='';

set3dkrigpara(1);				% set kriging parameter from variable struct


