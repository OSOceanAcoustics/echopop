function  [model,nugt,sill,L,p,b]=get3dvariopara_theo()
%% get parameter from the panel edit field (set by user) of semi-variogram/correlogram
%% to compute model-based variogram (theory)
%
%% Last Modified December 22, 2001
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global hdl para

%%% read in parameters
model=get(hdl.vario.model,'Value');
value=get(hdl.vario.correlogram,'value');
if value == 1
   model=-model;						% correlogram model
end
nugt=str2num(get(hdl.vario.nugt_edit,'String'));
sill=str2num(get(hdl.vario.sill_edit,'String'));
L=str2num(get(hdl.vario.lscl_edit,'String'));
p=str2num(get(hdl.vario.powr_edit,'String'));
b=str2num(get(hdl.vario.hole_edit,'String'));

if nugt > para.vario.max_nugt para.vario.max_nugt=1.2*nugt;end
if sill > para.vario.max_sill para.vario.max_sill=1.2*sill;end
if L > para.vario.max_lscl para.vario.max_lscl=1.2*L;end
if p > para.vario.max_powr para.vario.max_powr=1.2*p;end
if b > para.vario.max_hole para.vario.max_hole=1.2*b;end

para.vario.nugt=nugt;
para.vario.sill=sill;
para.vario.lscl=L;
para.vario.powr=p;
para.vario.hole=b;
para.vario.model=abs(model);
para.vario.range=str2num(get(hdl.vario.range_edit,'String'));
para.vario.res=str2num(get(hdl.vario.res,'String'));
para.vario.vario=get(hdl.vario.variogram,'value');
para.vario.corr=get(hdl.vario.correlogram,'value');

