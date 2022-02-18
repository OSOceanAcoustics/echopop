function  lsqfit()
%% function	lsqfit() performs the Nonlinear Least-Square-Fit function
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data

ndigits=5;
if para.status.variogram >= 1			% Data-based Semi-Variogram or Correlogram has been computed
para.vario.model=get(hdl.vario.model,'Value');

indx=find(isnan(data.out.vario.gammah));
if ~isempty(indx)
  for i=1:length(indx)
	 data.out.vario.gammah(indx(i))=mean_nan(data.out.vario.gammah);
  end
end
cnt=data.out.vario.cnt;
wgt=cnt/sum(cnt);
data_in=[data.out.vario.lag(:) data.out.vario.gammah(:) wgt(:)];
model_para=[para.vario.nugt, para.vario.sill, para.vario.lscl, ... 
				para.vario.powr, para.vario.hole];

%% get truncated data
range=str2num(get(hdl.vario.range_edit,'String'));
indx=find(data.out.vario.lag <= range);
truncated_data=data_in(indx,:);
option=optimset('fminbnd');
option.MaxIter=2000;
LB=zeros(1,5);
UB=[1 3 2 4 10];
option.LevenbergMarquardt=1;
%[fit_para,option]=leastsq('cost_function',model_para,[],[],truncated_data);
[fit_para,norm2]=lsqnonlin('cost_function',model_para,LB,[],option,truncated_data);
%% ignore the imaginary partof the fitting coefficeint
para.vario.nugt=abs(real(fit_para(1)));
para.vario.sill=abs(real(fit_para(2)));
para.vario.lscl=abs(real(fit_para(3)));
para.vario.powr=abs(real(fit_para(4)));
para.vario.hole=abs(real(fit_para(5)));
variogram_theo(2);
else
   krig_message(1,['Select ' '''' 'Compute' '''' ' Button to compute Data-based Semi-Variogram/Correlogram from the data']);  
end