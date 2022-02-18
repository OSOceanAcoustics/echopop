function Strout=modelString(index)
%% function Strout=modelString
%% builds variogram model help file strings
%%
%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl

switch index
  case 1
      model_name='Spherical:';
      expression='\gamma_h = C_0[1.5h/L-0.5(h/L)^3]+ \gamma_0';
      var4='L = length scale';
  case 2
      model_name='Exponential:';
      expression='\gamma_h = C_0(1-e^{-h/L}) + \gamma_0';
      var4='L = length scale';
  case 3
      model_name='Gaussian:';
      expression='\gamma_h = C_0(1-e^{-(h/L)^2}) + \gamma_0';
      var4='L = length scale';
  case 4
      model_name='Linear:';
      expression='\gamma_h = C_0 h + \gamma_0';
      var4=' ';
  case 5
      model_name='Sinc:';
      expression='\gamma_h = C_0 [1-sin(bh)] + \gamma_0';
      var4='b = length scale of hole effect';
  case 6
		model_name='Exponential-Cosine(type I):';
      expression='\gamma_h = C_0 [1 - cos(b h) e^{-h/L} ] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 7
		model_name='Exponential-Cosine(type II):';
      expression='\gamma_h = C_0 [1 +  cos(b h)e^{-h/L}] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 8
		model_name='Gaussian-Cosine:';
      expression='\gamma_h = C_0 [1- cos(b h)e^{-(h/L)^2}] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 9
		model_name='Bessel(J_0):';
      expression='\gamma_h = C_0 ( 1 - J_0(bh)) + \gamma_0';
      var4='b = length scale of hole effect';
  case 10
		model_name='Exponential-Bessel:';
      expression='\gamma_h = C_0 [1- J_0(b h)e^{-h/L}] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 11
		model_name='Gaussian-Bessel:';
      expression='\gamma_h = C_0 [1- J_0(b h)e^{-(h/L)^2}] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 12
		model_name='Gaussian-Linear:';
      expression='\gamma_h = C_0 [1- (1- b h) e^{-(h/L)^2}] + \gamma_0';
      var4='L = length scale';
      var5='b = length scale of hole effect';
  case 13
     	model_name='General Exponential-Bessel:';
     	expression='\gamma_h = C_0 [ 1 - J_0(bh)e^{-(h/L)^p}] + \gamma_0';
     	var4='L = length scale';
     	var5='b = length scale of hole effect,      p = power';
end

if isfield(hdl,'help')
 if ~isempty(hdl.help) 
   figure(hdl.help.h0);
   clf;
 else
   hdl.help.h0 = figure('Color',[0.8 0.8 0.8], ...
      'Units','normalized', ...
		'name','   VARIOGRAM  MODEL', ...
      'Position',[0.1 0.2 0.4 0.3]);
   set(gcf,'menubar','none');
 end
else
 hdl.help.h0 = figure('Color',[0.8 0.8 0.8], ...
      'Units','normalized', ...
		'name','   VARIOGRAM  MODEL', ...
      'Position',[0.1 0.2 0.4 0.3]);
 set(gcf,'menubar','none');
end

h1 = uicontrol('Parent',hdl.help.h0, ...
	'Units','normalized', ...
	'Callback','hdl.help=[];close', ...
	'FontSize',8, ...
	'FontWeight','bold', ...
	'ListboxTop',0, ...
	'Position',[0.80  0.05 0.1  0.08], ...
	'String','Quit', ...
	'Tag','Quit');

ht1=text(0.02,0.95,model_name,'sc');
ht2=text(0.2,0.85,expression,'sc');
ht3a=text(0.04,0.7,'where:','sc');
ht3b=text(0.3,0.7,'h = lag,      \gamma_0 = nugget','sc');
ht4a=text(0.3,0.6,['C_0 = sill - nugget,     ' var4],'sc');
if exist('var5')
  ht4b=text(0.3,0.5,var5,'sc');
  ht4=[ht4a ht4b];
else
  ht4=ht4a;
end
ht5=text(0.02,0.3,'Relation between semi-variogram and correlogram:','sc');
ht6=text(0.3,0.2,'C_h = 1 - \gamma_h','sc');

ht=[ht1 ht2 ht3a ht3b ht4 ht5 ht6];
set(ht,'fontweight','bold')
axis off