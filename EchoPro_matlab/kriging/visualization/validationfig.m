function validationfig()
%%  function validationfig() establishes the Validation GUI window
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para data 

grey=[0.75 0.75 0.75];
dark_grey=[0.65 0.65 0.65];

default_setting=0;
if ~isfield(data, 'in') 
   data.in.dim=3;		% default setting
   default_setting=1;
end

   
if ~isempty(findobj('type','figure','Tag','Validation')) 
  figure(hdl.validation.h0);
  return
end


EPS=0.001;

if ~isempty(findobj('type','figure','Tag','Validation'))
  figure(hdl.validation.h0);
  set(hdl.validation.h0,'Name','Validation');
end

%	
hdl.validation.h0= figure('Units','normalized', ...
	'Color',[0.8 0.8 0.8], ...
	'Name','Validation', ...
   'Renderer','zbuffer', ...
   'Tag','Validation', ...
   'NumberTitle','off',...
	'RendererMode','manual','Position',[0.3 0.3 0.2 0.2]);
   set(0, 'showhidden', 'on')
   ch=get(gcf, 'children');
%	delete(ch(1))								%Help
   wmhdl=findobj(ch,'Label','&Help');		
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%	delete(ch(3))								%Tools
   wmhdl=findobj(ch,'Label','&File');		
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%     new feature in V6.x delete(ch(6))			%Edit
   wmhdl=findobj(ch,'Label','&Edit');		
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%      new feature in V6.x                      %insert
   wmhdl=findobj(ch,'Label','&Insert');		
   if ~isempty(wmhdl)
     delete(wmhdl);
     ch(find(ch == wmhdl))=[];
   end
%      new feature in V6.x                      %View
   wmhdl=findobj(ch,'Label','&View');		
   if ~isempty(wmhdl)
     delete(wmhdl);
     ch(find(ch == wmhdl))=[];
   end
%      new feature of V7.0                      
   wmhdl=findobj(ch,'Label','&Desktop');	    %Desktop
   if ~isempty(wmhdl)
     delete(wmhdl);
     ch(find(ch == wmhdl))=[];
   end

hdl.validation.help=uimenu(hdl.validation.h0,'label','&Help','separator','off');
uimenu(hdl.validation.help,'label','   Q-1','Callback','visualization_help(5)','separator','on');
uimenu(hdl.validation.help,'label','   Q-2','Callback','visualization_help(5)','separator','off');
uimenu(hdl.validation.help,'label','   Double Kriging','Callback','visualization_help(7)','separator','off');
uimenu(hdl.validation.help,'label','   Leaving-one-out','Callback','visualization_help(8)','separator','off');
hdl_visualization_help_CMPT=uimenu(hdl.validation.help,'label','   Compute/Re-Compute','Callback','visualization_help(9)','separator','on');
hdl_quit=uimenu(hdl.validation.h0,'label','Quit','Callback','close');

x0=0.2;y0=0.7;txt_w=0.1;pmenu_w=0.1;
model_str={'1.  Q-1','1.  Q-2','3.  double kriging','4.  leaving-one-out',};
h1 = uicontrol('Parent',hdl.validation.h0, ...
	'Units','normalized', ...
	'BackgroundColor',grey, ...
	'HorizontalAlignment','center', ...
	'FontWeight','bold', ...
	'Position',[x0 y0 0.5 txt_w], ...
	'String',' Method', ...
	'Style','text');
hdl.validation.model = uicontrol('Parent',hdl.validation.h0, ...
	'Units','normalized', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','validation_proc(0)', ...
	'Position',[x0 y0-txt_w 0.5 pmenu_w], ...
	'String',model_str, ...
	'Style','popupmenu', ...
	'Tag','ValidationModel');

%% Push button
h1 = uicontrol('Parent',hdl.validation.h0, ...
	'Units','normalized', ...
	'Callback','validation_proc(0)', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
	'Position',[0.2 0.34 0.40 0.1], ...
	'String','Compute');
hdl.validation.recompute = uicontrol('Parent',hdl.validation.h0, ...
	'Units','normalized', ...
	'Callback','validation_proc(1)', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
   'Position',[0.2 0.2 0.40 0.1], ...
   'enable','off', ...
	'String','Re-Compute');
h1 = uicontrol('Parent',hdl.validation.h0, ...
	'Units','normalized', ...
	'Callback','close', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
	'Position',[0.65 0.27 0.3 0.1], ...
	'String','Quit');

if default_setting ==1 
   data.in.dim=[]; 		% default setting
   default_setting=0;
end


para.status.validation = 1;
hdl.status.validation = 1;


