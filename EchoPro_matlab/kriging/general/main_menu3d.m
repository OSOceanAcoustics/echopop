function fig = main_menu3d()
%%  base window
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl color

load window_position
hdl.window_position=win_pos;

if ~isempty(findobj('type','figure','Tag','Navigator'))
  figure(hdl.navigator.h0);
  return
end

hdl.navigator.h0 = figure('Units','normalized', ...
	'Color',color.background, ...
	'Name','Navigator', ...
	'Position',win_pos, ...
    'NumberTitle','off', ...
	'Tag','Navigator');

    set(0, 'showhidden', 'on');
    ch=get(gcf, 'children');
%   delete(ch(1));                              %File
   wmhdl=findobj(ch,'Label','&File');		    
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%	delete(ch(2))								%Help
   wmhdl=findobj(ch,'Label','&Help');		    
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%	delete(ch(4))								%Tools
   wmhdl=findobj(ch,'Label','&Tools');		
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%   delete(ch(5))								%Edit
   wmhdl=findobj(ch,'Label','&Edit');		
   delete(wmhdl);
   ch(find(ch == wmhdl))=[];
%      new feature of V6.x                      %insert
   wmhdl=findobj(ch,'Label','&Insert');	
   if ~isempty(wmhdl)
     delete(wmhdl);
     ch(find(ch == wmhdl))=[];
   end
%      new feature of V6.x                      %View
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

%%    Delete Toolbar
   toolbar_hdl=findobj(ch,'Type','uitoolbar');
   delete(toolbar_hdl);
   ch(find(ch == toolbar_hdl))=[];

h1 = uicontrol('Parent',hdl.navigator.h0, ...
	'Units','normalized', ...
	'Callback','copyright', ...
	'FontSize',16, ...
	'BackgroundColor',color.background, ...
	'FontWeight','bold', ...
	'ForegroundColor',[0 0 1], ...
	'Position',[0.26 0.9 0.5 0.06], ...
	'String','EasyKrig  V3.0', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',hdl.navigator.h0, ...
	'Units','normalized', ...
	'BackgroundColor',color.background, ...
	'FontSize',14, ...
	'FontWeight','bold', ...
	'ForegroundColor',[0 0 0], ...
	'ListboxTop',0, ...
	'Position',[0.1 0.85 0.8 0.046], ...
	'String','WOODS HOLE OCEANOGRAPHIC INSTITUTION', ...
	'Style','text', ...
	'Tag','StaticText2');


h1 = uicontrol('Parent',hdl.navigator.h0, ...
	'Units','normalized', ...
	'Callback','close_window(0:5);', ...
	'FontSize',10, ...
	'FontWeight','bold', ...
	'ListboxTop',0, ...
	'Position',[0.82 0.1 0.1 0.06], ...
	'String','Quit', ...
	'Tag','Quit');

h2(1)= uimenu('Parent',hdl.navigator.h0,'label','&Task');
h2(2)=uimenu(h2(1),'label','   &Load Data','callback','dataprep3dfig;','separator','off');
h2(3)=uimenu(h2(1),'label','   &Variogram','callback','variogram3dfig;','separator','off');
h2(4)=uimenu(h2(1),'label','   &Kriging','callback','kriging3dfig;','separator','off');
h2(5)=uimenu(h2(1),'label','   &Visualization','callback','dispkrig3dfig;','separator','off');
h2(6)= uimenu(h2(1),'label',' &Save Window Position','callback','save_window_pos(hdl.navigator.h0);','separator','on');


h3(1)=uimenu('Parent',hdl.navigator.h0,'label','&About');
h3(2)= uimenu(h3(1),'label','&EasyKrig3.0','callback','program_info(1);');
h3(3)= uimenu(h3(1),'label','&Kriging','callback','program_info(2);');
h3(4)= uimenu(h3(1),'label','&ReadMe','callback','program_info(3);');
h3(5)= uimenu(h3(1),'label',' &Variable Structure','callback','edit variables.txt','separator','on');

h4=uimenu('Parent',hdl.navigator.h0,'label','&Help');
h4a= uimenu(h4,'label','&Task');
uimenu(h4a,'label','   &Load Data','callback','navigator_help(1)','separator','off');
uimenu(h4a,'label','   &Variogram','callback','navigator_help(1)','separator','off');
uimenu(h4a,'label','   &Kriging','callback','navigator_help(1)','separator','off');
uimenu(h4a,'label','   &Visualization','callback','navigator_help(1)','separator','off');
uimenu(h4a,'label','   &Save Window Position','callback','navigator_help(2)','separator','on');
uimenu(h4a,'label','   &Quit','callback','navigator_help(2)','separator','off');

h5= uimenu('Parent',hdl.navigator.h0,'label','&Quit','callback','close_window(0:5);');

if nargout > 0, fig = hdl.navigator.h0; end

