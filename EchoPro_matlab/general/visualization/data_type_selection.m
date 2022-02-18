function data_type_selection(hObject, eventdata, handles)
%% change enabled variable display based on the data type selection:un-kriged, kriged, and biological
%% GUI -based operation
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       11/30/2013

global hdl 

if get(hdl.vis.radio_visual_un_kriged,'value') == 1  % using values prior to kriging
   set(handles.radio_visual_func_NASC,'Enable','on')
   if get(handles.radio_visual_func_length,'value')
       set(handles.radio_visual_func_length,'value',0);
       set(handles.radio_visual_func_NASC,'value',1);
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_age,'value')
       set(handles.radio_visual_func_age,'value',0);
       set(handles.radio_visual_func_NASC,'value',1);
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_NASC,'value')
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_number,'value')
       hObject=handles.radio_visual_func_number;
   elseif get(handles.radio_visual_func_biomass,'value')
       hObject=handles.radio_visual_func_biomass;
   end
elseif get(hdl.vis.radio_visual_kriged,'value') == 1  % using the kriged values
   set(handles.radio_visual_func_NASC,'Enable','on')
   if get(handles.radio_visual_func_length,'value')
       set(handles.radio_visual_func_length,'value',0);
       set(handles.radio_visual_func_NASC,'value',1);
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_age,'value')
       set(handles.radio_visual_func_age,'value',0);
       set(handles.radio_visual_func_NASC,'value',1);
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_NASC,'value')
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_number,'value')
       hObject=handles.radio_visual_func_number;
   elseif get(handles.radio_visual_func_biomass,'value')
       hObject=handles.radio_visual_func_biomass;
   end      
elseif get(hdl.vis.radio_visual_biological,'value') == 1  % using biological trawl data
   if get(handles.radio_visual_func_length,'value')
       hObject=handles.radio_visual_func_length;
   elseif get(handles.radio_visual_func_age,'value')
       hObject=handles.radio_visual_func_age;
   elseif get(handles.radio_visual_func_NASC,'value')
       set(handles.radio_visual_func_NASC,'value',0);
       set(handles.radio_visual_func_number,'value',1);
       hObject=handles.radio_visual_func_NASC;
   elseif get(handles.radio_visual_func_number,'value')
       hObject=handles.radio_visual_func_number;
   elseif get(handles.radio_visual_func_biomass,'value')
       hObject=handles.radio_visual_func_biomass;
   end      
end
set_enabled_selection(hObject, eventdata, handles);
