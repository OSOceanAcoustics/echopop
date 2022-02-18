function set_enabled_selection(hObject, eventdata, handles )
% set displayed variables 'enabled' property to either ON or OFF
% depending on whether 'kriged' is ON or OFF
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global hdl para data

if ~isfield(data,'final') | ~isfield(data.final,'table')
    return
end
if get(hdl.vis.radio_visual_un_kriged,'value') == 1  % using values prior to kriging
    switch get(hObject,'String')
        case 'Length'
            set(handles.radio_visual_var_age,'enable','off')
            set(handles.radio_visual_var_bot_depth,'enable','off')
            set(handles.radio_visual_var_layer_depth,'enable','off')
            set(handles.radio_visual_var_strata,'enable','off')
            set(handles.radio_visual_var_trawl,'enable','off')
            set(handles.radio_visual_var_lat,'enable','off')
            set(handles.radio_visual_var_transect,'enable','off')
            set(handles.radio_visual_var_gender,'enable','off')
            set(handles.radio_visual_var_weight,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','off')
            set(handles.radio_visual_var_length,'enable','off')
            set(handles.radio_visual_var_lat_lon,'enable','of')
            set(handles.radio_visual_var_len_age,'enable','off')
            %             data.visual.var1=data.final.table.trawl(:,9);
            %             data.visual.var1_str='Length (cm)';
        case 'Age'
            set(handles.radio_visual_var_age,'enable','off')
            set(handles.radio_visual_var_bot_depth,'enable','off')
            set(handles.radio_visual_var_layer_depth,'enable','off')
            set(handles.radio_visual_var_strata,'enable','off')
            set(handles.radio_visual_var_trawl,'enable','off')
            set(handles.radio_visual_var_lat,'enable','off')
            set(handles.radio_visual_var_transect,'enable','off')
            set(handles.radio_visual_var_gender,'enable','off')
            set(handles.radio_visual_var_weight,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','off')
            set(handles.radio_visual_var_length,'enable','off')
            set(handles.radio_visual_var_lat_lon,'enable','of')
            set(handles.radio_visual_var_len_age,'enable','off')
            %             data.visual.var1=data.final.table.trawl(:,11);
            %             data.visual.var1_str='Age';
        otherwise
            set(handles.radio_visual_var_bot_depth,'enable','on')
            set(handles.radio_visual_var_layer_depth,'enable','on')
            set(handles.radio_visual_var_strata,'enable','on')
            set(handles.radio_visual_var_weight,'enable','off')
            set(handles.radio_visual_var_lat,'enable','on')
            set(handles.radio_visual_var_transect,'enable','on')
            set(handles.radio_visual_var_gender,'enable','on')
            set(handles.radio_visual_var_survey_region,'enable','on')
            set(handles.radio_visual_var_length,'enable','on')
            set(handles.radio_visual_var_lat_lon,'enable','on')
            set(handles.radio_visual_var_age,'enable','on')
            set(handles.radio_visual_var_len_age,'enable','on')
            switch get(hObject,'String')
                case 'NASC'
                    set(handles.radio_visual_var_age,'enable','off')
                    set(handles.radio_visual_var_len_age,'enable','off')
                    set(handles.radio_visual_var_gender,'enable','off')
                    set(handles.radio_visual_var_length,'enable','off')
                    data.visual.var1=data.final.table.biomass(:,9);
                    data.visual.var1_str='NASC (x1000 m^2 nmi^2)';
                case 'Abundance'
                    data.visual.var1=data.final.table.biomass(:,10:12);
                    data.visual.var1_str='Abundance';
                case 'Biomass'
                    data.visual.var1=data.final.table.biomass(:,13:15);
                    data.visual.var1_str='Biomass (x 1000 mt)';
            end
    end
elseif get(hdl.vis.radio_visual_kriged,'value') == 1  % using the kriged values
    set(handles.radio_visual_var_bot_depth,'enable','off')
    set(handles.radio_visual_var_layer_depth,'enable','off')
    set(handles.radio_visual_var_strata,'enable','on')
    set(handles.radio_visual_var_trawl,'enable','off')
    set(handles.radio_visual_var_lat,'enable','on')
    set(handles.radio_visual_var_transect,'enable','off')
    set(handles.radio_visual_var_gender,'enable','on')
    set(handles.radio_visual_var_survey_region,'enable','on')
    set(handles.radio_visual_var_length,'enable','on')
    set(handles.radio_visual_var_lat_lon,'enable','on')
    set(handles.radio_visual_var_age,'enable','on')
    set(handles.radio_visual_var_len_age,'enable','on')
    set(handles.radio_visual_var_weight,'enable','off')
    switch get(hObject,'String')
        case 'Length'
            set(handles.radio_visual_var_strata,'enable','off')
            set(handles.radio_visual_var_length,'enable','off')
            set(handles.radio_visual_var_lat,'enable','off')
            set(handles.radio_visual_var_lat_lon,'enable','off')
            set(handles.radio_visual_var_gender,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','off')
            set(handles.radio_visual_var_age,'enable','off')
            set(handles.radio_visual_var_len_age,'enable','off')
    %            data.visual.var1=data.final.table.trawl(:,9);
            data.visual.var1_str='Length (cm)';
        case 'Age'
            set(handles.radio_visual_var_strata,'enable','off')
            set(handles.radio_visual_var_length,'enable','off')
            set(handles.radio_visual_var_lat,'enable','off')
            set(handles.radio_visual_var_lat_lon,'enable','off')
            set(handles.radio_visual_var_gender,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','off')
            set(handles.radio_visual_var_age,'enable','off')
            set(handles.radio_visual_var_len_age,'enable','off')
            %            data.visual.var1=data.final.table.trawl(:,11);
            data.visual.var1_str='Age';
        otherwise
            set(handles.radio_visual_var_length,'enable','on')
            set(handles.radio_visual_var_age,'enable','on')
            switch get(hObject,'String')
                case 'NASC'
                    set(handles.radio_visual_var_length,'enable','off')
                    set(handles.radio_visual_var_age,'enable','off')
                    set(handles.radio_visual_var_len_age,'enable','off')
                    set(handles.radio_visual_var_gender,'enable','off')
                    data.visual.var1=data.final.table.biomass(:,4);
                    data.visual.var1_str='NASC (x1000 m^2 nmi^2)';
                case 'Abundance'
                    data.visual.var1=data.final.table.biomass(:,5:7);
                    data.visual.var1_str='Abundance';
                case 'Biomass'
                    data.visual.var1=data.final.table.biomass(:,8:10);
                    data.visual.var1_str='Biomass (x 1000 mt)';
            end
    end
elseif get(hdl.vis.radio_visual_biological,'value') == 1  % using biological trawl data
    set(handles.radio_visual_var_age,'enable','on')
    set(handles.radio_visual_var_bot_depth,'enable','on')
    set(handles.radio_visual_var_layer_depth,'enable','on')
    set(handles.radio_visual_var_strata,'enable','on')
    set(handles.radio_visual_var_trawl,'enable','on')
    set(handles.radio_visual_var_lat,'enable','on')
    set(handles.radio_visual_var_transect,'enable','on')
    set(handles.radio_visual_var_gender,'enable','on')
    set(handles.radio_visual_var_weight,'enable','on')
    set(handles.radio_visual_var_survey_region,'enable','off')
    set(handles.radio_visual_var_length,'enable','on')
    set(handles.radio_visual_var_lat_lon,'enable','on')
    set(handles.radio_visual_var_len_age,'enable','off')
    set(handles.radio_visual_func_NASC,'enable','off')
    construct_catch_table
    switch get(hObject,'String')
        % construct catch table
        case 'Length'
            set(handles.radio_visual_var_length,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','on')
            data.visual.var1=data.final.table.trawl(:,9);
            data.visual.var1_str='Length (cm)';
        case 'Age'
            set(handles.radio_visual_var_age,'enable','off')
            set(handles.radio_visual_var_survey_region,'enable','on')
            data.visual.var1=data.final.table.trawl(:,11);
            data.visual.var1_str='Age';
        otherwise
            set(handles.radio_visual_var_length,'enable','on')
            set(handles.radio_visual_var_age,'enable','on')
            switch get(hObject,'String')
                case 'NASC'
                    set(handles.radio_visual_var_age,'enable','off')
                    set(handles.radio_visual_var_bot_depth,'enable','off')
                    set(handles.radio_visual_var_layer_depth,'enable','off')
                    set(handles.radio_visual_var_strata,'enable','off')
                    set(handles.radio_visual_var_trawl,'enable','off')
                    set(handles.radio_visual_var_lat,'enable','off')
                    set(handles.radio_visual_var_transect,'enable','off')
                    set(handles.radio_visual_var_gender,'enable','off')
                    set(handles.radio_visual_var_weight,'enable','off')
                    set(handles.radio_visual_var_length,'enable','off')
                    set(handles.radio_visual_var_lat_lon,'enable','of')
                    set(handles.radio_visual_var_len_age,'enable','off')
             %                     data.visual.var1=data.final.table.biomass(:,4);
             %                     data.visual.var1_str='NASC (x1000 m^2 nmi^2)';
                case 'Abundance'
                    set(handles.radio_visual_var_age,'enable','off')
                    set(handles.radio_visual_var_length,'enable','off')
                    set(handles.radio_visual_var_weight,'enable','off')
                    set(handles.radio_visual_var_len_age,'enable','on')
           %         data.visual.var1=data.final.catch;
                    data.visual.var1_str='Abundance';
                case 'Biomass'
                    set(handles.radio_visual_var_age,'enable','off')
                    set(handles.radio_visual_var_length,'enable','off')
                    set(handles.radio_visual_var_weight,'enable','off')
                    set(handles.radio_visual_var_len_age,'enable','on')
            %        data.visual.var1=data.final.table.catch;
                    data.visual.var1_str='Biomass (kg)';
            end
    end
end


