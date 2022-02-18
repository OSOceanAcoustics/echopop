function   get_operation_window(hObject,handles)
%  get operation window depending on the operation selection
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para data hdl

%% get survey year if survey year has not beem defined or processed
if ~isfield(para,'survey_year')   
    h_survey_year=findobj(hdl.main,'tag','popup_survey_year');
    available_survey_years=get(h_survey_year,'string');
    ind=get(h_survey_year,'value');
    tmp=char(available_survey_years(ind));
    para.survey_year=tmp(1:4);
end

if ~isfield(para.proc,'KS_stratification')
    load_proc_parameters(para.survey_year);
end

%% re-load filenames & parameters defined in the input_file 
%% based on the updated survey year

Tag_name=get(hObject,'Tag');
switch Tag_name
    case 'pb_redefine_filenames'  % load all necessary files
        varargout_hdl=redefine_filenames;
        hdl.redefine_filenames=varargout_hdl;
    case 'pb_processing' % process acoustic data
        varargout_hdl=processing;
        hdl.processing=varargout_hdl;
        if str2num(para.survey_year) < 2003
%              set(findobj(varargout_hdl,'Tag','popup_proc_scale'),'value',1);
             set(findobj(varargout_hdl,'Tag','popup_proc_scale'),'enable','off');
%              set(findobj(varargout_hdl,'Tag','popup_proc_stratification'),'value',1);
%              set(findobj(varargout_hdl,'Tag','popup_proc_stratification'),'enable','off');
             set(findobj(varargout_hdl,'Tag','popup_proc_variable'),'value',para.proc.kriging_input);
%              set(findobj(varargout_hdl,'Tag','radio_proc_mix_region'),'enable','off');
%              set(findobj(varargout_hdl,'Tag','edit_start_transect_num'),'enable','off');
%              set(findobj(varargout_hdl,'Tag','edit_end_transect_num'),'enable','off');
             set(findobj(varargout_hdl,'Tag','pb_proc_advanced'),'enable','off');
             switch para.proc.source
                 case 1
                     set(findobj(varargout_hdl,'Tag','radio_US'),'value',1);
                 case 2
                     set(findobj(varargout_hdl,'Tag','radio_CAN'),'value',1);
                 case 3
                     set(findobj(varargout_hdl,'Tag','radio_ALL'),'value',1);
             end
            set(findobj(varargout_hdl,'Tag','edit_start_transect_num'),'string',para.proc.start_transect);
            set(findobj(varargout_hdl,'Tag','edit_end_transect_num'),'string',para.proc.end_transect);
        else
            set(findobj(varargout_hdl,'Tag','popup_proc_variable'),'value',para.proc.kriging_input);
            set(findobj(varargout_hdl,'Tag','radio_proc_mix_region'),'enable','off');
            set(findobj(varargout_hdl,'Tag','popup_proc_scale'),'enable','off');
            switch para.proc.source
                case 1
                    set(findobj(varargout_hdl,'Tag','radio_US'),'value',1);
                case 2
                    set(findobj(varargout_hdl,'Tag','radio_CAN'),'value',1);
                case 3
                    set(findobj(varargout_hdl,'Tag','radio_ALL'),'value',1);
            end
%             set(findobj(varargout_hdl,'Tag','popup_proc_stratification'),'value',para.proc.stratification_index+1);
            set(findobj(varargout_hdl,'Tag','edit_start_transect_num'),'string',para.proc.start_transect);
            set(findobj(varargout_hdl,'Tag','edit_end_transect_num'),'string',para.proc.end_transect);
       end
    case 'pb_show_parameter'   % show parameters
        hdl.show_parameter=show_parameters;
    case 'pb_visualization'  % results visualization
        hdl.visualization=visualization;
        set(findobj(hdl.visualization,'Tag','radio_visual_kriged'),'value',para.visual.kriged)
    case 'pb_report'  % generate reports and/or files for Oracle database
        hdl.report=report;
end
hdl.main_handles=handles;  
return