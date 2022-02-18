function proc_visualization
% visualize the trawl and acoustic data & results
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global para hdl data

lat_dis_inc=15;                                    % distance for each lat bin (nmi)
nlat_inc=60/lat_dis_inc;                           % number of latitude samples per lat deg
lat_inc=1/nlat_inc;                                % latitude increment in deg for each lat bin (nmi)

if get(hdl.vis.radio_visual_un_kriged,'value') == 1  % using values prior to kriging
    para.proc.disp_method=[1 2];                     % display method
    if get(hdl.vis.radio_visual_func_NASC,'value') == 1                  %% NASC
        data.visual.func=data.final.table.biomass(:,9);
        data.visual.func_str='<NASC>';
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            data.visual.var=data.final.table.biomass(:,1);
            data.visual.var_str='Transect';
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.biomass(:,5);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.biomass(:,7);
            data.visual.var_str='Stratum Index';
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            data.visual.var=data.final.table.biomass(:,22);
            ind=find(data.visual.var ==0);
            data.visual.func(ind)=[];
            data.visual.var(ind)=[];
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            data.visual.var=data.final.table.biomass(:,8);
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='NASC (x1000)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.biomass(:,5:6);
            data.visual.var_str='lat_lon';
            data.visual.func_str='NASC';
            para.proc.disp_method=3;
        end
    elseif get(hdl.vis.radio_visual_func_number,'value') == 1            %% Abundance
        data.visual.func_str='Abundance(x1000)';
        data.visual.func=1e-3*data.final.table.biomass(:,10:12);
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            data.visual.var=data.final.table.biomass(:,1);
            data.visual.var_str='Transect';
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.biomass(:,5);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.biomass(:,7);
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_length,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,3:end),2)  ...
                    nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,3:end),2) nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,3:end),2)];
            else
                data.visual.var=1e-3*[nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end),2)  ...
                    nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end),2) nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end),2)];
            end
%             data.visual.var1=data.final.table.Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var_str='Length (cm)';
%             data.visual.func1=data.visual.func;
            data.visual.var=data.final.table.Len_Age_Matrix_AcoustM(2:end,1);
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_age,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[[0 nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,3:end))]'  ...
                    [0 nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,3:end))]' [0 nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,3:end))]'];
            else
                
                data.visual.func=1e-3*[nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end),2)'  ...
                    nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end))' nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end))'];
            end
            data.visual.var=data.final.table.Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var_str='Age';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,3:end)));  ...
                    nansum(nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,3:end))); nansum(nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,3:end),2))];
            else
                data.visual.var=1e-3*[nansum(nansum(data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end)));  ...
                    nansum(nansum(data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end))); nansum(nansum(data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end)))];
            end
            data.visual.var=[1 2 3];                  % M F non-sexed
            data.visual.var_str='Gender';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            data.visual.func=data.visual.func(:,end);
            data.visual.var=data.final.table.biomass(:,22);
            ind=find(data.visual.var ==0);
            data.visual.func(ind)=[];
            data.visual.var(ind)=[];
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            data.visual.func=data.visual.func(:,end);
            data.visual.var=data.final.table.biomass(:,8);
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Abundance(x1000)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.biomass(:,5:6);
            data.visual.var_str='lat_lon';
            para.proc.disp_method=3;
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            data.visual.var_str='len_age';
            data.visual.var1=data.final.table.Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var1_str='Length (cm)';
            data.visual.var2_str='Age';
            data.visual.func=1e-3*[data.final.table.Len_Age_Matrix_AcoustM(2:end,2:end) ; ...
                data.final.table.Len_Age_Matrix_AcoustF(2:end,2:end); data.final.table.Len_Age_Matrix_AcoustALL(2:end,2:end)];
            if para.proc.exclude_age1 == 1
                data.visual.func(:,1)=0;
            end
            para.proc.disp_method=[4 5];
        end
    elseif get(hdl.vis.radio_visual_func_biomass,'value') == 1           %% Biomass
        data.visual.func_str='Biomass(tons)';
        data.visual.func=1e-3*data.final.table.biomass(:,13:15);         % in metric tons
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            data.visual.var=data.final.table.biomass(:,1);
            data.visual.var_str='Transect';
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.biomass(:,5);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.biomass(:,7);
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_length,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,3:end),2)  ...
                    nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,3:end),2) nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end),2)];
            else
                data.visual.var=1e-3*[nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2)  ...
                    nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2) nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2)];
            end
            data.visual.var=data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var_str='Length (cm)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_age,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[[0 nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,3:end))]'  ...
                    [0 nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,3:end))]' [0 nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end))]'];
            else
                data.visual.var=1e-3*[nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2)'  ...
                    nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end))' nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end))'];
            end
            data.visual.var=data.final.table.Wgt_Len_Age_Matrix_AcoustM(1,2:end);;
            data.visual.var_str='Age';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,3:end)));  ...
                    nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,3:end))); nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end),2))];
            else
                data.visual.var=1e-3*[nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end)));  ...
                    nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end))); nansum(nansum(data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end)))];
            end
            data.visual.var=[1 2 3];                  % M F non-sexed
            data.visual.var_str='Gender';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            data.visual.func=data.visual.func(:,end);
            data.visual.var=data.final.table.biomass(:,22);
            ind=find(data.visual.var ==0);
            data.visual.func(ind)=[];
            data.visual.var(ind)=[];
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            data.visual.func=data.visual.func(:,end);
            data.visual.var=data.final.table.biomass(:,8);
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Biomass(mmt)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.biomass(:,5:6);
            data.visual.var_str='lat_lon';
            para.proc.disp_method=3;
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            data.visual.var_str='len_age';
            data.visual.var1=data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.Wgt_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var1_str='Length (cm)';
            data.visual.var2_str='Age';
            data.visual.func=1e-3*[data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,2:end) ; ...
                data.final.table.Wgt_Len_Age_Matrix_AcoustF(2:end,2:end); data.final.table.Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end)];
            if para.proc.exclude_age1 == 1
                data.visual.func(:,1)=0;
            end
            para.proc.disp_method=[4 5];
        end
    end
elseif get(hdl.vis.radio_visual_kriged,'value') == 1  %%%%%%%%%%%%%%%%%% using the kriged values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    para.proc.disp_method=[1 2];                     % display method
    if get(hdl.vis.radio_visual_func_NASC,'value') == 1                  %% NASC
        data.visual.func=data.final.table.kriged_biomass0(:,4);
        data.visual.func_str='<NASC>';
        if get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,3);
            data.visual.var_str='Stratum Index';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='NASC (x1000)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1:2);
            data.visual.var_str='lat_lon';
            data.visual.func_str='NASC';
            para.proc.disp_method=3;
        end
    elseif get(hdl.vis.radio_visual_func_number,'value') == 1            %% Abundance
        data.visual.func_str='Abundance(x1000)';
        data.visual.func=1e-3*data.final.table.kriged_biomass0(:,5:7);
        if  get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,3);
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_length,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,3:end),2)  ...
                    nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,3:end),2) nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,3:end),2)];
            else
                data.visual.var=1e-3*[nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end),2)  ...
                    nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end),2) nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end),2)];
            end
            data.visual.var=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var_str='Length (cm)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_age,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[[0 nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,3:end))]'  ...
                    [0 nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,3:end))]' [0 nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,3:end))]'];
            else
                data.visual.var=1e-3*[nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end),2)'  ...
                    nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end))' nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end))'];
            end
            data.visual.var=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var_str='Age';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,3:end)));  ...
                    nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,3:end))); nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,3:end),2))];
            else
                data.visual.var=1e-3*[nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end)));  ...
                    nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end))); nansum(nansum(data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end)))];
            end
            data.visual.var=[1 2 3];                  % M F non-sexed
            data.visual.var_str='Gender';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Abundance(x1000)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1:2);
            data.visual.var_str='lat_lon';
            para.proc.disp_method=3;
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            data.visual.var1=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var_str='len_age';
            data.visual.var1_str='Length (cm)';
            data.visual.var2_str='Age';
            data.visual.func=1e-3*[data.final.table.kriged_Num_Len_Age_Matrix_AcoustM(2:end,2:end) ; ...
                data.final.table.kriged_Num_Len_Age_Matrix_AcoustF(2:end,2:end); data.final.table.kriged_Num_Len_Age_Matrix_AcoustALL(2:end,2:end)];
            if para.proc.exclude_age1 == 1
                data.visual.func(:,1)=0;
            end
            para.proc.disp_method=[4 5];
        end
    elseif get(hdl.vis.radio_visual_func_biomass,'value') == 1           %% Biomass
        data.visual.func_str='Biomass(mmt)';
        data.visual.func=1e-3*data.final.table.kriged_biomass0(:,8:10);
        if  get(hdl.vis.radio_visual_var_lat,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1);
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,3);
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_length,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,3:end),2)  ...
                    nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,3:end),2) nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end),2)];
            else
                data.visual.var=1e-3*[nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end),2)  ...
                    nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end),2) nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end),2)];
            end
            data.visual.var=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var_str='Length (cm)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_age,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=[[0 nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,3:end))]'  ...
                    [0 nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,3:end))]' [0 nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end))]'];
            else
                data.visual.func=[nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end))'  ...
                    nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end))' nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end))'];
            end
            data.visual.var=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var_str='Age';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            if para.proc.exclude_age1 == 1
                data.visual.func=1e-3*[nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,3:end)));  ...
                    nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,3:end))); nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,3:end),2))];
            else
                data.visual.var=1e-3*[nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end)));  ...
                    nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end))); nansum(nansum(data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end)))];
            end
            data.visual.var=[1 2 3];                  % M F non-sexed
            data.visual.var_str='Gender';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Biomass(mmt)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            data.visual.var=data.final.table.kriged_biomass0(:,1:2);
            data.visual.var_str='lat_lon';
            para.proc.disp_method=3;
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            data.visual.var1=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.var_str='len_age';
            data.visual.var1_str='Length (cm)';
            data.visual.var2_str='Age';
            data.visual.func=1e-3*[data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,2:end) ; ...
                data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustF(2:end,2:end); data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustALL(2:end,2:end)];
            if para.proc.exclude_age1 == 1
                data.visual.func(:,1)=0;
            end
            para.proc.disp_method=[4 5];
        end
    end
elseif get(hdl.vis.radio_visual_biological,'value') == 1  %%%%******************* biological trawl data ***************************
    if get(hdl.vis.radio_visual_func_length,'value') == 1                  %%%%%%% Length %%%%%%%
        data.visual.func_str='Length (cm)';
        para.proc.disp_method=6;                     % display method
        data.visual.func=data.final.table.trawl(:,9);
        bins=para.bio.hake_len_bin;
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            var0=data.final.table.trawl(:,2);
            ind=find(~isnan(var0) == 1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Transect';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            var0=data.final.table.trawl(:,3);
            ind=find(~isnan(var0) == 1);
            var=var0(ind);
            func=data.visual.func(ind);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
%             data.visual.var3=sort(unique(round(var*nlat_inc)/nlat_inc));
            data.visual.var3=min(round(var)):lat_inc:max(round(var));
            %% for 3D histogram
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            % for boxplot
            var_sel=round(var*nlat_inc)/nlat_inc;
            for i=1:length(data.visual.var3);
%                 ind3=find(var_sel == data.visual.var3(i));
                ind3=find(var_sel == data.visual.var3(i));
                if isempty(ind3)
                   data.visual.func3(1,i)=nan;
                else
                   data.visual.func3(1:length(ind3),i)=func(ind3);
                end
            end
            data.visual.var_str='Latitude';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            var0=data.final.table.trawl(:,5);
            ind=find(~isnan(var0) == 1);
            var=var0(ind);
            func=data.visual.func(ind);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
           end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Strata';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_trawl,'value') == 1
            var0=data.final.table.trawl(:,1);
            ind=find(~isnan(var0) == 1);
            var=var0(ind);
            func=data.visual.func(ind);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Trawl';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_age,'value') == 1
            var0=data.final.table.trawl(:,11);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var)
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Age';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            var0=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                if i == length(data.visual.var)
                    ind=find(var == 1 | var == 2);
                else
                    ind=find(var == data.visual.var(i));
                end
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Gender';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_weight,'value') == 1
            data.visual.func=data.bio.len_wgt_ALL(:);
            data.visual.var=para.bio.hake_len_bin(:);
            para.proc.disp_method=1;                % scatter plot and curve
            data.visual.var_str='Length (cm)';
            data.visual.func_str='Weight (kg)';
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            var0=round(data.final.table.trawl(:,13));
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            var0=round(data.final.table.trawl(:,6));
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Length (cm)';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            para.proc.disp_method=7;                % bubble plot
            trawl_no=data.final.table.trawl(:,1);
            biomass=data.final.table.trawl(:,14);
            len=data.visual.func;
            unique_trawl=unique(trawl_no);
            data.visual.func=[];
            data.visual.var=[];
            for i=1:length(unique_trawl);
                ind=find(trawl_no == unique_trawl(i));
                data.visual.func(i,1)=nanmean(len(ind));
                data.visual.func_biomass(i)=nansum(biomass(ind));
                data.visual.var(i,1)=nanmean(data.final.table.trawl(ind,3));
                data.visual.var(i,2)=nanmean(data.final.table.trawl(ind,4));
            end
            data.visual.var_str='lat_lon';
        end
    elseif get(hdl.vis.radio_visual_func_age,'value') == 1                  %%%%%% Age  %%%%%%%
        data.visual.func_str='Age';
        para.proc.disp_method=6;                     % display method
        data.visual.func=data.final.table.trawl(:,11);
        bins=para.bio.hake_age_bin;
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            var0=data.final.table.trawl(:,2);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Transect';
            para.proc.disp_method=[para.proc.disp_method 9];
       elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            var0=data.final.table.trawl(:,3);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
 
            data.visual.var3=min(round(var)):lat_inc:max(round(var));
            %% for 3D histogram
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            % for boxplot
            var_sel=round(var*nlat_inc)/nlat_inc;
            for i=1:length(data.visual.var3);
%                 ind3=find(var_sel == data.visual.var3(i));
                ind3=find(var_sel == data.visual.var3(i));
                if isempty(ind3)
                   data.visual.func3(1,i)=nan;
                else
                   data.visual.func3(1:length(ind3),i)=func(ind3);
                end
            end

            
%             
%             data.visual.var3=sort(unique(round(var)));
%             for i=1:length(data.visual.var);
%                 ind=find(var == data.visual.var(i));
%                 data.visual.func(i,bins)=hist(func(ind),bins);
%             end
%             var_sel=round(var);
%             for i=1:length(data.visual.var3);
%                 ind3=find(var_sel == data.visual.var3(i));
%                 data.visual.func3(1:length(ind3),i)=func(ind3);
%             end            
            
            
            data.visual.var_str='Latitude';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            var0=data.final.table.trawl(:,5);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Strata';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_trawl,'value') == 1
            var0=data.final.table.trawl(:,1);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.func3=[];
            data.visual.var3=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
                data.visual.func3(1:length(ind),i)=func(ind);
            end
            data.visual.var3=data.visual.var;
            data.visual.var_str='Trawl';
            para.proc.disp_method=[para.proc.disp_method 9];
        elseif get(hdl.vis.radio_visual_var_length,'value') == 1  % length vs age
            para.proc.disp_method=[para.proc.disp_method  8];  % add a scatter plot: length vs. age
            var0=data.final.table.trawl(:,9);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            data.visual.func3=data.final.table.trawl(ind,11);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);           
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Length (cm)';
            data.visual.var3=var;
            data.visual.func3_str=data.visual.func_str;
            data.visual.var3_str=data.visual.var_str;
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            var0=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Gender';
        elseif get(hdl.vis.radio_visual_var_weight,'value') == 1  
            var0=data.final.table.trawl(:,11);
            func0=data.final.table.trawl(:,12);
            ind=find(~isnan(var0)==1 & ~isnan(func0)==1);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=[];
            data.visual.func=zeros(length(data.visual.var),1);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i)=nanmedian(func(ind));
            end
            para.proc.disp_method=1;
            data.visual.var_str='Age';
            data.visual.func_str='Weight (kg)';
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            var0=round(data.final.table.trawl(:,13));
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            var0=round(data.final.table.trawl(:,6));
            ind=find(~isnan(var0)==1);
            var=var0(ind);
            func=data.visual.func(ind,:);
            data.visual.func=[];
            data.visual.var=unique(var);
            for i=1:length(data.visual.var);
                ind=find(var == data.visual.var(i));
                data.visual.func(i,bins)=hist(func(ind),bins);
            end
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1
            data.visual.var=nan;
            data.visual.var_str='Age';
            para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            para.proc.disp_method=7;                % bubble plot
            trawl_no=data.final.table.trawl(:,1);
            biomass=data.final.table.trawl(:,12);
            age=data.visual.func;
            unique_trawl=unique(trawl_no);
            data.visual.func=[];
            data.visual.var=[];
            for i=1:length(unique_trawl);
                ind=find(trawl_no == unique_trawl(i));
                data.visual.func(i,1)=nanmean(age(ind));
                data.visual.func_biomass(i)=nansum(biomass(ind));
                data.visual.var(i,1)=nanmean(data.final.table.trawl(ind,3));
                data.visual.var(i,2)=nanmean(data.final.table.trawl(ind,4));
            end
            data.visual.var_str='lat_lon';
        end
    elseif get(hdl.vis.radio_visual_func_number,'value') == 1           %%%%%%% Abundance %%%%%%
        data.visual.func_str='Abundance';
        trawl0_ref=data.final.table.trawl(:,1);
        para.proc.disp_method=[2];                % histogram
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            var0=data.final.table.trawl(:,2);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Transect';
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            var0=data.final.table.trawl(:,3);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            var0=data.final.table.trawl(:,5);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);     
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_trawl,'value') == 1
            var0=data.final.table.trawl(:,1);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [~,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Trawl';
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            var0=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var)+1,1);           
            for j=1:length(data.visual.var)-1
                %% individually measured samples
                ind1=find(var == data.visual.var(j));  %Male or Female
                data.visual.func(j)=length(ind1);
            end
            %% total weight-based measures
           data.visual.func(3)=nansum(data.final.table.catch(:,2));
           data.visual.var_str='Gender';
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            var0=round(data.final.table.trawl(:,13));
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            var0=round(data.final.table.trawl(:,6));
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1 & var0 ~= 0);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=length(ind1);
                data.visual.func(j,2)=length(ind2);
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,2));
            end
            data.visual.var_str='Bottom Depth (m)';
%         elseif get(hdl.vis.radio_visual_var_survey_region,'value') == 1            
%             data.visual.var=nan;
%             data.visual.var_str='Abundance';
%             para.proc.disp_method=2;
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            para.proc.disp_method=7;                % bubble plot
            trawl_no=data.final.table.trawl(:,1);
            unique_trawl=unique(trawl_no);
            data.visual.func=[];
            data.visual.var=[];
            for i=1:length(unique_trawl);
                ind=find(trawl_no == unique_trawl(i));
                %% total weight-based measures
                ind1=find(unique_trawl(i) == data.final.table.catch(:,1));
                data.visual.func(i,1)=data.final.table.catch(ind1,2);
                data.visual.var(i,1)=nanmean(data.final.table.trawl(ind,3));
                data.visual.var(i,2)=nanmean(data.final.table.trawl(ind,4));
            end
            data.visual.var_str='lat_lon';
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            para.proc.disp_method=6;                % 3D-histogram
            func0=data.final.table.trawl(:,11);
            len0=round(data.final.table.trawl(:,9));
            age0=round(data.final.table.trawl(:,11));
            ind=find(~isnan(len0) == 1 & ~isnan(age0) == 1);
            len=len0(ind);
            age=age0(ind);
            func=func0(ind);
            data.visual.func=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
            for i=1:length(para.bio.hake_age_bin)
                for j=1:length(para.bio.hake_len_bin)
                   if j == 1 
                       ind=find(age == para.bio.hake_age_bin(i) & len <= para.bio.hake_len_bin(j));
                   else
                       ind=find(age == para.bio.hake_age_bin(i) & len > para.bio.hake_len_bin(j-1) & len <= para.bio.hake_len_bin(j) );
                   end
                   data.visual.func(j,i)=nansum(func(ind));
                end
            end
            data.visual.var0_str='len_age';
            data.visual.func0_str=data.visual.func_str;
            data.visual.var=para.bio.hake_len_bin;           
            data.visual.var1=data.final.table.Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.Wgt_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.func1=data.visual.func;
            data.visual.var_str='Length (cm)';
            data.visual.func_str='Age';   
            
            
            
        end
    elseif get(hdl.vis.radio_visual_func_biomass,'value') == 1           %%%%%%% Biomass %%%%%%%
        data.visual.func=data.final.table.trawl(:,12);
        data.visual.func_str='Biomass(kg)';
        trawl0_ref=data.final.table.trawl(:,1);
        para.proc.disp_method=[2];                % bubble plot
        if get(hdl.vis.radio_visual_var_transect,'value') == 1
            var0=data.final.table.trawl(:,2);
            var0_gender=data.final.table.trawl(:,10);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Transect';
        elseif get(hdl.vis.radio_visual_var_lat,'value') == 1
            var0=data.final.table.trawl(:,3);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Latitude';
        elseif get(hdl.vis.radio_visual_var_strata,'value') == 1
            var0=data.final.table.trawl(:,5);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Strata';
        elseif get(hdl.vis.radio_visual_var_trawl,'value') == 1
            var0=data.final.table.trawl(:,1);
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Trawl';
        elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
            var0=data.final.table.trawl(:,10);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var)+1,1);           
            func=func0(ind);
            for j=1:length(data.visual.var)-1
                %% individually measured samples
                ind1=find(var == data.visual.var(j));  %Male or Female
                data.visual.func(j)=nansum(func(ind1));
            end
            %% total weight-based measures
           data.visual.func(3)=nansum(data.final.table.catch(:,3));
           data.visual.var_str='Gender';
        elseif get(hdl.vis.radio_visual_var_layer_depth,'value') == 1
            var0=round(data.final.table.trawl(:,13));
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Layer Depth (m)';
        elseif get(hdl.vis.radio_visual_var_bot_depth,'value') == 1
            var0=round(data.final.table.trawl(:,6));
            var0_gender=data.final.table.trawl(:,10);
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            func0=data.visual.func;
            ind=find(~isnan(var0)==1);
            var_gender=var0_gender(ind);
            trawl_ref=trawl0_ref(ind);
            var=var0(ind);
            func=func0(ind);
            data.visual.var=unique(var);
            data.visual.func=zeros(length(data.visual.var),3);           
            for j=1:length(data.visual.var);
                %% individually measured samples
                ind1=find(var == data.visual.var(j) & var_gender == 1);  %Male
                ind2=find(var == data.visual.var(j) & var_gender == 2);  %FeMale
                data.visual.func(j,1)=nansum(func(ind1));
                data.visual.func(j,2)=nansum(func(ind2));
                %% total weight-based measures
                ind0=find(var == data.visual.var(j));
                [Tnum,ind]=intersect(data.final.table.catch(:,1),unique(trawl_ref(ind0)));
                data.visual.func(j,3)=nansum(data.final.table.catch(ind,3));
            end
            data.visual.var_str='Bottom Depth (m)';
        elseif get(hdl.vis.radio_visual_var_lat_lon,'value') == 1
            para.proc.disp_method=7;                % bubble plot
            trawl0_no=data.final.table.trawl(:,1);
            lat0=data.final.table.trawl(:,3);
            lon0=data.final.table.trawl(:,4);
            ind=find(~isnan(trawl0_no) == 1 & ~isnan(lat0) == 1 & ~isnan(lon0) == 1);
            trawl_no=trawl0_no(ind);
            lat=lat0(ind);
            lon=lon0(ind);
            unique_trawl=unique(trawl_no);
            data.visual.func=[];
            data.visual.var=[];
            for i=1:length(unique_trawl);
                ind=find(trawl_no == unique_trawl(i));
                %% total weight-based measures
                ind1=find(unique_trawl(i) == data.final.table.catch(:,1));
                data.visual.func(i,1)=data.final.table.catch(ind1,3);
                data.visual.var(i,1)=nanmean(lat(ind));
                data.visual.var(i,2)=nanmean(lon(ind));
            end
            data.visual.var_str='lat_lon';
        elseif get(hdl.vis.radio_visual_var_len_age,'value') == 1
            para.proc.disp_method=6;                % 3D-histogram
            func0=data.final.table.trawl(:,12);
            len0=round(data.final.table.trawl(:,9));
            age0=round(data.final.table.trawl(:,11));
            ind=find(~isnan(len0) == 1 & ~isnan(age0) == 1);
            len=len0(ind);
            age=age0(ind);
            func=func0(ind);
            data.visual.func=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin));
            for i=1:length(para.bio.hake_age_bin)
                for j=1:length(para.bio.hake_len_bin)
                   if j == 1 
                       ind=find(age == para.bio.hake_age_bin(i) & len <= para.bio.hake_len_bin(j));
                   else
                       ind=find(age == para.bio.hake_age_bin(i) & len > para.bio.hake_len_bin(j-1) & len <= para.bio.hake_len_bin(j) );
                   end
                   data.visual.func(j,i)=nansum(func(ind));
                end
            end
            data.visual.var0_str='len_age';
            data.visual.func0_str=data.visual.func_str;
            data.visual.var=para.bio.hake_len_bin;           
            data.visual.var1=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(2:end,1);
            data.visual.var2=data.final.table.kriged_Wgt_Len_Age_Matrix_AcoustM(1,2:end);
            data.visual.func1=data.visual.func;
            data.visual.var_str='Length (cm)';
            data.visual.func_str='Age';   
        end
    end  % end of Function Option
end   % end of Data Type Option

plot_results;

return

