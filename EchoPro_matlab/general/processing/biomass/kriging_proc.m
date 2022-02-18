function  kriging_proc
%% perform kriging on biomass density/NASC/number density 
%% to obtain length-age-sex structured biomass and variance estimates
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global data para hdl MAP_PROJECTION MAP_VAR_LIST


%% set niminal grid_cell area corresponding to the grid_cell file
if ~isempty(findstr(data.in.filename.grid_cell,'2_5'))
    para.proc.kriging_A0=2.5*2.5;
else ~isempty(findstr(data.in.filename.grid_cell,'1_25'));
    para.proc.kriging_A0=1.25*1.25;
end
A0=para.proc.kriging_A0;   % nominal area of each grid

%% load US & CAN Mesh files off the West Coast of US/CAN 
d=xlsread(data.in.filename.grid_cell);

if para.krig.proc_opt == 1      
%% before kriging process is called
    data.in.US_CAN_Mesh=d(:,3:4);
%% whether to exclude the extrapolated region
% para.proc.extrapolation=1;
    if para.proc.extrapolation == 0
        %% find the kriged region without extrapolation
        fprintf(' Determine the kriging region that is bounded by the transect...\n')
        cmd=['[sel_grid_indx, data]=krig_transect_region(para,data);'];
        eval(cmd)
    else
        sel_grid_indx=1:size(d,1);
    end
    data.in.sel_grid_indx=sel_grid_indx;   % selected grid index of the extrapolation removed kriging grids
    US_CAN_Mesh=d(sel_grid_indx,3:4);
    Area=A0*d(sel_grid_indx,5);                % actual area of each mesh grid
    data.in.US_CAN_Mesh=US_CAN_Mesh;
    data.in.US_CAN_Mesh_area=Area;
else
%% after kriging process is called
    Area=data.in.US_CAN_Mesh_area;
end

data.in.US_CAN_Mesh_stratum=get_geographic_strata(data.in.US_CAN_Mesh(:,1),para.proc.stratification_filename,para.proc.stratification_index+1);
if para.proc.extrapolation == 0 
   acc_area = 0;
   for  i = 1:3
       if para.proc.source == i | para.proc.source == 3
           data.in.US_CAN_mesh_reg(i).area = A0*d(data.in.US_CAN_mesh_reg(i).ind,5);
           acc_area = acc_area + sum(data.in.US_CAN_mesh_reg(i).area);
       end
   end
   data = unkrig_mesh_dat(data, para); 
   data.in.US_CAN_unkrig_dat = unkrig_mesh_dat(data, para); 
   data = compare_krig_biomass(data);
end
%  disp([acc_area sum(Area)])
if para.krig.loop == 1  % kriging_loop through more flexible variables such as western bound extension with fake hake data 
    if ~exist('dat4krig')   % load the data file
        %%%%%%%%%%%%%%%%%%%% test western bound extension with fake hake data %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%    load file
        fac=sum(data.final.table.biomass(:,15))/sum(data.final.table.biomass(:,21));
        fprintf('load data file ..... \n');
        load('C:\Projects\EchoPro\EchoProGUI_Current version\input_data\2012\Western Extent Simulations\dat4krig.mat');  % load kriging-input data file
        fprintf('data has been loaded !! \n');
        [krig_i,krig_j]=size(dat4krig);
    end
    for ii=1:krig_i                % loop through realizations
        fprintf(' ********************** realization # = %d *********************************\n',ii);
        for jj=1:krig_j            % loop through extended distance (nmi)
            fprintf(' Extended %d (nmi):  \t',jj);
            Area=A0*d(:,5);                     % actual area of each mesh grid
            data0=dat4krig{ii,jj};
            
            if isfield(hdl,'krig') & isfield(hdl.krig,'msg')
                if ishandle(hdl.krig.msg)
                    delete(hdl.krig.msg)
                end
            end
            
            %% smooth line
            US_BCw_border=100;    % sample index of the boarder line
            if strmatch(para.survey_year,'1998')
                BCw_BCe_border=149;  % CAN west and east border index
            else
                BCw_BCe_border=146;  % CAN west and east border index
            end
            Lsmooth=xlsread(data.in.filename.smoothed_contour);
            Ref_lon=para.krig.shelf_break_Ref_lon; % -124.78338;
            switch para.proc.kriging_input
                case 1   % Biomass density
                    data1=data.final.table.biomass(:,[5 6 21]);
                case 2   % NASC
                    data1=data.final.table.biomass(:,[5 6 9]);
                case 3   % Number density
                    data1=data.final.table.biomass(:,[5 6 18]);
            end
            B2_unkriged=sum(data1(:,3))*fac*1e-9;   % in mmt
            data.final.table.krig.bootstrap.w(ii,jj)=B2_unkriged;  % un-kriged biomass

        %% longitude coordinate transformation: longitude vaues (function of lat) at ~200-isobath are the referrance longitudes
            [data.in.US_CAN_dat,data.in.US_CAN_mesh_cor]=convert_data(data1,Lsmooth,US_BCw_border,BCw_BCe_border,Ref_lon,data.in.US_CAN_Mesh);
            if para.proc.bootstrap.cnt > 1
                delete(hdl.krig.msg);
            end
            startkrig_biomass;
            %% kriged biomass estimate
            tic
            linear_scale=para.dataprep.transform_index;
            data.out.krig.lon=data.in.US_CAN_Mesh(:,2);
            data.out.krig.lat=data.in.US_CAN_Mesh(:,1);
            x=data.out.krig.lon;         % US west coast grids in x direction
            y=data.out.krig.lat;         % US west coast grids in y direction
            
            %% Biomass and Variance
            %  data.final.table
            [var, var1, cv_s, C0, data]=get_kriged_biomass(data,para);
            
            Cov=nansum(var1)*C0;
            CV=sqrt(var1.*C0)./var;             % kriging CV
            CV_s=cv_s;              % sample CV
            %% display
            clr_max=length(colormap);
            %% Plot kriging resluts on the US west coast grids
            ind=find(var == 0 | isnan(var) == 1);
            var(ind)=[];
            x(ind)=[];
            y(ind)=[];
            CV(ind)=[];
            Area(ind)=[];
            var_i=data.out.krig.gv;   % for kriging_input == 2 (gv is not biomass density)
            var_i(ind)=[];        
            
            
            %%%%%%%%% Biomass of the whole survey region (mmt) %%%%%%%%%%%%%%
            B2=nansum(var.*Area)*1e-6;                       % mmt: all mesh grids
            
            %%%CAN Data %%%%%%%
            ind_CAN=find(y > 48.5);
            B2_CAN=nansum(var(ind_CAN).*Area(ind_CAN))*1e-6;
            ind_CAN_strict=find(y > 48.5 & y <= 55.0);
            B2_CAN_strict=nansum(var(ind_CAN_strict).*Area(ind_CAN_strict))*1e-6;
            
            
            y_scale=60*para.dataprep.y_norm;
            x_scale=para.dataprep.x_norm*60*cos(min(data.in.var1_raw)*pi/180);
            aspect_ratio=y_scale/x_scale*para.vario.ytox_ratio;
            
            % normalized variance - CV
            Bn=nansum(var_i.*Area)*1e-9;
            CV_var=mean(Area0)*sqrt(var1.*C0)*1e-9/Bn*sqrt(length(var1));
            CVn=mean(Area)*sqrt(Cov)*1e-9/Bn;
            %% CAN CV
            Bn_CAN=nansum(var_i(ind_CAN).*Area(ind_CAN))*1e-9;
            Cov_CAN=nansum(var1(ind_CAN))*C0;
            CVn_CAN=A0*sqrt(Cov_CAN)*1e-9/Bn_CAN;           
           
            B2_kriged=B2*1e-3;
            data.final.table.krig.bootstrap.w(ii,jj+krig_j)=B2;
            fprintf(' un-kriged biomass = %10.5f\t    kriged biomass = %10.5f\n',B2_unkriged,B2_kriged);
            SD_var = SD_map(CV_var, var, data, para);

            data.final.table.kriged_biomass=[y x var.*Area*1e-6 CV_var CV_SD];
            data.final.table.kriged_biomass_description={'Latitude','Longitude','Biomass (mt)', 'krig_CV','krig_SD'};
        end  % end loop through kriging realization
    end     % end of loop through extended distance

else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% no kriging loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if para.krig.proc_opt == 1      % perform kriging
        %% smooth line
        US_BCw_border=100;    % sample index of the boarder line
        if strmatch(para.survey_year,'1998')
            BCw_BCe_border=149;  % CAN west and east border index
        else
            BCw_BCe_border=146;  % CAN west and east border index
        end
        Lsmooth=xlsread(data.in.filename.smoothed_contour);
        Ref_lon=para.krig.shelf_break_Ref_lon; % -124.78338;
        switch para.proc.kriging_input
            case 1   % Biomass density
                data1=data.final.table.biomass(:,[5 6 21]);
            case 2   % NASC
                data1=data.final.table.biomass(:,[5 6 9]);
            case 3   % Number density
                data1=data.final.table.biomass(:,[5 6 18]);
        end
        %% longitude coordinate transformation: longitude vaues (function of lat) at ~200-isobath are the referrance longitudes
%         data1=xlsread('C:\Projects\EchoPro\EchoProGUI_Current\input_data\2009\test_data\NASC_tbl_z3_nz3.xlsx');
        [data.in.US_CAN_dat,data.in.US_CAN_mesh_cor]=convert_data(data1,Lsmooth,US_BCw_border,BCw_BCe_border,Ref_lon,data.in.US_CAN_Mesh);
%%%%%%%%%%%%%% 2009 kriging data from previous version  %%%%%%%%%%%%%%%%%%%
%        data1=load('G:\Office Desktop Computer\My Documents\Projects\EK60 DataProc\Geostatistics\inputs\2009\2009aUS_CAN_merge_biom_dat_zero_cor-30-Dec-2010.txt');
%        data1=load('G:\Office Desktop Computer\My Documents\Projects\EK60 DataProc\Geostatistics\inputs\2009\v3_US_CAN_merge_biom_dat_zero.txt');
%        data.in.US_CAN_dat=data1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if para.proc.bootstrap.cnt > 1
            if isfield(hdl,'krig') & isfield(hdl.krig,'hmsg')& isfield(hdl.krig,'hmsg')
                delete(hdl.krig.hmsg);
            end
        end
        fprintf('\n Kriging ...\n')
        startkrig_biomass;
%         data.data1=data1;
%         data.data2=data2;
    else  % kriged biomass estimate
        tic
        linear_scale=para.dataprep.transform_index;
        data.out.krig.lon=data.in.US_CAN_Mesh(:,2);
        data.out.krig.lat=data.in.US_CAN_Mesh(:,1);
        x=data.out.krig.lon;         % US west coast grids in x direction
        y=data.out.krig.lat;         % US west coast grids in y direction
        
        %% Biomass and Variance
        [var, var1, cv_s, C0, data]=get_kriged_biomass(data,para);
            %% C0 = variance C(0), a scalar
            %% var = normalized kriging varaince
            %% var_s = sample variance
        Area0=data.in.US_CAN_Mesh_area;
        CV=sqrt(var1.*C0)./var;                 % kriging CV
        CV_s=cv_s;                              % sample CV
        %% display
        clr_max=length(colormap);
        n=length(x);
        ind_valid=1:n;
        %% remove bad points and zero value for display kriging maps
%         ind=find(var == 0 | isnan(var) == 1);
        ind=find(isnan(var) == 1);
        ind_valid(ind)=[];
        var(ind)=[];
        x(ind)=[];
        y(ind)=[];
        var1(ind)=[];
        CV(ind)=[];
        CV_s(ind)=[];
        Area(ind)=[];
        var_i=data.out.krig.gv;   % for kriging_input == 2 (gv is not biomass density)
        var_i(ind)=[];
        Cov=nansum(var1)*C0;                    % normalized var (dimensionless) to actual var
        
        if para.proc.bootstrap.cnt >= para.proc.bootstrap.limit
%             var_name='Biomass (mt)';
%             display_NASC_number_biomass(x,y,var,var_name,para.visual.disp_clr_fac,8)
            n=length(x);

            % log domain
%             var_disp=log(var+eps);
%             para.visual.disp_clr_fac = 1.0;
%             vmax=max(var);
            
            % linear domain
            var_disp=var;
            vmax=0.8*max(var);
            var_disp(var_disp> vmax)=vmax;
            clr=max(0,floor(clr_max*(abs(var_disp)-min(abs(var_disp)))/(para.visual.disp_clr_fac*max(abs(var_disp))-min(abs(var_disp)))))+1;
%             clr=min(clr_max, max(0,floor(clr_max*(var_disp-min(var_disp))/(para.visual.disp_clr_fac*max(var_disp)-min(var_disp))))+1);
            
            indx_nan=find(isnan(abs(var_disp)) == 1);
            var_disp(indx_nan)=0;
            clr=min(clr,clr_max);
            
            figure(8)     % color-coded kriged biomass map
            clf
            m_proj('Mercator','lat',[34 60],'long',[-140 -120]);
            m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
            cmap=colormap('jet'); colormap(cmap)
            hold on
            i0=1;
            for i=1:n
                if isempty(indx_nan) | min(abs(indx_nan-i)) ~= 0
                    if clr(i) <= 0
                        m_plot(x(i),y(i),'o','color',[1 1 1],'markersize',2);
                    else
                        m_plot(x(i),y(i),'o','color',cmap(clr(i),:),'markersize',2);
                    end
                    if i == i0;hold on,end
                else
                    i0=i0+1;
                end

            end
            xlabel('LONGITUDE','fontsize',16,'fontweight','bold')
            ylabel('LATITUDE','fontsize',16,'fontweight','bold')
            title(['BIOMASS DENSITY (' para.survey_year ' Hake Survey)'],'fontsize',14,'fontweight','bold');
            m_grid('linest','none','linewidth',2,'tickdir','in','xaxisloc','bottom');
            set(gca,'dataaspectratio',[1 2 1]);

            %   colormap(cmap)
            pos = get(gca,'Position');
            stripe = 0.075; edge = 0.02;
            [az,el] = view;
            if all([az,el]==[0 90]), space = 0.05; else space = .1; end
            set(gca,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
            rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
            ax= axes('Position', rect);
            colormap(cmap)
            image([0 1],[min(abs(var_disp)) max(abs(var_disp))],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')');
            set(ax,'Ydir','normal')
            set(ax,'YAxisLocation','right')
            set(ax,'xtick',[])
            title(ax,'tons/nm^2','fontsize',14,'fontweight','bold')
            drawnow
        end
        
        %%%%%%%%% Biomass of the whole survey region (mmt) %%%%%%%%%%%%%%
        B2=nansum(var.*Area)*1e-6;                       % kmt: all mesh grids
        
        %%%CAN Data %%%%%%%
        ind_CAN=find(y > para.proc.WA_CAN_border);
        B2_CAN=nansum(var(ind_CAN).*Area(ind_CAN))*1e-6;
        ind_CAN_strict=find(y > para.proc.WA_CAN_border & y <= para.proc.AK_CAN_border);
        B2_CAN_strict=nansum(var(ind_CAN_strict).*Area(ind_CAN_strict))*1e-6;
        
        fprintf('\n\n Survey Year  = %s\n',para.survey_year);
        fprintf('Mesh Grids: Biomass_{krig} = %6.2f (kmt)\n',B2);
        fprintf('Mesh Grids: Biomass_{krig} (CAN) = %6.2f (kmt)\n',B2_CAN);
        fprintf('Mesh Grids: Biomass_{krig} (CAN < 55 deg) = %6.2f (kmt)\n',B2_CAN_strict);
        
        disp(' ------------------------------')
        y_scale=60*para.dataprep.y_norm;
        x_scale=para.dataprep.x_norm*60*cos(min(data.in.var1_raw)*pi/180);
        aspect_ratio=y_scale/x_scale*para.vario.ytox_ratio;
        fprintf('Lat Scale = %8.4f (nmi)\n',y_scale);
        fprintf('Lon Scale = %8.4f (nmi)\n',x_scale);
        fprintf('Aspect Ratio = %8.4f \n',aspect_ratio);
        fprintf('Kriging Length Scale = %8.4f (norm)\n',para.vario.lscl);
        fprintf('Kriging Length Scale (along isobath) = %8.4f (nmi)\n',para.vario.lscl*y_scale);
        fprintf('Kriging Length Scale (across isobath) = %8.4f (nmi)\n',para.vario.lscl*x_scale);
        fprintf('Kriging Correlation Length (along isobath) = %8.4f (nmi)\n',para.vario.lscl*y_scale*log(sqrt(2)));
        fprintf('Kriging Correlation Length (across isobath) = %8.4f (nmi)\n',para.vario.lscl*x_scale*log(sqrt(2)));
        
        % normalized variance - coefficient of variance CV
        Bn=nansum(var_i.*Area)*1e-9;
        CV_var=A0*sqrt(var1.*C0)*1e-9/Bn*sqrt(length(var1));
%         CVn=A0*sqrt(Cov)*1e-9/Bn;
        % more accurate   10/26/2017
        CVn = sqrt(nansum(var1.*Area.^2)*C0)*1e-9/Bn;

        %% CAN CV
        Bn_CAN=nansum(var_i(ind_CAN).*Area(ind_CAN))*1e-9;
        Cov_CAN=nansum(var1(ind_CAN))*C0;
%         CVn_CAN=A0*sqrt(Cov_CAN)*1e-9/Bn_CAN;
        % more accurate   10/26/2017
        CVn_CAN = sqrt(nansum(var1(ind_CAN).*Area(ind_CAN).^2)*C0)*1e-9/Bn_CAN;
               
                
        fprintf('CV = %8.6f  \n',CVn)
        fprintf('CV_CAN = %8.6f  \n',CVn_CAN)
        fprintf('=======================================\n\n');
        
        if para.proc.bootstrap.cnt >= para.proc.bootstrap.limit
        %% Plot kriged CV map
%             var_name='CV';
%             x=data.out.krig.lon;
%             y=data.out.krig.lat;
%             display_NASC_number_biomass(x,y,CV_var,var_name,1)
            cv_max=max(CV_var);
            CV_disp=CV_var;
            CV_disp(var_disp> vmax)=cv_max;
            % CV_disp(end)=cv_max;
            clr=max(0,floor(clr_max*(CV_disp-min(CV_disp))/(max(CV_disp)-min(CV_disp))))+1;
            indx_nan=find(isnan(CV_disp) == 1);
            CV_disp(indx_nan)=0;
            clr=min(clr,clr_max);
            cmap=colormap;
            
            figure(9)  % color-coded kriged CV map
            clf
            m_proj('Mercator','lat',[34 60],'long',[-140 -120]);
            m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
            hold on
            m_plot(data.out.krig.lon,data.out.krig.lat,'.','markersize',1)
            i0=1;
            ncv=length(CV_var);
            for i=1:ncv
                if isempty(indx_nan) | min(abs(indx_nan-i)) ~= 0
                    m_plot(data.out.krig.lon(i),data.out.krig.lat(i),'.','color',cmap(clr(i),:),'markersize',2);
                    if i == i0;hold on,end
                else
                    i0=i0+1;
                end
            end
            xlabel('LONGITUDE','fontsize',16,'fontweight','bold')
            ylabel('LATITUDE','fontsize',16,'fontweight','bold')
            title(['COEFFICIENT OF VARIANCE (' para.survey_year ' Hake Survey)'],'fontsize',14,'fontweight','bold');
            m_grid('linest','none','linewidth',2,'tickdir','in','xaxisloc','bottom');
            set(gca,'dataaspectratio',[1 2 1]);
            
            colormap('jet')
            pos = get(gca,'Position');
            stripe = 0.075; edge = 0.02;
            [az,el] = view;
            if all([az,el]==[0 90]), space = 0.05; else space = .1; end
            set(gca,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
            rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
            ax= axes('Position', rect);
            image([0 1],[min(CV_disp) max(CV_disp)+0.03],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')');
            set(ax,'Ydir','normal')
            set(ax,'YAxisLocation','right')
            set(ax,'xtick',[])
            title(ax,'CV','fontsize',14,'fontweight','bold')
            drawnow
            SD_var = SD_map(CV_var, var, data, para);
        end
        
 %       data.final.table.kriged_biomass=[y x var.*Area*1e-6 CV_var(ind_valid)];
        data.final.table.kriged_biomass=[y x var.*Area*1e-6 CV_var SD_var];
        data.final.table.kriged_biomass_description={'Latitude','Longitude','Biomass (mt), krig_CV', 'krig_SD'};
        data.final.table.kriged_biomass0(:,12)=CV_var;
        data.final.table.kriged_biomass0(:,13)=SD_var;
        fprintf(' -------------------------------------------------------------------\n');
        
        if para.proc.bootstrap.cnt <= para.proc.bootstrap.limit
            cnt=para.proc.bootstrap.cnt;
            if para.proc.bootstrap.cnt == 1
                data.final.bootstrap.out_description={'count','no. of transects','Total Biomass','CAN Biomass','Kriging CV', 'Jolly-Hampton CV'};
            end
            Tnew=sort(unique(data.final.table.biomass(:,1)));
            n=length(Tnew);
            data.final.bootstrap.transects(cnt,1:n)=Tnew;
            CVjh = jolly_hampton_CV(para, data);
            CVjh_mean = nanmean(CVjh);
            data.final.bootstrap.out(cnt,1:6)=[cnt n B2 B2_CAN CVn CVjh_mean];
            fprintf('Jolly-Hampton (kriged): CV = %6.4f \n\n',CVjh_mean)
            fprintf('Elapse Time = %10.3f (seconds)\n',etime(clock,data.proc.start_time))
            para.proc.total_time=etime(clock,data.proc.start_time);
            if para.proc.bootstrap.cnt == para.proc.bootstrap.limit
                fprintf('*** The Total Elapse Time = %10.2f (seconds)\n\n',etime(clock,data.proc.start_time))
            end
        end
        data.final.CV_s=nanmean(CV_s);
        %% degree of coverage
        d = data.final.table.biomass(:,4) - data.final.table.biomass(:,3); % VLend - VLstart  in nmi
        D = nansum(d);          % total transect length in nmi
        A = D/sqrt(nansum(Area));   % degree of coverage
        data.final.table.degree_of_coverage = A;
        fprintf('Degree of Coverage = %4.2f\n', A);
    end
    if isfield(para.krig,'hmsg') & ishandle(para.krig.hmsg)
        close(para.krig.hmsg)
    end
end
return