function  [indx_uniq, data]=krig_transect_region(para,data)
%% determine the region of kriged grids that emncompass only the transect region with no extrapolation

plot_flag=0;                          % display intermidate grid/transect position flag 

%% lat/lon of the original kriging grids
lat_k=data.in.US_CAN_Mesh(:,1);
lon_k=data.in.US_CAN_Mesh(:,2);
%% lat/lon of the sample points on transects
lat_t=data.final.table.biomass(:,5);
lon_t=data.final.table.biomass(:,6);
%% change bad lat & lon data (999) to interpolated values
ind=find(lat_t > 90 | lon_t > 180);
if ~isempty(ind)
    indx=(min(ind)-1):(max(ind)+1);
    for i=1:length(ind)
        indx(indx == ind(i))=[];
    end
    lat_t(ind)=interp1(indx,lat_t(indx),ind);
    lon_t(ind)=interp1(indx,lon_t(indx),ind);
end

% transect number
tx=data.final.table.biomass(:,1);
uniq_tx=unique(tx);             % unique transect numbers
nt=length(uniq_tx);             % number of unique transects

% lat_b1=51.94;                   % display starting latitude of the lower bound of QC Island and north 
lat_b1=49;                        % display starting latitude of the lower bound of QC Island and north
lat_b1=34.3;                      % display starting latitude of all survey region

if plot_flag == 1
    figure(1)
    clf
end
% find the lat/lon bounds of each transect 
legend_flag=1;
i0=1;
for i=i0:nt
%     if uniq_tx(i) == 88
%         disp(i)
%     end
    ind=find(tx == uniq_tx(i));
    lat=lat_t(ind);
    lon=lon_t(ind);
    if max(lat)-min(lat) < max(lon)-min(lon)   
        %% more likely E-W transects, find E-W longitude first
        [val,indx1]=min(lon);     
        lon1(i)=val;
        lat1(i)=lat(indx1);
        [val,indx2]=max(lon);
        lon2(i)=val;
        lat2(i)=lat(indx2);
    else
        %% more likely S-N transects, find S-N latitude first
        [val,indx1]=min(lat);
        lat1(i)=val;
        lon1(i)=lon(indx1);
        [val,indx2]=max(lat);
        lat2(i)=val;
        lon2(i)=lon(indx2);
    end
    if plot_flag == 1  % plot one transect at a time
        if mean(lat) >= lat_b1  % north of SC/QC IS
            plot(lon,lat,'.-b',lon1(i),lat1(i),'^m',lon2(i),lat2(i),'og',lon(1),lat(1),'sr','linewidth',2)
            title(sprintf('Transect = %d',uniq_tx(i)))
            if rem(uniq_tx(i),1) == 0
                text(mean(lon),mean(lat)+0.02,sprintf('%d',uniq_tx(i)),'color','k','fontweight','bold');
            else
                text(mean(lon),mean(lat)+0.02,sprintf('%3.1f',uniq_tx(i)),'color','k','fontweight','bold');
            end
%             if legend_flag == 0
%                 legend('Transect','Transect SW Bound','Transect NE Bound','Transect Start Pos');
%                 legend_flag=1;
%             end
            if legend_flag == 1 & i == nt
                legend('Transect','Transect SW Bound','Transect NE Bound','Transect Start Pos');
            end
%             if uniq_tx(i) >= 1240 & uniq_tx(i) <= 1380
%                 axis([-134.5 -129 51 55])
%                 grid on
%                 disp(i)
%             end
        end
        if  i== i0
            hold on
        end
    end
end
if plot_flag == 1
    grid on
end

dlat=1.25/60;               % resolution in latitude direction, half of the kriging grid cell

%% finding kriging grids in different regions

if para.proc.source ~= 2
    % ===========================================================================================
    %% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    %% get transect region definitions for region 1 of para.survey_year
    eval(['[tx0,tx1,tx_l,tx_r]=transect_region_def_' para.survey_year '(1);'])
    % 2011 example
    % tx0=1;    % southern most transect number
    % tx1=118;  % northern most transect number
    % tx_l=[tx0:tx1]+0.1;
    % tx_r=[tx0:tx1]+0.4;
    
    % find upper and lower bounds of latitudes for the region
    ind0=find(tx == tx0);
    ind1=find(tx == tx1);
    lat_1a=min(lat_t(ind0));  % minimum latitude for region 1
    lat_1b=max(lat_t(ind1));  % maximum latitude for region 1
    
    
    %% west bound (left side) of each transect
    k=1;
    for i=1:length(tx_l)
        ind0=find(floor(tx) == floor(tx_l(i)));
        if ~isempty(ind0)
            ind1=round(rem(tx_l(i),1)*10);  % check is 0.1, 0.4, 0.6, or 0.9?
            switch ind1
                case 1  % west end of transect
                    [val, ind]=min(lon_t(ind0));
                case 4  % east end of transect
                    [val, ind]=max(lon_t(ind0));
                case 6  % south end of transect
                    [val, ind]=min(lat_t(ind0));
                    val=lon_t(ind0(ind));
                case 9  % north end of transect
                    [val, ind]=max(lat_t(ind0));
                    val=lon_t(ind0(ind));
            end
            lon_w(k)=val;
            lat_w(k)=lat_t(ind0(ind));
            %% overwrite kth transect if it has a duplicated latitude
            if k> 1 & lat_w(k) == lat_w(k-1)
                %             disp([i k floor(tx_l(i)) lat_w(k)])
            else
                k=k+1;
            end
        end
    end
    lat_1=max(lat_1a, min(lat_w)):dlat:min(lat_1b, max(lat_w)); % latitude array for region 1
    
    %% region #1 left side (west) bound in longitude
    % sort lat_w array
    [var, ind_sort]=sort(lat_w);
    lon_l1=interp1(lat_w(ind_sort),lon_w(ind_sort),lat_1);
    
    %% east bound (right side) of each transect
    k=1;
    for i=1:length(tx_r)
        ind0=find(floor(tx) == floor(tx_r(i)));
        if ~isempty(ind0)
            ind1=round(rem(tx_r(i),1)*10);
            switch ind1
                case 1  % west end of transect
                    [val, ind]=min(lon_t(ind0));
                case 4  % east end of transect
                    [val, ind]=max(lon_t(ind0));
                case 6  % south end of transect
                    [val, ind]=min(lat_t(ind0));
                    val=lon_t(ind0(ind));
                case 9  % north end of transect
                    [val, ind]=max(lat_t(ind0));
                    val=lon_t(ind0(ind));
            end
            lon_e(k)=val;
            lat_e(k)=lat_t(ind0(ind));
            %% overwrite kth transect if it has a duplicated latitude
            if k> 1 & lat_e(k) == lat_e(k-1)
                %               disp([i k floor(tx_r(i)) lat_e(k)])
            else
                k=k+1;
            end
        end
        if i >= 4000
            figure(1)
            hold on
            plot(lon_e(i),lat_e(i),'sk','markersize',12)
            hold off
            disp(i)
        end
    end
    
    lat_1=max(lat_1a, min(lat_e)):dlat:min(lat_1b, max(lat_e)); % latitude array for region 1
    
    %% region #1 right side (east) bound in longitude
    lon_r1=interp1(lat_e,lon_e,lat_1);
    
    % figure(1);hold on;plot(lon_r1,lat_1,'g',lon_l1,lat_1,'g');hold off
    
    %% find the indecies of the krig grids bounded by the transects from S to N
    indx_i=[];
    for i=1:length(lat_1)
        dlon(i)=dlat*cosd(lat_1(i));    % resolution in longitude direction
        ind=find(lon_k >= lon_l1(i)-dlon(i) & lon_k <= lon_r1(i)+dlon(i) ...
            & lat_k >= lat_1(i)-dlat & lat_k < lat_1(i)+dlat);
        %     figure(2);plot(lon_l1(i)+[-dlon(i) dlon(i)], lat_1(i)+[-dlat dlat], 'sb', lon_k, lat_k, '+r')
        if ~isempty(ind)
            indx_i=[indx_i; ind];      % concatenated selected krig grid index array
        end
    end
    
    indx_i_uniq=unique(indx_i);
    indx = indx_i_uniq;
    
    data.in.US_CAN_mesh_reg(1).ind = indx_i_uniq;
    data.in.US_CAN_mesh_reg(1).lat = lat_k(indx_i_uniq);
    data.in.US_CAN_mesh_reg(1).lon = lon_k(indx_i_uniq);
    data.in.US_CAN_mesh_reg(1).tx_min = tx0;
    data.in.US_CAN_mesh_reg(1).tx_max = tx1;
    
end   % end of "if para.proc.source ~= 2"

if para.proc.source ~= 1    
    % ===========================================================================================
    %% region 2: transects paralell to longitudes north of QCI
    %% get transect region definitions for region 2 of para.survey_year
    eval(['[tx0,tx1,tx_l,tx_u]=transect_region_def_' para.survey_year '(2);'])
    if ~isempty(tx0)
        
        % 2011 example
        % tx0=118;  % east most transect number
        % tx1=124;  % west most transect number
        % tx_l=[114.4 114.1 118.1 120.6 122.6 124.6];
        % tx_u=[114.4 118.4 120.9 122.9 124.9];
        
        % find upper and lower bounds of longitude for the region
        ind0=find(tx == tx0);
        ind1=find(tx == tx1);
        lon_2a=min(lon_t(ind1));    % minimum longitude for region 2
        lon_2b=max(lon_t(ind0));    % maxmum longitude for region 2
        lat_m=mean(lat_t(ind1));    % mean latitude for region 2
        dlon=dlat*cosd(lat_m);      % scaled logitude increment
        
        lon_2=lon_2a:dlon:lon_2b;   % longitude array for region 2

        
        %% specifies lower (south) and upper (north) region boundaries based on the transects
        %% #.1 = west end of transect
        %% #.4 = east end of transect
        %% #.6 = south end of transect
        %% #.9 = north end of transect
        % 2011 example
        % tx_l=[114.4 114.1 118.1 120.6 122.6 124.6];
        % tx_u=[114.4 118.4 120.9 122.9 124.9];
        %% lat/lon of the south (lower) bound
        k = 1;
        for i=1:length(tx_l)
            ind0=find(floor(tx) == floor(tx_l(i)));
            if floor(tx_l(i)) == 122
                disp(ind0)
            end
            if ~isempty(ind0)
                ind1=round(rem(tx_l(i),1)*10);
                switch ind1
                    case 1   % west end of transect
                        [val, ind]=min(lon_t(ind0));
                        val=lat_t(ind0(ind));
                    case 4   % east end of transect
                        [val, ind]=max(lon_t(ind0));
                        val=lat_t(ind0(ind));
                    case 6   % south end of transect
                        [val, ind]=min(lat_t(ind0));
                    case 9   % north end of transect
                        [val, ind]=max(lat_t(ind0));
                end
                lat_s(k)=val;
                lon_s(k)=lon_t(ind0(ind));
                if k> 1 & lon_s(k) == lon_s(k-1)
                    %               disp([i k floor(tx_l(i)) lon_s(k)])
                else
                    k=k+1;
                end
            end
        end
        
        lon_2=max(lon_2a, min(lon_s)):dlon:min(lon_2b, max(lon_s));   % longitude array for region 2
        %% region #2 lower (south) bound in latitude
        
        lat_l2=interp1(lon_s,lat_s,lon_2);
        
        %% lat/lon of the north (upper) bound
        k = 1;
        for i=1:length(tx_u)
            ind0=find(floor(tx) == floor(tx_u(i)));
            if ~isempty(ind0)
                ind1=round(rem(tx_u(i),1)*10);
                switch ind1
                    case 1  % west end of transect
                        [val, ind]=min(lon_t(ind0));
                        val=lat_t(ind0(ind));
                    case 4  % east end of transect
                        [val, ind]=max(lon_t(ind0));
                        val=lat_t(ind0(ind));
                    case 6  % south end of transect
                        [val, ind]=min(lat_t(ind0));
                    case 9  % north end of transect
                        [val, ind]=max(lat_t(ind0));
                end
                lat_n(k)=val;
                lon_n(k)=lon_t(ind0(ind));
                if k> 1 & lon_n(k) == lon_n(k-1)
                    %               disp([i k floor(tx_u(i)) lon_n(k)])
                else
                    k=k+1;
                end
            end
        end
        lon_2=max(lon_2a, min(lon_n)):dlon:min(lon_2b, max(lon_n));   % longitude array for region 2
        
        %% region #2 upper (north) bound in latitude
        lat_u2=interp1(lon_n,lat_n,lon_2);
        
        %% find the indecies of the krig grids bounded by the transects from E to W
        indx_i = [];
        for i=1:length(lon_2)
            cond_flag=0;
            if ~isnan(lat_l2(i)) | ~isnan(lat_u2(i))
                if isnan(lat_l2(i))          % longitude at the lower end (lower latitude) is lower than that at the upper end (higher latitude)
                    cond_flag=1;
                    %             [var_max_n,max_ind_n]=max(lat_n);
                    %             [var_max_s,max_ind_s]=max(lat_s);
                    %             slp=(var_max_n-var_max_s)/(lon_n(max_ind_n)-lon_s(max_ind_s));
                    %             lat_l2_i=slp*(lon_2(i)-var_max_s)+lat_s(max_ind_s);
                    
                    [var_n2, lon2_n]=min(abs(lon_2(i)-lon_n));
                    [var_s2, lon2_s]=min(abs(lon_2(i)-lon_s));
                    slp=(lat_n(lon2_n)-lat_s(lon2_s))/(lon_n(lon2_n)-lon_s(lon2_s));
                    lat_l2_i=slp*(lon_2(i)-lon_s(lon2_s))+lat_s(lon2_s);
                    
                    
                    ind=find(lon_k >= lon_2(i)-dlon & lon_k <= lon_2(i)+dlon ...
                        & lat_k >= lat_l2_i-dlat & lat_k < lat_u2(i)+dlat);
                elseif isnan(lat_u2(i))      % longitude at the lower end (lower latitude) is higher than that at the upper end (higher latitude)
                    cond_flag=2;
                    [var_n2, lon2_n]=min(abs(lon_2(i)-lon_n));
                    [var_s2, lon2_s]=min(abs(lon_2(i)-lon_s));
                    slp=(lat_n(lon2_n)-lat_s(lon2_s))/(lon_n(lon2_n)-lon_s(lon2_s));
                    lat_u2_i=slp*(lon_2(i)-lon_s(lon2_s))+lat_s(lon2_s);
                    ind=find(lon_k >= lon_2(i)-dlon & lon_k <= lon_2(i)+dlon ...
                        & lat_k >= lat_l2(i)-dlat & lat_k < lat_u2_i+dlat);
                else
                    cond_flag=3;
                    ind=find(lon_k >= lon_2(i)-dlon & lon_k <= lon_2(i)+dlon ...
                        & lat_k >= lat_l2(i)-dlat & lat_k < lat_u2(i)+dlat);
                end
                indx_i = [indx_i; ind];   % concatenated selected krig grid index array
                %         if lon_2(i) >= -133.15 & lon_2(i) <= -133.05
                %             disp(i)
                %         end
                %         if  isempty(ind)
                %             disp([i slp])
                %         end
            end
        end
        indx_i_uniq=unique(indx_i);
        indx = [indx; indx_i_uniq];
        
        data.in.US_CAN_mesh_reg(2).ind = indx_i_uniq;
        data.in.US_CAN_mesh_reg(2).lat = lat_k(indx_i_uniq);
        data.in.US_CAN_mesh_reg(2).lon = lon_k(indx_i_uniq);
        data.in.US_CAN_mesh_reg(2).tx_min = tx0;
        data.in.US_CAN_mesh_reg(2).tx_max = tx1;
    else
        lon_2=[];
        lat_l2=[];
        lon_2=[];
        lat_u2=[];
    end
        
    % ===========================================================================================
    %% region 3: paralell transects to latitudes west of QC IS
    %% get transect region definitions for region 3 of para.survey_year
    eval(['[tx0,tx1,tx_l,tx_r]=transect_region_def_' para.survey_year '(3);'])
    
    if ~isempty(tx0)
        % 2011 example
        % tx0=106;    % southern most transect number
        % tx1=128;  % northern most transect number
        % tx_l=[105.1 138.1 136.1 132.1 130.1 128.1];
        % tx_r=[105.1 144.4 138.4 136.4 134.4 124.6 124.9];
        
        % find upper and lower bounds of latitudes for the region
        ind0=find(tx == tx0);
        ind1=find(tx == tx1);
        lat_3a=min(min(lat_t(ind0)),min(lat_t(ind1)));       % minimum latitude for region 3
        lat_3b=max(max(lat_t(ind1)),max(lat_t(ind0)));
        %     lat_3a=min(lat_t(ind0));       % minimum latitude for region 3
        %     lat_3b=max(lat_t(ind1));
        
        %% specifies left (west) and right (east) region boundaries based on the transects
        %% #.1 = west end of transect
        %% #.4 = east end of transect
        %% #.6 = south end of transect
        %% #.9 = north end of transect
        
        % 2011 example
        % tx_l=[105.1 138.1 136.1 132.1 130.1 128.1];
        % tx_r=[105.1 144.4 138.4 136.4 134.4 124.6 124.9];
        %% west bound (left side) of each transect
        lat_w=[];
        lon_w=[];
        k = 1;
        for i=1:length(tx_l)
            ind0=find(floor(tx) == floor(tx_l(i)));
            if ~isempty(ind0)
                ind1=round(rem(tx_l(i),1)*10);
                switch ind1
                    case 1  % west end of transect
                        [val, ind]=min(lon_t(ind0));
                    case 4  % east end of transect
                        [val, ind]=max(lon_t(ind0));
                    case 6  % south end of transect
                        [val, ind]=min(lat_t(ind0));
                        val=lon_t(ind0(ind));
                    case 9  % north end of transect
                        [val, ind]=max(lat_t(ind0));
                        val=lon_t(ind0(ind));
                end
                lon_w(k)=val;
                lat_w(k)=lat_t(ind0(ind));
                if k> 1 & lat_w(k) == lat_w(k-1)
                    %               disp([i k floor(tx_l(i)) lat_w(k)])
                else
                    k=k+1;
                end
            end
        end
        %     lat_3=lat_3a:dlat:lat_3b;
        lat_3=max(lat_3a, min(lat_w)):dlat:min(lat_3b, max(lat_w)); % latitude array for region 1
        
        %% region #3 left side (west) bound in longitude
        lon_l3=interp1(lat_w,lon_w,lat_3);
        %     figure(1)
        %     hold on
        %     plot(lon_w,lat_w,'ok','markersize',12,'linewidth',2)
        
        
        %% east bound (right side) of each transect
        lat_e=[];
        lon_e=[];
        k = 1;
        for i=1:length(tx_r)
            ind0=find(floor(tx) == floor(tx_r(i)));
            if ~isempty(ind0)
                ind1=round(rem(tx_r(i),1)*10);
                switch ind1
                    case 1  % west end of transect
                        [val, ind]=min(lon_t(ind0));
                    case 4  % east end of transect
                        [val, ind]=max(lon_t(ind0));
                    case 6  % south end of transect
                        [val, ind]=min(lat_t(ind0));
                        val=lon_t(ind0(ind));
                    case 9  % north end of transect
                        [val, ind]=max(lat_t(ind0));
                        val=lon_t(ind0(ind));
                end
                lon_e(k)=val;
                lat_e(k)=lat_t(ind0(ind));
                if k> 1 & lat_e(k) == lat_e(k-1)
                    %               disp([i k floor(tx_r(i)) lat_e(k)])
                else
                    k=k+1;
                end
            end
            if i >= 4000
                figure(1)
                hold on
                plot(lon_e(i),lat_e(i),'sk','markersize',12)
                hold off
                disp(i)
            end
        end
        lat_3=max(lat_3a, min(lat_e)):dlat:min(lat_3b, max(lat_e)); % latitude array for region 1
        %% region #3 right side (west) bound in longitude
        lon_r3=interp1(lat_e,lon_e,lat_3);
        %     plot(lon_e,lat_e,'ok','markersize',12,'linewidth',2)
        %
        %     plot(lon_l3,lat_3,lon_r3,lat_3,'-k','linewidth',1.5)
        %     hold off
        indx_i = [];
        for i=1:length(lat_3)
            dlon(i)=dlat*cosd(lat_3(i));
            if i > length(lon_l3)
                lon_l3(i) = lon_l3(i-1);
            end
            if isnan(lon_l3(i))          % latitude at the west end is lower than that at the east end
                [var_max_w,max_ind_w]=max(lat_w);
                [var_max_e,max_ind_e]=max(lat_e);
                slp=(lon_w(max_ind_w)-lon_e(max_ind_e))/(var_max_w-var_max_e);
                lon_l3_i=slp*(lat_3(i)-var_max_e)+lon_e(max_ind_e);
                ind=find(lon_k >= lon_l3_i-dlon(i) & lon_k <= lon_r3(i)+dlon(i) ...
                    & lat_k >= lat_3(i)-dlat & lat_k < lat_3(i)+dlat);
            elseif isnan(lon_r3(i))      % latitude at the east end is lower than that at the west end
                [var_max_w,max_ind_w]=max(lat_w);
                [var_max_e,max_ind_e]=max(lat_e);
                slp=(lon_w(max_ind_w)-lon_e(max_ind_e))/(var_max_w-var_max_e);
                lon_r3_i=slp*(lat_3(i)-var_max_e)+lon_e(max_ind_e);
                ind=find(lon_k >= lon_l3(i)-dlon(i) & lon_k <= lon_r3_i+dlon(i) ...
                    & lat_k >= lat_3(i)-dlat & lat_k < lat_3(i)+dlat);
            else
                ind=find(lon_k >= lon_l3(i)-dlon(i) & lon_k <= lon_r3(i)+dlon(i) ...
                    & lat_k >= lat_3(i)-dlat & lat_k < lat_3(i)+dlat);
            end
            indx_i = [indx_i; ind];    % concatenated selected krig grid index array
            %     if  isempty(ind)
            %         disp([i slp])
            %     end
        end
        indx_i_uniq=unique(indx_i);
        indx = [indx; indx_i_uniq];
        
        data.in.US_CAN_mesh_reg(3).ind = indx_i_uniq;
        data.in.US_CAN_mesh_reg(3).lat = lat_k(indx_i_uniq);
        data.in.US_CAN_mesh_reg(3).lon = lon_k(indx_i_uniq);
        data.in.US_CAN_mesh_reg(3).tx_min = min(tx0, tx1);
        data.in.US_CAN_mesh_reg(3).tx_max = max(tx0, tx1);
    else
        lon_l3=[];
        lat_3=[];
        lon_r3=[];
        lat_3=[];
    end
end    % end of "if para.proc.source ~= 1"

acc_ind = [];
for i = 1:3
    if para.proc.source == i | para.proc.source == 3
        acc_ind = [acc_ind; data.in.US_CAN_mesh_reg(i).ind];
    end
end
data.in.US_CAN_mesh_acc_ind =acc_ind;

indx_uniq=unique(indx);

%figure(5); plot(sort(indx),'.');hold on;plot(indx_uniq,'or');hold off
if plot_flag == 1
    figure(2)
    %     clf
    plot(lon_k,lat_k,'.m','markersize',3)
    hold on
    plot(lon_l1,lat_1,'-g',lon_r1,lat_1,'-r', ...
        lon_2,lat_l2,'--g',lon_2,lat_u2,'--r', ...
        lon_l3,lat_3,'p-g',lon_r3,lat_3,'p-r', ...
        lon_k(indx_uniq),lat_k(indx_uniq),'.k','markersize',5,'linewidth',2)
    xlabel('Longitude')
    ylabel('Latitude')
    legend('Original Krig Grids','Regon 1 West Boundary','             East Boundary', ...
        'Regon 2 South Boundary','             North Boundary', ...
        'Regon 3 West Boundary','             East Boundary', ...
        'Non-extrapolated Krig Grids','location','northeast')
    hold off
    grid on
    figure(3)
    clf
    plot(lon_k,lat_k,'.m','markersize',3)
    hold on
    plot(lon_l1,lat_1,'-g',lon_r1,lat_1,'-r', ...
        lon_2,lat_l2,'--g',lon_2,lat_u2,'--r', ...
        lon_l3,lat_3,'p-g',lon_r3,lat_3,'p-r', ...
        lon_k(indx_uniq),lat_k(indx_uniq),'.k','markersize',5,'linewidth',2)
    xlabel('Longitude')
    ylabel('Latitude')
    grid on
    for i=1:nt
        ind=find(tx == uniq_tx(i));
        lat=lat_t(ind);
        lon=lon_t(ind);
        plot(lon,lat,'.-',lon1(i),lat1(i),'^m',lon2(i),lat2(i),'og',lon(1),lat(1),'sr','linewidth',2)
        if rem(uniq_tx(i),1) == 0
            text(mean(lon),mean(lat)+0.05,sprintf('%d',uniq_tx(i)),'color','r','fontweight','bold');
        else
            text(mean(lon),mean(lat)+0.05,sprintf('%3.1f',uniq_tx(i)),'color','r','fontweight','bold');
        end
        if legend_flag == 1 & i == nt
            legend('Original Krig Grids','Regon 1 West Boundary','             East Boundary', ...
                'Regon 2 South Boundary','             North Boundary', ...
                'Regon 3 West Boundary','             East Boundary', ...
                'Non-extrapolated Krig Grids', ...
                'Transect','Transect SW Bound','Transect NE Bound','Transect Start Pos','location','northeast');
        end
    end
    hold off
    figure(4)
    clr ='brg';
    marker = '.os';
    for i = 1:3
        if para.proc.source == i | para.proc.source == 3
            plot(data.in.US_CAN_mesh_reg(i).lon, data.in.US_CAN_mesh_reg(i).lat, [marker(i) clr(i)])
            if i == 1 hold on, end
        end
    end
    grid on
    hold off
end
drawnow
return