function  CV = jolly_hampton_CV(para, data)
% provide CV based on Jolly & Hampton Variance or CV analysis
% Jolly, G. M., and I. Hampton. 1990. "A stratified random transect design
% for acoustic surveys of fish stocks." Can. J. Fish. Aquat. Sci. 47: 1282-1 291.
%
% 

year = para.survey_year;
extrapolation = para.proc.extrapolation;  % 0 = no-extrapolation; 1 = extrapolation
% stratification=2;      % 1 = Kolmogorov-Smirnov goodness-of-fit startification
%                        % 2 = INPFC

nmi2m=1852;            % nautical mile to meters

fac=para.proc.JH_fac;  % fraction of the number of transects randomly selected in each stratum
if fac == 1
    nr=1;              % number of realization
else
    nr=10000;          % number of realization
end

% data_root = 'N:\Survey.Acoustics\acousticsData_nwcfs2\Historical Summary (for Kriging)\Outputs\';

%% add an additional stratum
% lat_INPFC=[36  38.0000 40.5  43.0000   45.7667   48.5000  51   55.0000]; % more strata
lat_INPFC = [36         40.5      43.000      45.7667     48.5       55.0000];  % INPFC
% lat_INPFC = [36          40.5            45.7667            51   55.0000]; % fewer strata
% lat_INPFC = [36          40.5           ];
ns=length(lat_INPFC);

fprintf('\n ------- Jolly-Hampton Variance (CV) Analysis .... \n');
if para.proc.kriging ~= 1
    fprintf('Analysis on Un-kriged Biomass\n')
    % exclude/include zag transects (T > 1000)
    d0 = data.final.table.biomass;
    ind0=find(d0(:,1) < 2000);  % exclude zag transects
    d=d0(ind0,:);
    T=d(:,1);               % transect number
    lat=d(:,5);             % latitude
    lon=d(:,6);             % longitude
    if size(d,2) > 21
        spacing=d(:,24);        % transect spacing
    else
        spacing = 10*ones(size(d,1));
    end
    %% unique transect
    uniq_T=unique(T);
    nT=length(uniq_T);      % number of unique transect spacing
    % longitude increment alone the transect
    lon_inc=0.05;
    for i=1:nT
        ind=find(T == uniq_T(i));       
        lon2 = lon(ind);                % indeices of longitudes on the ith transect
        biomass_j = d(ind, 15);         % biomass of samples on the ith transect
        
        biomass(i)=sum(biomass_j);      % total biomass on the ith transect
        %% distance of the processed transect
        dist(i)=m_idist(min(lon2),mean(lat(ind)),max(lon2),mean(lat(ind)))/nmi2m;
        area(i)=dist(i)*mean(spacing(ind));   % area of the ith transect covers (+- transect spaning x transect length)
        latitude(i)=mean(lat(ind));           % mean latitude of the ith transect 
        longitude_min(i)=min(lon2);           % minimum longitide of the ith transect
        longitude_max(i)=max(lon2);           % maximum longitide of the ith transect
    end
    if min(latitude) >= lat_INPFC(1)   % some years, the southernmost transect is north of the INPFC(1)
        lat_INPFC(1)=[];
        ns = ns-1;
    end    
    tot_area=0;
    for i=1:ns
        % transect indices within the ith stratum
        if i == 1
          strata(i).ind=find(latitude < lat_INPFC(i));
        else
          strata(i).ind=find(latitude >= lat_INPFC(i-1) & latitude < lat_INPFC(i));
        end
        if ~isempty(strata(i).ind)
            strata(i).n=length(strata(i).ind);  % number of transects in ith stratum
            strata(i).area=sum(area(strata(i).ind));
            lon_INPFC(i,1:2)=[min(longitude_min(strata(i).ind))-2 max(longitude_max(strata(i).ind))+2];
            tot_area=tot_area+strata(i).area;
        else
            strata(i).n=0;
            strata(i).area=0;
            lon_INPFC(i,1:2)=nan;
        end
    end
else
    %% vitual "transects" based on the latitude
    fprintf('Analysis on Kriged Biomass\n')
    %% un-kriged
%     if extrapolation == 0   % no extrapolation
%         fprintf('no extrapolation\n')
%         if stratification == 1  % no extrapolation & KS - stratification
%             file_dir = [data_root 'Historical Outputs (KS Stratification with aged data)\without extrapolation\' year '\'];
%             fprintf('KS \n')
%         else   % no-extrapolation & INPFC
%             file_dir = [data_root 'Historical Outputs (INPFC with aged data)\without extrapolation\' year '\'];
%             fprintf('INPFC \n')
%         end
%     else   % extrapolation
%         fprintf('with extrapolation\n')
%         if stratification == 1  % no extrapolation & KS - stratification
%             file_dir = [data_root 'Historical Outputs (KS Stratification with aged data)\with extrapolation\' year '\'];
%             fprintf('KS \n')
%         else
%             file_dir = [data_root 'Historical Outputs (INPFC with aged data)\with extrapolation\' year '\'];
%             fprintf('INPFC \n')
%         end
%     end
%     files = dir([file_dir 'EchoPro_kriged_output*0.xlsx']);
%     filename = [file_dir files.name];
%     d=xlsread(filename);
    
%     toc
    d = data.final.table.kriged_biomass0;
    [m,n]=size(d);
    lat=d(:,1);
    lon=d(:,2);
    %% latitude array with equal increment
    n_transect_per_lat = 5;      % number of "virtual transects" within a latitude degree 
    lat_eq_inc=round(d(:,1)*n_transect_per_lat+0.5)/n_transect_per_lat;   % 
    lon_inc=0.05;
    uniq_lat_eq_inc=unique(lat_eq_inc);  % unique equal-spacing transects 
    dis=zeros(length(lat),1);
    d1=zeros(m,n+2);
    nlat=length(uniq_lat_eq_inc);
    
    lat_inc=round(d(:,1)*5+0.5)/5;
    lon_inc=0.02;
    uniq_lat_inc=unique(lat_inc);
    dis=zeros(length(lat),1);
    d1=zeros(m,n+2);
    nlat=length(uniq_lat_inc);
    m2nmi = nmi2m;
    index_ij=[];
    %% loop through latitudes or "transects"
    for i=1:nlat
        ind=find(lat_eq_inc == uniq_lat_eq_inc(i));
        lon2 = lon(ind);                % indeices of longitudes on the ith transect
        biomass(i)=nansum(d(ind,10));
        % transect length
        dist(i)=m_idist(min(lon2),uniq_lat_eq_inc(i),max(lon2),uniq_lat_eq_inc(i))/nmi2m;
        if i == 1 | i == nlat
            area(i)=dist(i)*mean(diff(uniq_lat_eq_inc))*60;
        else
            area(i)=dist(i)*(uniq_lat_eq_inc(i)-uniq_lat_eq_inc(i-1))*60/2;
        end
        latitude(i)=uniq_lat_eq_inc(i);
        longitude_min(i)=min(lon2);
        longitude_max(i)=max(lon2);      
    end
    if min(latitude) >= lat_INPFC(1)   % some years, the southernmost transect is north of the INPFC(1)
        lat_INPFC(1)=[];
        ns = ns-1;
    end
    tot_area=0;
    for i=1:ns
        if i == 1
          strata(i).ind=find(latitude < lat_INPFC(i));
        elseif i ~= ns
          strata(i).ind=find(latitude >= lat_INPFC(i-1) & latitude < lat_INPFC(i));
        else
          strata(i).ind=find(latitude >= lat_INPFC(i-1) & latitude <= lat_INPFC(i));
        end
        strata(i).n=length(strata(i).ind);  % number of transects in the ith stratum
        strata(i).area=sum(area(strata(i).ind));  % area of the ith stratum
        lon_INPFC(i,1:2)=[min(longitude_min(strata(i).ind))-2 max(longitude_max(strata(i).ind))+2];
        tot_area=tot_area+strata(i).area;   % accumulative area from 1st to the ith stratum
    end
%     lat_ind=1;
%     lon_ind=2;
end


%% Jolly-Hampton algorithm
for ii=1:nr   % loop through realizations
    index_i=[];
    for i=1:ns % loop through strata
        %% number of transects within the current stratum
        ni(i)=round(fac*strata(i).n);
        %% ni(i) transect indeices selected randomly within the ith stratum
        indx_i=sort(min(max(1,round(strata(i).n*rand(1,ni(i)))),strata(i).n));
        % if same transects are selelected, randomly select different one(s)
        % until the number of the selected transects = ni(i)
        ind=find(diff(indx_i) == 0);
        indx_i(ind)=[];
        while ~isempty(ind)
            for j=1:length(ind)
                indx_ij=sort(min(max(1,round(strata(i).n*rand(1,1))),strata(i).n));
                indx_i=sort([indx_i indx_ij]);
            end
            ind=find(diff(indx_i) == 0);
            indx_i(ind)=[];
        end
        %     disp([indx_i(1) indx_i(end) length(ind) length(indx_i) strata(i).n])
        index_i=[index_i strata(i).ind(indx_i)];
        Lij=dist(strata(i).ind(indx_i));        % length of the jth transect in ith stratum
        latitude_ij=latitude(strata(i).ind(indx_i));
        biomass_ij=biomass(strata(i).ind(indx_i));  % biomass of the selected transects in the ith stratum
        %     hdl=plot(-125*ones(size(latitude_ij)),latitude_ij,'og');
        %     pause(1)
        %     delete(hdl)
        strata(i).biomass=biomass_ij;           % biomass of the selected transects in ith stratum
        strata(i).L=Lij;                        % transect length of the jth transect in ith stratum
        strata(i).latitude=mean(latitude_ij);
        strata(i).wgt=Lij/mean(Lij);             % transect length weithing factor of the jth transect in ith stratum (2-7-2018)
%         strata(i).wgt=Lij/sum(Lij);             % transect length weithing factor of the jth transect in ith stratum
        strata(i).rhom_ij=biomass_ij./Lij;      % normalized biomass of the jth transect in ith stratum
        strata(i).rhom_i=sum(biomass_ij.*Lij)/sum(Lij);    % transect-length-normalized mean density in the ith stratum (2-7-2018)
 %       strata(i).rhom_i=sum(biomass_ij)/(ni(i)*sum(Lij));    % transect-length-normalized mean density in the ith stratum
        area_i(i)=strata(i).area;               % total area in the ith stratum
        rhom_i(i)=strata(i).rhom_i;             
        %% variance of the transect-length weighted biomass within the stratum
        if ni(i) ~= 1
            var_rhom_i(i)=nansum(strata(i).wgt.^2.*(strata(i).rhom_ij-strata(i).rhom_i).^2)/(ni(i)*(ni(i)-1));
        else
            var_rhom_i(i)=nansum(strata(i).wgt.^2.*(strata(i).rhom_ij-strata(i).rhom_i).^2)/(ni(i)*(ni(i)));
        end
        biomass_i(i)=nansum(biomass_ij);
        nt(i)=strata(i).n;                      % number of transects in the ith stratum
        X_length(i)=sum(Lij);                   % total length of the selected ni(i) transects in the ith stratum
    end
    var_rho=sum(area_i.^2.*var_rhom_i)/sum(area).^2;
    biomass_m(ii)=sum(biomass_i.*area_i)/sum(area_i);
    %% area weighted variance of the "transect-length weighted biomass"
    CV(ii)=sqrt(nansum(var_rhom_i.*area_i.^2))/nansum(area_i.*rhom_i);
    rhom_all=nansum(area_i.*rhom_i)/sum(area);
    biomass_m_ave(ii)=rhom_all*sum(area_i);
end

return
