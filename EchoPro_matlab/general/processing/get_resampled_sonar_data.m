function [Tnew, accu_ind]=get_resampled_sonar_data(data,para,VL_BiomassInt0)
% transect reduction operation
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013
%% Modification:       03/10/2020   transect reduction will include trawl reduction


%% original transects
T0=VL_BiomassInt0(:,1);
T=sort(unique(T0),'ascend');
n0=length(T);
n=round((100-para.proc.transect_reduction_fraction)*n0/100);
if para.proc.start_transect == 2
    T = T(2:end);
end

if isfield(para.proc,'transect_reduction_mode') & para.proc.transect_reduction_mode == 1  % regular spaced transects
    if para.proc.transect_reduction_fraction ~= 50
        inc=round(n0/n);
        Tnew=T(1:inc:n0);
    else
        d0=median(VL_BiomassInt0(:,8));   % base transect spacing
        R1=6378.1;          % Equator radius
        R2=6356.8;          % Polar radius
        Rm=(R1^2*R2)^(1/3)/1.852;  % mean radius in nmi
        %%% double the transect spacing
        Tnew(1)=T(1);
        k=2;
        for i=3:length(T)-1   % 
%             d=mode(VL_BiomassInt0(find(T0 == T(i)),8));
            lat0_i=VL_BiomassInt0(find(T0 == Tnew(k-1)),5);   % previously selected new transect
            lat1_i=VL_BiomassInt0(find(T0 == T(i)),5);        % current transect
            lon0_i=VL_BiomassInt0(find(T0 == Tnew(k-1)),6);   % previously selected new transec
            lon1_i=VL_BiomassInt0(find(T0 == T(i)),6);        % current transect
            lat0=median(lat0_i);   % median latitude of previously selected new transect Tnew(k-1)
            lat1=median(lat1_i);     % median latitude of the current transect T(i)
            lon0=median(lon0_i);   % median longitude of previously selected new transect Tnew(k-1)
            lon1=median(lon1_i);     % median longitude of the current transect T(i)
            lat0_sd=sqrt(nanmean((lat0_i-lat0).^2));
            lat1_sd=sqrt(nanmean((lat1_i-lat1).^2));
            lon0_sd=sqrt(nanmean((lon0_i-lon0).^2));
            lon1_sd=sqrt(nanmean((lon1_i-lon1).^2));
            d_lat=Rm*(lat1-lat0)*pi/180;
            d_lon=Rm*cosd((lat0+lat1)/2)*(max(lon1_i)-min(lon1_i))*pi/180;
%             if i >= 50
%                 disp(i)
%             end
%             fprintf('i = %d\t T = %d\t  Lon_min = %6.4f\t Lon_max =%6.4f\n',i,T(i),min(lon1_i),max(lon1_i));
%% determine whether is parallel to latitude or parallel to longitude
            if abs((d0-d_lat)/d0) > 0.2  | abs((d0-d_lon)/d0) > 0.2
                if lat0_sd < lon0_sd & lat1_sd < lon1_sd
                    D=d_lat;
%                     fprintf('T1 = %d\t, T2 = %d\t, Lat0_sd = %6.4f\t  Lat1_sd = %6.4f\t Lon0_sd = %6.4f\t Lon1_sd = %6.4f\t D_lat = %6.3f (nm)\n', ...
%                         T(i-1),T(i),lat0_sd,lat1_sd,lon0_sd,lon1_sd,d_lat);
                elseif lon0_sd < lat0_sd & lon1_sd < lat1_sd
                    D=d_lon;
%                     fprintf('T1 = %d\t, T2 = %d\t, Lat0_sd = %6.4f\t  Lat1_sd = %6.4f\t Lon0_sd = %6.4f\t Lon1_sd = %6.4f\t D_lon = %6.3f (nm)\n', ...
%                         T(i-1),T(i),lat0_sd,lat1_sd,lon0_sd,lon1_sd,d_lon);
                else
                    D=d0;
%                     fprintf('T1 = %d\t, T2 = %d\t, Lat0_sd = %6.4f\t  Lat1_sd = %6.4f\t Lon0_sd = %6.4f\t Lon1_sd = %6.4f\t D_lat = %6.3f (nm)\t D_lon = %6.3f (nm)\n', ...
%                         T(i-1),T(i),lat0_sd,lat1_sd,lon0_sd,lon1_sd,d_lat,d_lon);
                end
            end
            if T(i) <= 130   % 2013 Transects
                if D > 1.5*d0  %| i == 100
                    Tnew(k)=T(i);
%                     fprintf('k = %d\t i = %d\t  Selected T = %d\n',k, i, Tnew(k));
                    k=k+1;
                end
            else
                Tnew(k)=T(i);
%                 fprintf('### k = %d\t i = %d\t  Selected T = %d\n',k, i, Tnew(k));
                k=k+1;
            end
        end
    end
else                                        % randomized transects
   ind=round(n0*rand(1,n));
   ind(ind > n0)=n0;
   ind(ind < 1)=1;
   ind(ind > length(T)) = length(T);
   Tindx=T(ind);
   Tnew=sort(unique(Tindx));
   while length(Tnew) < n
       ind=round(n0*rand(1,n-length(Tnew)));
       Tindx_i=setdiff(ind,Tnew);
       Tindx_ii = intersect(Tindx_i, T);
       Tnew=[Tnew; Tindx_ii(:)];
   end
end
n1=length(Tnew);
accu_ind=[];
for i=1:n1
    ind=find(T0 == Tnew(i));
    accu_ind=[accu_ind; ind];
end
fprintf('number of transects = %d out of %d\n',n1,n0);
VL_BiomassInt=VL_BiomassInt0(accu_ind,:);
return