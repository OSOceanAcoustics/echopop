function [US_CAN_dat,US_CAN_mesh_cor]=convert_data(dat1,Lsmooth,US_BCw_border,BCw_BCe_border,Ref_lon,US_CAN_mesh)
%% lat and lon coordinates tranformation
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013


lat1=dat1(:,1);
lon1=dat1(:,2);
indx=find(lon1 > 0);    % Longitude west is negative by convention
lon1(indx)=-lon1(indx);

biom1=dat1(:,3);

%%%%%%%%%%%%%%%%%%% corrections  %%%%%%%%%%%%%%%%%%%%%%%%%
US_lat_cor=Lsmooth(1:US_BCw_border,1)+eps;
US_lon_cor=Lsmooth(1:US_BCw_border,2)-Ref_lon;
CAN_BCw_lat_cor=Lsmooth(US_BCw_border-1:BCw_BCe_border-1,1);
CAN_BCw_lon_cor=Lsmooth(US_BCw_border-1:BCw_BCe_border-1,2)-Ref_lon;
CAN_BCe_lat_cor=Lsmooth(BCw_BCe_border-1:end,1);
CAN_BCe_lon_cor=Lsmooth(BCw_BCe_border-1:end,2)-Ref_lon;

%% interval-based biomass location correction in transects
% US & CAN
US_CAN_lat_cor=[US_lat_cor; CAN_BCw_lat_cor(3:end)];
US_CAN_lon_cor=[US_lon_cor; CAN_BCw_lon_cor(3:end)];
%% longitude correction values based in interpolation
US_CAN_bio_lon_del=interp1(US_CAN_lat_cor,US_CAN_lon_cor,lat1);
%% corrected longitudes for both US & CAN transects
US_CAN_bio_lon_cor=lon1-US_CAN_bio_lon_del;

US_CAN_dat=[lat1 US_CAN_bio_lon_cor biom1];
indx=find(isnan(US_CAN_bio_lon_cor) == 1);
US_CAN_dat(indx,:)=[];



%% mesh location correstion for biomass estimates
% US & CAN: longitude correction values for mesh grids based on interpolation
US_CAN_mesh_lon_del=interp1(US_CAN_lat_cor,US_CAN_lon_cor,US_CAN_mesh(:,1));
%% corrected longitudes for both US & CAN meshes (grid cells)
US_CAN_mesh_lon_cor=US_CAN_mesh(:,2)-US_CAN_mesh_lon_del;
%% [lat lon]
US_CAN_mesh_cor=[US_CAN_mesh(:,1)  US_CAN_mesh_lon_cor];


return