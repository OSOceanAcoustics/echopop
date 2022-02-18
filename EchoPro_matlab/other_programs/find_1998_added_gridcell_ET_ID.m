% find 1998 added grids & associate ET_ID

clear

filename1='V:\acousticsData_nwcfs2\Historical Summary (for Kriging)\Kriging files & parameters\Kriging grid files\krig_grid2_5nm_cut_centroids_1998.xlsx';
filename2='V:\acousticsData_nwcfs2\Historical Summary (for Kriging)\Outputs\Historical Outputs (KS Stratification with aged data)\with extrapolation\1998\EchoPro_kriged_output-26-Jan-2016_0.xlsx';
d1=xlsread(filename1);
d2=xlsread(filename2);
lat1=d1(:,3);
lon1=d1(:,4);
lat2=d2(:,1);
lon2=d2(:,2);

% ind1998=1:length(ET_ID1998);
% ind2013=1:length(ET_ID2013);
% 
% [C,IA]=setdiff(ET_ID2013,ET_ID1998);


figure(1)
plot(lon1,lat1,'.',lon2,lat2,'or')

figure(2)
plot(lat1-lat2,'.-')

figure(3)
plot(lon1-lon2,'.-')