%% distance

clear
a=6378.1;          % Equator radius in km
b=6356.8;          % Polar radius in km

n=100000;
th=linspace(0,90,n);
th0=30:1/6:60;
km2nm=1/1.852;
R0=(a*a*b)^(1/3);
dth=th(2)-th(1);
dth0=th0(2)-th0(1);
r=a*b./sqrt((a*sind(th)).^2+(b*cosd(th)).^2);
dr=r*dth*pi/180;
for i=1:length(th0)
    ind=find(th > th0(i)-dth0/2 & th < th0(i)+dth0/2);
    D_lat(i)=sum(dr(ind))*km2nm;
end
D0_lat=R0*dth0*pi/180*km2nm;
figure(1)
plot(th0,D_lat,'.-',th0,D0_lat,'or')
xlabel('Lat (deg)')
ylabel('Transect Spacing (km)')
% axis([min(th0) max(th0) 9 11])
    
fprintf(' Mean D = %6.4f (nm),\t s.d = %6.4f (nm) \n',mean(D_lat),std(D_lat))
fprintf('Max diff = %6.4f (nm),  Rel Diff = %4.2f%%\n',max(abs(D_lat-D0_lat)),100*(D0_lat/mean(D_lat)-1))

