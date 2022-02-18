% plot standard deviation and biomass 

figure(2)
indx = 1:length(var);
subplot(3, 1, 1)
plot(indx, var(indx).*Area0(indx)*1e-3, '.')
ylabel('Biomass (mt)')
subplot(3, 1, 2)
SD = sqrt(var1(indx)).*Area0(indx)*1e-3/sqrt(length(var1));
plot(indx, SD , '.')
ylabel('Standard Deviation')
subplot(3, 1, 3)
CV = sqrt(var1(indx)).*Area0(indx)./var(indx);
plot(indx, CV, '.')
xlabel('Sample Index')
ylabel('CV')

cv_bin =  0.001:0.001:0.5;

cnt = hist(CV, cv_bin);

figure(3)
indx_cv = 2:30;
subplot(2,1,1)
bar(cv_bin, cnt)
xlabel('CV')
ylabel('Counts')

subplot(2,1,2)
bar(cv_bin(indx_cv), cnt(indx_cv))
xlabel('CV')
ylabel('Counts')

figure(4)
ind = find(CV > 0.5);
plot(data.out.krig.lon, data.out.krig.lat, 'o', data.out.krig.lon(ind), data.out.krig.lat(ind), '.r')
grid on
xlabel('LONGITUDE','fontsize',16,'fontweight','bold')
ylabel('LATITUDE','fontsize',16,'fontweight','bold')

