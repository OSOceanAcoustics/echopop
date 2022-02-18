function [xb,yb]=find_west_bound(data)
% function find_west_bound is to find the west boundary of the transects
% as a function of latitude (normalized distance in y-direction)

% nb=200;                 % number of samples along y-axis
% yb=linspace(min(data.in.y),max(data.in.y),nb);

T=data.final.table.biomass(:,1);
Lat=data.final.table.biomass(:,5);
Lon=data.final.table.biomass(:,6);
Tuniq=unique(T);

for i=1:length(Tuniq)
    ind=find(T == Tuniq(i) & Lat < 51);
    ind(ind > length(data.in.x)) = length(data.in.x);
    if ~isempty(ind)
        [xb(i), indi]=nanmin(data.in.x(ind));
        yb(i)=data.in.y(ind(indi));
    end
end
ind=find(xb == 0);
xb(ind)=[];
yb(ind)=[];