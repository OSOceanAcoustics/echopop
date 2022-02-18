function get_fish_school_characteristics
%% get parameters associated with 2D fish schools 
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013
%% Last Modification:       3/13/2020
global data para

%% detect regions where there are fish or NASC > 0
NASC=data.final.table.biomass(:,9);
p=sign(NASC);
ind1=find(diff(p)>0);
ind2=find(diff(p)<0);
if length(ind1) < length(ind2) 
   if ind1(1) > ind2(1)
       ind1=[1; ind1];
   else
       disp('either start or end region time coincides with either start or end transect time!!')
   end
elseif length(ind1) > length(ind2)
   if ind2(end) < ind1(end)
       ind2=[ind2; ind1(end)+1];
   else
       disp('either start or end region time coincides with either start or end transect time!!')
   end    
else
   if ind1(1) > ind2(1)   % NASCs are > 0 for the first and the last cells
       ind1 = [1 ind1];
       ind2 = [ind2 length(p)];
   end
end
n=length(ind1);
totW=nansum(data.final.table.biomass(:,15));
for i=1:n
    ind=find(data.final.table.biomass(:,1) == data.final.table.biomass(ind2(i),1) & data.final.table.biomass(:,2) == data.final.table.biomass(ind2(i),2));
    ind=ind1(i)+1:ind2(i);
    data.final.table.school(i,1)=data.final.table.biomass(ind2(i),1);        % transect number
    data.final.table.school(i,2)=data.final.table.biomass(ind2(i),2);        % Region ID
    data.final.table.school(i,3)=nanmean(data.final.table.biomass(ind,5));   % lat
    data.final.table.school(i,4)=nanmean(data.final.table.biomass(ind,6));   % lon
    data.final.table.school(i,6)=data.final.table.biomass(ind2(i),4)-data.final.table.biomass(ind1(i)+1,3);   % fish school length (nm)
    data.final.table.school(i,9)=nanmean(data.final.table.biomass(ind,9));   % fish school mean NASC
    tot_fish_no(i)=nansum(data.final.table.biomass(ind,12));                 % total number of fish in the fish school
    data.final.table.school(i,5)=nanmean(data.final.table.biomass(ind,22));  % fish school depth (m)
    data.final.table.school(i,7)=nanmean(data.final.table.biomass(ind,23));  % fish school height (m)
    transect_spacing(i)=nanmean(data.final.table.biomass(ind,24));           % trasect spacing
    data.final.table.school(i,11)=nansum(data.final.table.biomass(ind,15));  % total biomass of fish in the fish school (kg)
    data.final.table.school(i,12)=nanmean(data.final.table.biomass(ind,8));  % fish school bottom depth (m)
    sub_tot=nansum(data.final.table.school(1:i,11));
    sub_W=nansum(data.final.table.biomass(1:ind(end),15));
    if abs(sub_tot-sub_W) > 1e-5
        fprintf('i = %d\t Wschool = %f\t  W = %f\t  ind1 = %d\t  ind2 = %d\n',i,sub_tot,sub_W,ind(1),ind(end));
        fprintf('\n')
    end
end
% depth is deeper than the maximum recording depth
ind=find(data.final.table.school(:,12)<0);
data.final.table.school(ind,12)=max(data.final.table.school(:,12));
D=data.final.table.school(:,5);                                                     % fish layer depth
H=data.final.table.school(:,7);                                                     % fish layer height (thickness)
W=D*para.acoust.bw(para.acoust.freq_ind);                                           % mean width at depth D
data.final.table.school(:,8)=1852*data.final.table.school(:,6).*H;                  % fish school 2D area (m^2)
Vinso=W.*data.final.table.school(:,8);                                              % total insonified volume fraction
Vschool=1852*transect_spacing(:).*data.final.table.school(:,8);                     % total estimated school volume within transect spacing
data.final.table.school(:,10)=tot_fish_no(:)./Vschool;                              % fish school numerical density (no. fish/m^3)

data.final.table.school_description={'Transect','Region no.','Lat','Lon','Layer Depth','Layer Length','Height','Area','NASC','Fish Density','Weight','Bottom Depth'};
