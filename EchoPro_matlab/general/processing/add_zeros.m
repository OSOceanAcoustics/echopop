function   add_zeros(transect_num)
%% added zeros to the offshore end of the US transects
%% Dezhang Chu
%% Created on 12-20-2013
%% Modified on 12-22-2013

global data para

% data.final.table.biomass=data.final.table.biomass0;
indx_end=1:para.proc.end_index;
na=para.proc.add_zero_num;                       % number of added zeros on the west extend

X=data.final.table.biomass(indx_end,1);
VLst=data.final.table.biomass(indx_end,3);
VLed=data.final.table.biomass(indx_end,4);
lat=data.final.table.biomass(indx_end,5);
lon=data.final.table.biomass(indx_end,6);
NASC=data.final.table.biomass(indx_end,9);

data.final.table.biomass0=data.final.table.biomass;
Xuniq=unique(X);
VLint=median(diff(VLst));

for i=1:length(Xuniq)
    indx=find(X == Xuniq(i));
    var=2*length(indx)*VLint;
    ind=find(VLst(indx) > median(VLst(indx))+ var | VLst(indx) < median(VLst(indx))- var);  % remove bad VL data
    if ~isempty(ind)
%         figure(1)
%         plot(VLst(indx),'.')
        for j=1:length(ind)
            indx_sm=max(1,ind(j)-2):min(ind(j)+2,length(indx));
            VLst(indx(ind(j)))= median(VLst(indx(indx_sm)));
        end
%         hold on
%         plot(VLst(indx),'or')
%         title(sprintf('T = %d',Xuniq(i)))
%         grid
%         hold off
    end
    ind=find(VLed(indx) > median(VLed(indx))+ var | VLed(indx) < median(VLed(indx))- var);  % remove bad VL data
    if ~isempty(ind)
%         figure(1)
%         plot(VLed(indx),'.')
        for j=1:length(ind)
            indx_sm=max(1,ind(j)-2):min(ind(j)+2,length(indx));
            VLed(indx(ind(j)))= median(VLed(indx(indx_sm)));
        end
%         hold on
%         plot(VLed(indx),'or')
%         title(sprintf('T = %d',Xuniq(i)))
%         grid
%         hold off
    end
    if median(diff(lon(indx))) < 0
        VLinc(i)=-VLint;                        % VL increment
        VL_st_add(i)=max(VLst(indx));           % VL start sample index on transect Xuniq(i)
        VL_ed_add(i)=max(VLed(indx));           % VL end sample index on transect Xuniq(i)
    else
        VLinc(i)=VLint;                         % VL increment
        VL_st_add(i)=min(VLst(indx));
        VL_ed_add(i)=min(VLed(indx));
    end
    [lon_min,ind_lon_min]=sort(lon(indx));
    lat_m(i)=median(lat(indx));                 % median latitude on transect Xuniq(i)
    lon_m(i)=min(lon_min);                      % longitude of the westmost samples on transect Xuniq(i)
    lon_dif(i)=abs(median(diff(lon(indx))));    % median longitude difference between samples on transect Xuniq(i)
    lat_add(i,1:na)=lat_m(i)*ones(1,na);        % added latitudes for added samples on transect Xuniq(i)
    if VLinc(i) > 0   % West --> East
        lon_add(i,1:na)=lon_m(i)-lon_dif(i)*[1:na];  % added longitudes for added samples on transect Xuniq(i)
        ins(i)=min(indx);                            % transect start sample index on original table
    else               % East --> West
        lon_add(i,1:na)=lon_m(i)-lon_dif(i)*[na:-1:1];
        ins(i)=max(indx);                            % transect end sample index on original table
    end
end

%% construct the modified table with added zeros on the west extend of the transects
if nargin == 0
   biomass=ones(size(data.final.table.biomass,1)+length(Xuniq)*na,size(data.final.table.biomass,2));
else
   biomass=ones(size(data.final.table.biomass,1)+length(transect_num)*na,size(data.final.table.biomass,2));
end


for i=1:length(transect_num)
    if i > 1
        indx_Ti=find(X < Xuniq(transect_num(i)) & X > Xuniq(transect_num(i-1)));
    else
        indx_Ti=find(X < Xuniq(transect_num(i)));
    end
    if ~isempty(indx_Ti)
        ind_st=(i-1)*na;
        indx_new=ind_st+indx_Ti;
        biomass(indx_new,:)=data.final.table.biomass0(indx_Ti,:);
    end
    
    indx0=find(X == Xuniq(transect_num(i)));
    % construc matrix for added zeros 
    dat_add=repmat(data.final.table.biomass0(ins(i),:),na,1);
    dat_add(:,1)=Xuniq(transect_num(i));
    dat_add(:,5)=lat_add(transect_num(i),:)';                         % updated latitude 
    dat_add(:,6)=lon_add(transect_num(i),:)';                         % updated longitude 
    if VLinc(transect_num(i)) > 0  % West --> East
        if  transect_num(i) ==1
            indx_orig=1:ins(i);                         % indices for original table
            indx_new=indx_orig;                         % indices for modified table
            indx_add=length(indx0)+1:length(indx0)+na;  % indices for added zeros
        else
            ind_st=(i-1)*na+ins(transect_num(i));
            indx_orig=ins(transect_num(i)):ins(transect_num(i))+length(indx0)-1;
            indx_new=ind_st+na:ind_st+na+length(indx0)-1;
            indx_add=ind_st+[0:na-1];
        end
        dat_add(:,3)=VL_st_add(transect_num(i))-[na:-1:1]'*VLint;
        dat_add(:,4)=VL_ed_add(transect_num(i))-[na:-1:1]'*VLint;
        index_all(i,1:4)=[indx_add([1 end])  indx_new([1 end])];
    else   % Esat --> West
        if  transect_num(i) == 1
            indx_orig=1:ins(i);
            indx_new=indx_orig;
        else
            ind_st=(i-1)*na+ins(transect_num(i));
            indx_orig=ins(transect_num(i))-length(indx0)+1:ins(transect_num(i));
            indx_new=ind_st-length(indx0)+1:ind_st;
        end
        dat_add(:,3)=VL_st_add(transect_num(i))+[1:na]'*VLint;
        dat_add(:,4)=VL_ed_add(transect_num(i))+[1:na]'*VLint;
        indx_add=indx_new(end)+[1:na];
        index_all(i,1:4)=[indx_new([1 end]) indx_add([1 end]) ];
    end
    biomass(indx_new,:)=data.final.table.biomass0(indx_orig,:);
    biomass(indx_add,:)=dat_add;
end
[n,m]=size(biomass);
n1=max(index_all(end,1:4));
biomass(n1+1:n,1:m)=data.final.table.biomass0(max(indx_orig)+1:end,:);

data.final.table.biomass=biomass;
return