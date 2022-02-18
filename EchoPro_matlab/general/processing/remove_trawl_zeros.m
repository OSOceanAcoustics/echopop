function out=remove_trawl_zeros(transect_info_filename,NASC_struct0)
%% INPUTS:
%% EV_export_filepath = directory for the exported  (interval, layer, cell, analysis) files from EV files
%% NASC_struct = NASC data structure including all interval points
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       07/14/2015

global para

%% data from EchoView Export
%% check if there is a VL discontinuity due to different VL between Canadian and US Ships
ind=find(diff(NASC_struct0(:,3)) < -4000);
nL=length(ind);
if ~isempty(ind)
    ind1=find(diff(NASC_struct0(:,3)) > 4000);
    if nL == length(ind1)+1
        for i=1:nL
            if i == nL
               VL_add(i)=NASC_struct0(ind(i),4)+500;
               NASC_struct0(ind(i)+1:end,3:4)=NASC_struct0(ind(i)+1:end,3:4)+VL_add(i);
            else
               VL_add(i)=NASC_struct0(ind(i),4)+500;
               NASC_struct0(ind(i)+1:ind1(i),3:4)=NASC_struct0(ind(i)+1:ind1(i),3:4)+VL_add(i);
            end
        end
    else
        fprintf('problems with VL and transects numbering!! \n No action is taken !!\n ')
        out=NASC_struct0;
    end
end

[vl, sort_ind_nasc]=sort(NASC_struct0(:,3));
NASC_struct=NASC_struct0(sort_ind_nasc,:);
ind=find(diff(NASC_struct(:,3)) == 0);
NASC_struct(ind,:)=[];

VL0=NASC_struct(:,3);
VL1=NASC_struct(:,4);
out=NASC_struct;

figure(4)
plot(NASC_struct(:,6),NASC_struct(:,5),'.')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
axis([-136 -120 35 56])
grid
hold on

%% Ecent log VL: ST, BT, RT, and ET
d00=xlsread(transect_info_filename);
%% check if there is a VL discontinuity due to different VL between Canadian and US Ships

n=size(d00,1);
for i=1:n
    tmp=num2str(d00(i,2));
    sel_year_ind_flag(i)=strcmp(tmp(1:4),para.survey_year);
end
sel_year_ind= find(sel_year_ind_flag ==1);

d0=d00(sel_year_ind,:);

ind=find(diff(d0(:,7)) < -4000);
nL1=length(ind);
if nL1 ~= nL
    fprintf('VLs from Echoview export are not consistent with those from event log !!\n');
end
if ~isempty(ind)
    ind1=find(diff(d0(:,7)) > 4000);
    if nL == length(ind1)+1
        for i=1:nL
            if i == nL
               d0(ind(i)+1:end,6:7)=d0(ind(i)+1:end,6:7)+VL_add(i);
            else
               d0(ind(i)+1:ind1(i),6:7)=d0(ind(i)+1:ind1(i),6:7)+VL_add(i);
            end
        end
    else
        fprintf('problems with VL and transects numbering!! \n No action is taken !!\n ')
    end
end

% if ~isempty(ind)
%     d0(ind(1)+1:end,6:7)=d0(ind(1)+1:end,6:7)+VL_add;
% end

[vl, sort_ind_event_log]=sort(d0(:,6));
d=d0(sort_ind_event_log,:);

T=d(:,3);     % transect number
VLs=d(:,6);   % Start event
VLe=d(:,7);   % End event

Tuniq=unique(T);     % unique transects

nT=length(Tuniq);  % total number of transects
fprintf('\n -------------------------------------------------\n')
ind_zeros=[];
for i=1:nT
    if i == 1
        ind=find(VL1 < VLs(1));
    elseif i == nT
        ind=find( VL0 > VLe(end));
    else
        ni=find(T == Tuniq(i));
        nTi=length(ni);
        ind=[];
%         if nTi > 1
%             disp([i nTi]) 
%         end
        for j=1:nTi
            if j == 1
                indj=find( VL0 > VLe(ni(1)-1) & VL1 < VLs(ni(1)));
%                 disp(num2str([ni(1) VLe(ni(1)-1) VLs(ni(1))]))
%             elseif j == nTi
%                 indj=find( VL0 > VLe(ni(nTi)) & VL1 < VLs(ni(nTi)+1));
%                 disp(num2str([ni(nTi) VLe(ni(nTi)) VLs(ni(nTi)+1)]))
            else
                indj=find( VL0 > VLe(ni(j-1)) & VL1 < VLs(ni(j)));
%                 disp(num2str([ni(j-1) ni(j) VLe(ni(j-1)) VLs(ni(j))]))
            end
            ind=[ind;indj];
        end
    end
    ind_zeros=[ind_zeros ; ind];
    fprintf('T = %d\t, nn = %d\t  tot nn = %d\n',Tuniq(i),length(ind),length(ind_zeros))
    plot(out(ind,6),out(ind,5),'or')

end
        
out(ind_zeros,:)=[];

plot(out(:,6),out(:,5),'.r')

hold off
return