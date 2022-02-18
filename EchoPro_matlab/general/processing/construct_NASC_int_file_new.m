function out=construct_NASC_int_file(species_name,haul_strata_filename,EV_export_filepath,transect_reg_haul_filename,max_spacing)
%% INPUTS:
%% specise_name = name of the species to be processed
%% haul_strata_filename = filename of the table file: [ haul #,  strata #, weight]
%% EV_export_filepath = directory for the exported  (interval, layer, cell, analysis) files from EV files
%% transect_reg_haul_filename = filename of the table fiel: [transect #, region #, assigned haul #]
%% max_spacing = maximum spacing (kt) between the two transects
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Modified:       11/29/2013
%% Midified:       1/6/2021 to handle: Newer Echoview exports "class" without double quotes 

global para

%% for 2017 - IVC overlap transects
US_Ship_ID = para.US_Ship_ID;
ex_US_T = [];
% ex_CAN_T = [95:99 1096 1098];
ex_CAN_T = [];

NASC_threshold=20000;
if para.proc.stratification_index ~= 10
   haul_strata=xlsread(haul_strata_filename,para.proc.stratification_index+1);
else
   haul_strata=xlsread(haul_strata_filename,'length strata byhaul_1stratum');
end
T_reg_haul=xlsread(transect_reg_haul_filename);

files0=dir([EV_export_filepath '\*(intervals).csv']);
file_cells0=dir([EV_export_filepath '\*(cells).csv']);
file_layers0=dir([EV_export_filepath '\*(layers).csv']);

if para.proc.ordered_by_transects == 1
    % order files according to the transect numbers for 2015 when US & CAN
    % interleave the transects
    file_indx=sort_files_by_transect_num(files0);
    files=files0(file_indx);
    file_cells=file_cells0(file_indx);
    file_layers=file_layers0(file_indx);
end

nf=length(files);
tot_ind=0;
for i=1:nf
    ind_st=regexp(files(i).name,'-T')+2;
    ind_ed=regexp(files(i).name,'-Z')-1;
    TX(i)=str2num(files(i).name(ind_st:ind_ed));
end
%%% sort Transects ascendingly
[TX_sort, sort_indx]=sort(TX,'ascend');
high_NASC_ind = 1;
for i=1:nf    % loop through files or transects
    fprintf('T = %d\n',TX_sort(i))
%     if TX_sort(i) >= 95 & TX_sort(i) < 100
%         disp(i)
%     end
    T_reg(i).reg=[];
    %% intervel file
    dat=xlsread([EV_export_filepath '\' files(sort_indx(i)).name]);
    %% cell file
    [dat_cell, cell_str]=xlsread([EV_export_filepath '\' file_cells(sort_indx(i)).name]);
    %% layer file
    dat_layer=xlsread([EV_export_filepath '\' file_layers(sort_indx(i)).name]);
    
%% remove IVC overlap transects     
    if strcmp(file_cells(sort_indx(i)).name(1:4), US_Ship_ID)
        ship_name = 'US';
    else
        ship_name = 'CAN';
    end    
    transect = 1;
    if strcmp(ship_name, 'US') & ~isempty(ex_US_T) &  ~isempty(find(ex_US_T == TX_sort(i)))
        transect = [];
    end
    if strcmp(ship_name,'CAN') & ~isempty(ex_CAN_T) &  ~isempty(find(ex_CAN_T == TX_sort(i)))
        transect = [];
    end   

    if size(dat_layer,2) >= 4 & ~isempty(transect)
        T=TX_sort(i);                             % transect number
        n=size(dat,1);                            % number of VL intervals
        out.dat(tot_ind+[1:n],1)=T;               % Transect number
        out.dat(tot_ind+[1:n],3:4)=dat(:,3:4);    % VL_start, VL_end
        out.dat(tot_ind+[1:n],5:6)=dat(:,7:8);    % Lat & Lon
%         out.dat(tot_ind+1,5:6)=dat(2,7:8);        % the first numbers of lat and lon are always 999
%        disp(tot_ind+1)
    %% wrong latitude & longitude data, using adjacent data to provide an interpolated value
        ind=find(dat(:,7) == 999);
        indx1=find(diff(ind) > 1);
        if length(ind) > 1
            disp([T ind'])
        end
        if (~isempty(ind) & isempty(indx1)) | any(isnan(dat(1,7:8))) == 1
            %% single bad data
            fprintf('  n0 (value = 999 within region, latitude correction) = %d\n',length(ind))
            for j=length(ind):-1:1
                if dat(j,7) == 999 % jth latitude datum correction
                    if j == length(ind)
                        out.dat(tot_ind+j,5)=2*dat(j+1,7)-dat(j+2,7);
                    else
                        out.dat(tot_ind+j,5)=2*out.dat(tot_ind+j+1,5)-out.dat(tot_ind+j+2,5);
                    end
                end
                if dat(j,8) == 999  % jth longitude datum correction
                    if j == length(ind)
                        out.dat(tot_ind+j,6)=2*dat(j+1,8)-dat(j+2,8);
                    else
                        out.dat(tot_ind+j,6)=2*out.dat(tot_ind+j+1,6)-out.dat(tot_ind+j+2,6);
                    end
                end
            end
        end
        %% multiple bad data, using previous data to extrapolate to the current value
        if ~isempty(indx1)
            fprintf('  n1 = %d (latitude correction)\n',length(ind(indx1+1:end)))
            for j=ind(indx1)+1:ind(end)
                if dat(j,7) == 999 % jth latitude datum correction
                    if j == ind(indx1)+1
                        out.dat(tot_ind+j,5)=2*dat(j-1,7)-dat(j-2,7);
                    else
                        out.dat(tot_ind+j,5)=2*out.dat(tot_ind+j-1,7)-out.dat(tot_ind+j-2,7);
                    end
                end
                if dat(j,8) == 999  % jth longitude datum correction
                    if j == ind(indx1)+1
                        out.dat(tot_ind+j,6)=2*dat(j-1,8)-dat(j-2,8);
                    else
                        out.dat(tot_ind+j,6)=2*out.dat(tot_ind+j-1,8)-out.dat(tot_ind+j-2,8);
                    end
                end
            end
        end
        if size(dat,2) == 8                    % wrong format
            out.dat(tot_ind+[1:n],9:12)=nan;
            fprintf('wrong format in (interval)!\n')
        else
            out.dat(tot_ind+[1:n],11)=dat(:,9);        % bottom depth
        end
        
        %% find specified species
        ind_s=[];
        for ij=1:size(species_name,2)      % loop through species
            % species read from the cell excel file and conver to lower case
%             species_name0=lower(char(cell_str(:,3)));
%             species_name1=[' "' lower(char(species_name(ij))) '"'];
%             ind_s_ij=strmatch(species_name1,species_name0);
% % %  Newer Echoview exports "class" without double quotes  1-6-2021            
            species_name0=lower(char(cell_str(:,3)));
            species_name1_a=[' "' lower(char(species_name(ij))) '"'];
            species_name1_b=lower(char(species_name(ij)));
            ind_s_ij_a=strmatch(species_name1_a,species_name0);
            ind_s_ij_b=strmatch(species_name1_b,species_name0);
            if isempty(ind_s_ij_a) & isempty(ind_s_ij_b)
                ind_s_ij = [];
            elseif ~isempty(ind_s_ij_a) & isempty(ind_s_ij_b)
                ind_s_ij = ind_s_ij_a;
            elseif isempty(ind_s_ij_a) & ~isempty(ind_s_ij_b)
                ind_s_ij = ind_s_ij_b;
            elseif ~isempty(ind_s_ij_a) & ~isempty(ind_s_ij_b)
                fprintf('============== With & w/o quotes all work! ==============\n')
            end
            
            
            
            ind_s=[ind_s; ind_s_ij];
        end
        ind_s=ind_s-1;               % first line (row) is the title line
                                     % dat_cell (used later) and cell_str (i.e. species_name0) have one row
                                     % offset
        
        for j=1:n  % loop through all intervals in the file
            if isempty(dat_cell)
                layer_ind_all=[];
            else
                layer_ind_all=find(dat(j,2) == dat_cell(ind_s,5));    % compare VL interval
            end
            
            if isempty(layer_ind_all)
                out.dat(tot_ind+j,2)=999;
                out.dat(tot_ind+j,10:12)=0;
                out.dat(tot_ind+j,7)=1;   % will leads to zero biomass since corresponding NASC number is zero
            else
                region_id_all=unique(dat_cell(ind_s(layer_ind_all),1));     % distinct region id on the ith transect
                ind_layer=[];
                %% Region ID from Cell file
                for ij=1:length(region_id_all) % loop through region ID
                    % multiple layers or regions in the same VL interval (jth)
%                     if T == 49 & region_id_all(ij) == 29
%                         disp([j T region_id_all(ij)])
%                     end
% fprintf('i,j, reg = %d\t %d \t  %d\n', i, j, T_reg(i).reg)

                    if isempty(intersect(region_id_all(ij),T_reg(i).reg)) % the first cell with the same interval, same region on the same transect
                        T_reg(i).reg=[T_reg(i).reg region_id_all(ij)];
                        fprintf('     Multiple layers: Interval = %d,\t ij = %d,\t Reg = %d\n',dat(j,2),ij,region_id_all(ij));
                        if ij > 1
                           fprintf('*** same Interval = %d\t T = %d,\t Reg = %d\n',dat(j,2),T,region_id_all(ij));
                        end
                    end
                end
                layer_ind=layer_ind_all;
                layer_indx=dat_cell(ind_s(layer_ind),6);
                layers=dat_layer(layer_indx,3:4);                       % begin and end layer depth
                layer_mean_depth=(max(layers(:,2))+min(layers(:,1)))/2;
                layer_height=max(layers(:,2))-min(layers(:,1));
                region_ID=dat_cell(ind_s(layer_ind(1)),1);
                out.dat(tot_ind+j,2)=region_ID;                         % region ID
                haul_ind=find(T == T_reg_haul(:,1) & region_ID == T_reg_haul(:,2));
%                 if T == 95 & region_ID == 8
%                    disp(file_cells(sort_indx(i)).name)
%                    disp([j T region_id_all(ij)])
%                end

                if ~isempty(haul_ind)
                    stratum_ID_ind=find(haul_strata(:,4) == T_reg_haul(haul_ind,3));
                    stratum_ID=haul_strata(stratum_ID_ind,2);
                    out.dat(tot_ind+j,7)=stratum_ID;                    % stratum ID
                    out.dat(tot_ind+j,13)=T_reg_haul(haul_ind,3);       % haul number
                else
                    fprintf('% no appropriate haul corresponding to transect-region-hual is found!\n  T = %d,\t Reg = %d\n',T,region_ID);
                end
                if size(dat_cell,2) < 8
                    fprintf('wrong format in (cell)!\n')
                    out.dat(tot_ind+j,12)=0;
                else
                    NASC=nansum(dat_cell(ind_s(layer_ind),8));
                    out.dat(tot_ind+j,9)=layer_mean_depth;               % layer_mean_depth
                    out.dat(tot_ind+j,10)=layer_height;                  % layer_height
                    out.dat(tot_ind+j,12)=NASC;                          % NASC
%                     if T == 59 & region_ID == 3 & NASC >0
%                         fprintf('\t T = %d,\t region_ID=%d,\t Intervel = %d,\t VL (start) = %6.1f\t NASC =%6.1f\n',T,region_ID,dat(j,2),(dat(j,2)-1)/2,out.dat(tot_ind+j,12));
%                     end
                    if out.dat(tot_ind+j,12) > NASC_threshold
%                         if T == 11 & region_ID == 9 & dat(j,2) == 2451
%                             disp([i,j])
%                         end
                        fprintf('  Large NASC: Region_ID=%d,\t Intervel = %d,\t VL (start) = %6.1f\t NASC =%6.1f\n',region_ID,dat(j,2),(dat(j,2)-1)/2,NASC);
                        TX_high(high_NASC_ind) = T;
                        Reg_high(high_NASC_ind) = region_ID;
                        NASC_high(high_NASC_ind) = NASC;
                        high_NASC_ind = high_NASC_ind + 1;
                    end
                    %               fprintf('\t T = %d,\t region_ID=%d,\t Intervel = %d,\t NASC =%6.1f\n',T,region_ID,dat(j,2),NASC);
                end
                out.dat_detail(tot_ind+j).layer_ind=layer_ind;
                out.dat_detail(tot_ind+j).ind_s=ind_s;
                out.dat_detail(tot_ind+j).layer_indx=layer_indx;
                out.dat_detail(tot_ind+j).NASC=dat_cell(ind_s(layer_ind),8);
%                 if tot_ind+j == 653 
%                     disp(tot_ind+j)
%                 end
            end
        end
        % determine transect spacing
        if i >= 3
            ind1=find(out.dat(:,1) == TX_sort(i-2) & out.dat(:,5) < 60);  % remove '999' entry (lat < 60 N)         --> TX(i-2)
            ind2=find(out.dat(:,1) == T & out.dat(:,5) < 60);    % remove '999' entry  (lat < 60 N)                 --> TX(i)
            ind=find(out.dat(:,1) == TX_sort(i-1));              %                                                  --> TX(i-1)
            lat1=nanmean(out.dat(ind1,5));     % mean latitude of the current transect - 2
            lat2=nanmean(out.dat(ind2,5));     % mean latitude of the current transect
            if isnan(lat1) | isnan(lat2)       % begining or end survey transect
                out.dat(ind,8)=10;             % default number
            elseif abs(lat2-lat1) <= 2*max_spacing*1.1/30 &  max(out.dat(ind1,5))- min(out.dat(ind1,5)) < 1/6 &  max(out.dat(ind2,5))- min(out.dat(ind2,5)) < 1/6   % less than 15 nmi (0.5 * 30 nmi)  & paralell transects
                out.dat(ind,8)=abs(lat2-lat1)*30;    % transect (T-1) spacing = 1/2 of the spacing between two adjacent transects
            else
                out.dat(ind,8)=max_spacing;           % if transect spacing more than 15 nmi or zig-sap transects
            end
            if  i == nf                                % the last region
                ind=find(out.dat(:,1) == T & out.dat(:,5) < 60);    % remove '999' entry  (lat < 60 N)
                out.dat(ind,8)=max_spacing;
            end
        elseif i == 1                                 % the first transect
            ind=find(out.dat(:,1) == T & out.dat(:,5) < 60);
            out.dat(ind,8)=max_spacing;
        elseif i == 2                                 % the 2nd transect
            ind= out.dat(:,1) == T & out.dat(:,5) < 60;
            out.dat(ind,8)=max_spacing;
        end
        tot_ind=tot_ind+n;
    else
        if isempty(transect)
            fprintf('Duplicated transect: %s\n', file_cells(sort_indx(i)).name);
        else
            errordlg(sprintf('File %s is missing columns !!!',file_layers(sort_indx(i)).name),'ERROR DATA FILE');
        end
    end
end   % end of file loop
ind=find(out.dat(:,5) > 90);  % lat = 999 at the end of the file
out.dat(ind,:)=[]; 
if high_NASC_ind > 1
    out.dat_high = [TX_high(:) Reg_high(:) NASC_high(:)];
    header ={'Transect', 'Region', 'NASC'};
    xlswrite('NASC_check_output.xlsx', header, 1, 'A1')
    xlswrite('NASC_check_output.xlsx', out.dat_high, 1, 'A2')
end
return