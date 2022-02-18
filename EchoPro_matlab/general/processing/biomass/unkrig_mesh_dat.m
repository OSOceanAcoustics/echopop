function data = unkrig_mesh_dat(data, para)
% unkrig_mesh_dat computes the biomass density at the center of each grid
% cell (mesh) based on interpolation operation using unkriged biomass
% density data
%   10/13/2017

%% original unkriged biomass data
lat0 = data.final.table.biomass(:,5);     % latitude
lon0 = data.final.table.biomass(:,6);     % longitude
nwgt0 = data.final.table.biomass(:,21);   % biomass density (kg/nm^2)

lon = data.in.US_CAN_Mesh(:,2);
lat = data.in.US_CAN_Mesh(:,1);

Gcell_size = 2.5;       % grid cell size (2.5 nm x 2.5 nm)

data_out = [];
tic
for i = 1:3         % three survey regions
    if para.proc.source == i | para.proc.source == 3
        nc = length(data.in.US_CAN_mesh_reg(i).ind); % number of grid cells
        switch i
            case {1, 3} % E-W transects
                for j = 1:nc
                    tx_all = find(data.final.table.biomass(:,1) >= data.in.US_CAN_mesh_reg(i).tx_min & data.final.table.biomass(:,1) <= data.in.US_CAN_mesh_reg(i).tx_max);
                    [lat_diff, ind] = min(abs(data.in.US_CAN_mesh_reg(i).lat(j) - data.final.table.biomass(tx_all,5))); % latitude of the closet transect
                    %                 if data.final.table.biomass(ind,1) == 124
                    %                     disp([i ind])
                    %                 end
                    data.in.US_CAN_mesh_reg(i).tx(j) = data.final.table.biomass(tx_all(ind),1);                  % corresponding transect of jth grid cell
                    Tind = find(data.final.table.biomass(:,1) == data.in.US_CAN_mesh_reg(i).tx(j));      % intervel indices of the chosen transect
                    VLcent = (data.final.table.biomass(Tind,3) +  data.final.table.biomass(Tind,4))/2;     % VL of the center of intervals
                    [lon_sorted, sort_ind ] = sort(data.final.table.biomass(Tind,6));
                    rand_num = 1e-6*randn(length(Tind),1);
                    VLgrid_cell = interp1(lon_sorted + rand_num, VLcent(sort_ind), data.in.US_CAN_mesh_reg(i).lon(j)); % interpolated VL at the grid cell center
                    %% find the intervals falling in the jth grid cell
                    int_ind = find(VLcent >= VLgrid_cell - Gcell_size/2 & VLcent <= VLgrid_cell + Gcell_size/2);
                    %% mean value of the biomass for each grid cell
                    data.in.US_CAN_mesh_reg(i).mean_var(j) = nanmean(data.final.table.biomass(Tind(sort_ind(int_ind)),15));
                    if isnan(data.in.US_CAN_mesh_reg(i).mean_var(j))
                        data.in.US_CAN_mesh_reg(i).mean_var(j) = 0;
                        %                     figure(4)
                        %                     plot(data.final.table.biomass(Tind,6), data.final.table.biomass(Tind,5), '.', data.in.US_CAN_mesh_reg(i).lon(j), data.in.US_CAN_mesh_reg(i).lat(j),'or')
                        %                     disp([i,j, length(Tind)])
                    end
                end
            case 2      % S-N transects
                for j = 1:nc
                    lat_ind=find(data.final.table.biomass(:,5) >= min(data.in.US_CAN_mesh_reg(i).lat));
                    [lon_diff, ind] = min(abs(data.in.US_CAN_mesh_reg(i).lon(j) - data.final.table.biomass(lat_ind,6))); % longitude of the closet transect
                    data.in.US_CAN_mesh_reg(i).tx(j) = data.final.table.biomass(lat_ind(ind),1);                  % corresponding transect of jth grid cell
                    Tind = find(data.final.table.biomass(:,1) == data.in.US_CAN_mesh_reg(i).tx(j));      % intervel indices of the chosen transect
                    VLcent = (data.final.table.biomass(Tind,3) +  data.final.table.biomass(Tind,4))/2;     % VL of the center of intervals
                    [lat_sorted, sort_ind ] = sort(data.final.table.biomass(Tind,5));
                    rand_num = 1e-6*randn(length(Tind),1);
                    VLgrid_cell = interp1(lat_sorted + rand_num, VLcent(sort_ind), data.in.US_CAN_mesh_reg(i).lat(j)); % interpolated VL at the grid cell center
                    
                    %% find the intervals falling in the jth grid cell
                    int_ind = find(VLcent >= VLgrid_cell - Gcell_size/2 & VLcent <= VLgrid_cell + Gcell_size/2);
                    %% mean value of the biomass for each grid cell
                    data.in.US_CAN_mesh_reg(i).mean_var(j) = nanmean(data.final.table.biomass(Tind(sort_ind(int_ind)),15));
                    if isnan(data.in.US_CAN_mesh_reg(i).mean_var(j))
                        data.in.US_CAN_mesh_reg(i).mean_var(j) = 0;
                    end
                end
        end
    end
end

switch para.proc.source
    case 1
        data.in.US_CAN_mesh_var = data.in.US_CAN_mesh_reg(1).mean_var(:);
        data.in.US_CAN_mesh_area = data.in.US_CAN_mesh_reg(1).area;
        data.in.US_CAN_total_mesh_area = sum(data.in.US_CAN_mesh_area);
    case 2
        data.in.US_CAN_mesh_var = data.in.US_CAN_mesh_reg(2).mean_var(:);
        data.in.US_CAN_mesh_area = data.in.US_CAN_mesh_reg(2).area;
        data.in.US_CAN_total_mesh_area = sum(data.in.US_CAN_mesh_area);
    case 3
        data.in.US_CAN_mesh_var = [data.in.US_CAN_mesh_reg(1).mean_var(:); data.in.US_CAN_mesh_reg(2).mean_var(:); data.in.US_CAN_mesh_reg(3).mean_var(:)];
        data.in.US_CAN_mesh_area = [data.in.US_CAN_mesh_reg(1).area; data.in.US_CAN_mesh_reg(2).area; data.in.US_CAN_mesh_reg(3).area];
        data.in.US_CAN_total_mesh_area = sum(data.in.US_CAN_mesh_area);
end

data.in.US_CAN_unkrig_dat = data.in.US_CAN_mesh_var;
toc
end   % end function

