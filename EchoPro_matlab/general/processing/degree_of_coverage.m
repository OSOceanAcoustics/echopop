function    A = degree_of_coverage(data, para)
% para = parametrers structure from EchoPro
% data = output structure from EchoPro

% Dezhang Chu, NOAA, NMFS, NWFSC 
% 2-15-2018

d = data.final.table.biomass(:,4) - data.final.table.biomass(:,3); % VLend - VLstart  in nmi
D = nansum(d);          % total transect length in nmi


%% set niminal grid_cell area corresponding to the grid_cell file
if ~isempty(findstr(data.in.filename.grid_cell,'2_5'))
    para.proc.kriging_A0=2.5*2.5;
else ~isempty(findstr(data.in.filename.grid_cell,'1_25'));
    para.proc.kriging_A0=1.25*1.25;
end
A0=para.proc.kriging_A0;   % nominal area of each grid

%% load US & CAN Mesh files off the West Coast of US/CAN 
d=xlsread(data.in.filename.grid_cell);

if para.krig.proc_opt == 1      
%% whether to exclude the extrapolated region
% para.proc.extrapolation=1;
    if para.proc.extrapolation == 0
        %% find the kriged region without extrapolation
        fprintf(' Determine the kriging region that is bounded by the transect...\n')
        cmd=['[sel_grid_indx, data]=krig_transect_region(para,data);'];
        eval(cmd)
    else
        sel_grid_indx=1:size(d,1);
    end
    data.in.sel_grid_indx=sel_grid_indx;   % selected grid index of the extrapolation removed kriging grids
    Area=A0*d(sel_grid_indx,5);                % actual area of each mesh grid
else
%% after kriging process is called
    Area=data.in.US_CAN_Mesh_area;
end

A = D/sqrt(nansum(Area));

return