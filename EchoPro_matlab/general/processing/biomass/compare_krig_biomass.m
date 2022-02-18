function   data = compare_krig_biomass(data)
% compare kriged and unkriged biomass estimates on each transect

uniq_transect = unique(data.final.table.biomass(:,1));

nt = length(uniq_transect);

plot_comparison_all_transect = 0;
plot_comparison_specified_transect = 0;

i0 = 0;
for i = 1:nt
    j = 1;
    while isempty(intersect(data.in.US_CAN_mesh_reg(j).tx, uniq_transect(i)))
        j = j +1;
        if j >= 3, break,end
    end
    %% mesh grid-based analysis
    ind = find(data.in.US_CAN_mesh_reg(j).tx == uniq_transect(i));
    data.in.US_CAN_mesh_tx_biomass(i) = nansum(data.in.US_CAN_mesh_reg(j).mean_var(ind));
    data.in.US_CAN_mesh_tx_area(i) = nansum(data.in.US_CAN_mesh_reg(j).area(ind));
    ind1 = find(data.final.table.biomass(:,1) ==  uniq_transect(i));
    data.in.US_CAN_tx_biomass(i) = nansum(data.final.table.biomass(ind1,15));
    if plot_comparison_specified_transect == 1
        fprintf('Reg# = %d\t TX = %d\t  # of Gcell = %d\n', j, uniq_transect(i), length(ind))
        if uniq_transect(i) == 12
            figure(2);plot(data.in.US_CAN_mesh_reg(j).lon(ind), data.in.US_CAN_mesh_reg(j).lat(ind), 'o-b');hold on
        end
    end
    %% transect-based analysis
    int_dist = data.final.table.biomass(ind1,4) - data.final.table.biomass(ind1,3);  % interval distance nm
    if size(data.final.table.biomass,2) > 21   % some information was missing before 2003  - 1/8/2017
        data.in.US_CAN_tx_area(i) = nansum(data.final.table.biomass(ind1,24).*int_dist); % transect spacing x interval distance
    end
end

if plot_comparison_all_transect == 1
    figure(5)
    plot(uniq_transect, data.in.US_CAN_tx_biomass*1e-6, '.-b', uniq_transect, data.in.US_CAN_mesh_tx_biomass*1e-6, 'o-r')
    ylabel('Biomass (kmt)')
    xlabel('Transect')
    grid
    if size(data.final.table.biomass,2) > 21   % some information was missing before 2003  - 1/8/2017
        figure(6)
        plot(uniq_transect, data.in.US_CAN_tx_area, '.-b', uniq_transect, data.in.US_CAN_mesh_tx_area, 'o-r')
        ylabel('Area (nm^2)')
        xlabel('Transect')
        grid
    end
end

end
