function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2013(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    tx1=117;  % northern most transect number
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[tx0:tx1]+0.1;
    %% right (east) bound
    tx_r=[tx0:tx1]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=117;    % west most transect number
    tx1=123;  % east most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[114.4 114.1 115.1 118.6 119.6 120.6 122.6 124.6];
    tx_u=[114.4 117.4 118.9 119.9 120.9 122.9 124.9];
    tx_out1=tx_l;
    tx_out2=tx_u;
else
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=101;    % southern most transect number
    tx1=125;    % northern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[101 137 135 133 130 127 126 125]+0.1;
    tx_r=[101.1 102.1 103.1 138.4 135.4 133.4 130.4 123.6 123.9];
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return