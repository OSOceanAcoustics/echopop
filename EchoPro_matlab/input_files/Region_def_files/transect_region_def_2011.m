function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2011(region)

if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
%     tx1=118;  % northern most transect number
    tx1=117;  % northern most transect number (new constructed export NASC file 11-04-2015)
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[tx0:77 83:tx1]+0.1;
    %% right (east) bound
    tx_r=[tx0:tx1]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
%     tx0=118;    % east most transect number
    tx0=117;    % east most transect number  (new constructed export NASC file 11-04-2015)
    tx1=124;    % west most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
%     tx_l=[114.4 114.1 118.1 120.6 122.6 124.6];
%     tx_u=[114.4 118.4 120.9 122.9 124.9];
    tx_l=[114.4 114.1 117.1 120.6 122.6 124.6];  % new constructed export NASC file 11-04-2015
    tx_u=[114.4 117.4 120.9 122.9 124.9];        % new constructed export NASC file 11-04-2015
    tx_out1=tx_l;
    tx_out2=tx_u;
else
%% region 3: paralell transects to latitudes west of QC IS
    tx0=106;    % southern most transect number
    tx1=128;    % northern most transect number
%% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[105.1 138.1 136.1 132.1 130.1 128.1];
    tx_r=[105.1 144.4 138.4 136.4 134.4 132.4 124.6 122.6 128.4];
    tx_out1=tx_l;
    tx_out2=tx_r;
end
return