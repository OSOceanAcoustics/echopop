function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2012(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    tx1=107;  % northern most transect number
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
    tx0=107;    % west most transect number
    tx1=117;  % east most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[105.4 105.1 111.6 113.6 115.6 117.6];
    tx_u=[105.4 111.9 113.9 115.9 117.9];
    tx_out1=tx_l;
    tx_out2=tx_u;
else
%% region 3: paralell transects to latitudes west of QC IS
    tx0=97;    % southern most transect number
    tx1=120.1;   % northern most transect number
%% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[97 135 131 129 127 125 123 121 120]+0.1;
    tx_r=[97.1 99.1 133.4 131.4  129.4 128.4 126.4 117.6 117.9 119.4];    tx_out1=tx_l;
    tx_out2=tx_r;
end
return