function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2003(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=3;    % southern most transect number
    tx1=105;  % northern most transect number
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[[tx0:94  117:119 98:105]+0.1 ];
    %% right (east) bound
    tx_r=[[tx0:96 117  98:105]+0.4];
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=105;  % east most transect number
    tx1=110;  % west most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[105.1 [106:108]+0.6 110.4];
    tx_u=[105.4 [106:108]+0.9 109.4];
    tx_out1=tx_l;
    tx_out2=tx_u;
elseif region == 3
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=117;     % southern most transect number
    tx1=109;    % northern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[109:117]+0.1;
    tx_r=[[109:116]+0.4 98.1 117.4];
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return