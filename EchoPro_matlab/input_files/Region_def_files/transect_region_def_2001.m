function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2001(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    tx1=139;  % northern most transect number
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[tx0:85  151:-1:147 144 141 140 139]+0.1;
    %% right (east) bound
    tx_r=[tx0:85  151:-1:147 144 141 140 139]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=139;  % east most transect number
    tx1=134;  % west most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[139.1 [138:-1:136]+0.6 134.4];
    tx_u=[139.4 [138:-1:136]+0.9 135.4];
    tx_out1=tx_l;
    tx_out2=tx_u;
elseif region == 3
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=83;     % southern most transect number
    tx1=135;    % northern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[135:-1:128 147 126:-1:123 121:-1:119 83]+0.1;
    tx_r=[135:-1:128 126:-1:123 121:-1:119 83]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return