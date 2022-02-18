function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2019(region)

global para

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    if para.proc.source == 1
        tx1=86;  % northern most transect number (US only)
    else
        tx1=119;  % northern most transect number (US & CAN)
    end
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
    tx0=121;    % west most transect number
    tx1=127;    % east most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[tx0:tx1] + 0.6;
    tx_u=[tx0:tx1] + 0.9;
    tx_out1=tx_l;
    tx_out2=tx_u;
else
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=129;    % northern most transect number
    tx1=145;    % southern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[tx0:tx1]+0.1;
    tx_r=[tx0:tx1]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return