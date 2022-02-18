function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_1995(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    tx1=152;  % northern most transect number
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[tx0 2 4:78 80:2:4 89:99 101:105 120 121 124:127 129:2:133]+0.1;
    %% right (east) bound
    tx_r=[tx0 2 4:88 91:105 152 149 145 143 142 140 135 134 133]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=116;  % east most transect number
    tx1=119;  % west most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[118 117 119 116] + 0.6;
    tx_u=[[118 117 119 116] + 0.9];
    tx_out1=tx_l;
    tx_out2=tx_u;
elseif region == 3
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=106;     % southern most transect number
    tx1=128;    % northern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[116.6 119.6 [115:-1:106]+0.1 126.1];
    tx_r=[119.6 [115:-1:106]+0.4 128.1];
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return