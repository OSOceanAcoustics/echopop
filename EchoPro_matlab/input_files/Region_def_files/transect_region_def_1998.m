function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_1998(region)

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    tx1=222;  % northern most transect number
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[[tx0:60 81:103 249 252]+0.1 253.6 253.9 [254 248 245 239 238 236 233 231 227]+0.1 226.6 211.6];
    %% right (east) bound
    tx_r=[[tx0:60 70:80 99 265 262 260 271 274 255 246 238 236 233 231 229]+0.4  226.9 223.9 222.4 211.9];
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=224;  % east most transect number
    tx1=319;  % west most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[[226 211 209 207 205 202]+0.6 320.4];
    tx_u=[224.6 [211 217 215 213]+0.9 212.1 325,4];
    tx_out1=tx_l;
    tx_out2=tx_u;
elseif region == 3
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=103;     % southern most transect number
    tx1=368;    % northern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[368.6 367.4 [365.6 363 360 358 356 353]+0.6 [351 347]+0.1 343.6 [340 338 331]+0.1 330.6 328.6 ...
         [325 322 320 318 316 313:-2:303 300 298 294]+0.1 293.6 [286 284 282 103]+0.1];
    tx_r=[368.9 366.4 364.9 [360 359 357 354 352 347 345 344 342 340 336 334]+0.4 330.9 [329 327 324]+0.4 212.1 322.4 321.6 ...
        [318 316 313:-2:307 304]+0.4 302.9 300.4 298.4 297.9 [295 292:-2:282]+0.4];
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return