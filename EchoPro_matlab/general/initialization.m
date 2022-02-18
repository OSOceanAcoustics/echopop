function  initialization(survey_year)
% initialization is to set all default parameters
%%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para data hdl

para.target_species_name='Hake';                % name of the target species
para.tasks.opr_indx=1;                          % 1 = Redefine Data Filenames, 2 = View Parameters, 3 = Data processing
                                                % 4 = Visualization, 5 = Report  
if nargin  == 0                                                
    para.survey_year='2021';
    para.proc.source=3;                         % 1 = US, 2 = CAN, 3 = All
else
    para.survey_year=survey_year;
    switch para.survey_year
        case {'1995','2005','2007'}
            para.proc.source=1;                 % 1 = US only
        case '2003'
            para.proc.source=2;                 % 2 = CAN only            
        otherwise
            para.proc.source=3;                 % 3 = Both US and CAN            
    end    
end
para.bio=[];

if str2num(para.survey_year) > 2009
    para.US_Ship_ID = 'V160';
elseif str2num(para.survey_year) > 2003 
    para.US_Ship_ID = 'V21';
else
    para.US_Ship_ID = 'V01';
end
%% processing parameters
% para.proc.scale=1;                              % 1 = interval; 2 = region, 3 = transect; 4 = Leg, 5 = Geographic region, 
%                                                 % 6 = country
if isfield(para, 'platform_name')
    if strcmp(para.platform_name, 'SD')  % SD data
    %%   default bio data type = Observer trawls
        para.bio_data_type = 3;                         % 1 = Acoustic survey; 2 = Bottom Trawl survey; 3 = Observer
    end
else
    %%   default bio data type = FSV trawls
    para.bio_data_type = 1;                         % 1 = Acoustic survey; 2 = Bottom Trawl survey; 3 = Observer
end
para.proc.age_data_status=1;                    % 0 = fake age data, 1 = actual age data, 2 = from age_vs_len
para.proc.exclude_age1=1;                       % 0 = include age 1 hake,   1 = exclude age 1 hake
para.bio.Lmin=15;                               % minimum length for age1+
para.proc.len_wgt_method=1;                     % 1 = using length PDF of whole region with no gender difference
                                                % 2 = using length PDF of whole region with gender difference  
                                                % 3 = using length PDF of each strtum
                                                % 4 = using averaged length of of each strtum
para.proc.survey_len_wgt_key_sex=0;             % 1 = both male & female over the whole survey to compute length-weight-key 
                                                % 0 = seperate male, female over the whole survey to compute length-weight-key
                                                % 
para.proc.max_spacing=10;                       % maximum transect spacing for biomass estimate
para.proc.int_for_zeros=5;                      % number of sample intervals for zeros along transect (nmi)
para.proc.disp_bins=120;                         % number of bins for histogram visualization
para.proc.bootstrap.limit=1;                    % set initial bootstrapping limit be unity
para.proc.bootstrap.cnt=0;                      % set initial bootstrapping count be zero
para.proc.hemisphere='NW';                      % north/south and west/east hemispheres
para.proc.end_index=9000;                       % end index for transects whose sample indces are less than this number and whose west extend will be added zeros
para.proc.add_zero_num=10;                      % number of samples added to the west extend of the transects
para.proc.extrapolation=0;                      % No extrapolation for kriging
para.proc.ordered_by_transects=1;               % flag of constructing NASC table based on the transect numbers
para.proc.WA_CAN_border=48.5;                   % latitute of WA/CAN border
para.proc.AK_CAN_border=55.4;                    % latitute of AK/CAN border
para.proc.JH_fac=0.75;                          % percent of Jolly-Hampton transects in each stratum

para.bio.species_code_ID=22500;                 % RACE & NWFSC adult hake Species ID
para.bio.hake_len_bin=2:2:80;                   % length sequence array: 1 - 80 cm 
para.bio.hake_age_bin=1:20;                     % age sequence array: year - age  1-20
para.bio.age1_min_len=10;                       % minimum length (cm) of age1 hake

%% acoustic data parameters
para.acoust.file_sys=1;                         % 1 = EK60, 2 = EK500,  3 = EK6
para.acoust.file_type=1;                        % 1 = processed, 2 = raw

%% process parameters
para.acoust.freq0=[18 38 70 120 200];              % all available frequencies for Shimada
para.acoust.freq=38000;                         % frequencies to be processed
para.acoust.bw=[11 7 7 7 7]*pi/180;             % beamwidth (radian)
[var,ind1,ind2]=intersect(floor(para.acoust.freq/1000),para.acoust.freq0);   % selected fequency index
para.acoust.freq_ind=ind2;                      % default = 2 --> 38 kHz 
para.acoust.sig_b_coef=10^(-6.8);               % coefficient of the sigma_bs = coef *L^2

%% default visualization parameters
para.visual.var1=5;                          % 1 = length, 2 = age, 3 = NASC, 4 =  number, 5 = biomass
para.visual.var2=2;                          % 1 = transect, 2 = Lat, 3 = strata, 4 = trawl, 5 = Length, 6 = Age, 7 = gender
                                             % 8 = survey region, 9 = lat/lon
para.visual.kriged=0;                        % display quantities: 1 = kriged quantities, 0 = non-kriged quantities
para.visual.disp_clr_fac=0.3;               % display color scale factor
                                             
%% kriging parameters
para.proc.kriging=1;                         % whether conduct kriging
para.proc.kriging_input=1;                   % 1 = biomass density; 2 = NASC;  3 = number density
para.proc.kriging_auto=1;                    % flag for automated kriging
para.proc.kriging_A0=2.5*2.5;                % base area of the grid cell
para.krig.shelf_break_Ref_lon=-124.78338;    % reference longitude
para.krig.x_res=0.02;                        % kriging resolution in lon (deg)
para.krig.y_res=0.05;                        % kriging resolution in lat (deg)
para.krig.loop=0;                            % kriging loop flag

para.proc.bio_info_disp=0;                   % flag for whether or not to display biological file information for debugging

load_proc_parameters(para.survey_year);      % popup FSV datafile names



