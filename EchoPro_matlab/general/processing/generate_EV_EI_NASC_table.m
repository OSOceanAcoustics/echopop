function    generate_EV_EI_NASC_table
% construct intervel-based echo integration NASC table
% based on the exported files from echoview
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global data para hdl

if isfield(para, 'platform_name') 
    if strcmp(para.platform_name, 'SD')  % SD data
        platform_ID = 2;
    else
        platform_ID = 1;                 % FSV data
    end
else
    platform_ID = 1;                     % FSV data
end
    
%% haul_strata filename
haul_strata_filename=para.acoust.filename.strata;
transect_reg_haul_filename=para.bio_acoust.filename.Transect_region_haul;

if get(hdl.proc.adv_proc.radio_output_filename,'value') == 1 & isfield(para.acoust.filename,'proc_EV_NASC') 
    output_full_filename=para.acoust.filename.proc_EV_NASC;
else
    disp('need to define an output filename!!')
    return
end

tic
if get(hdl.proc.adv_proc.radio_2yp,'value') == 1 & get(hdl.proc.adv_proc.radio_1yp,'value') ~= 1 & isfield(para.acoust.filename,'EV_export_input_filepath') 
    if get(hdl.proc.adv_proc.radio_mixed,'value') == 1 & ~isempty(get(hdl.proc.adv_proc.text_mixed_file_path,'string'))
        disp('Hake Age 2+ & Hake Mix exported from EchoView simultanenously ...')
        species_name={'Hake','Hake Mix'};
    else
        disp('Hake Age 2+ ...')
        species_name={'Hake'};
    end
    %%%%%%%%%%%%%%% process Echoview export files  %%%%%%%%%%%%%%%%%%%%%%%
    %% Transect_region_haul_filename (year 2+)
    EV_export_input_filepath=para.acoust.filename.EV_export_input_filepath;
    para.proc.hake=2;     % age2+
    if platform_ID == 1   % FSV
        out=construct_NASC_int_file_new(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    else                  % SD
        out=SD_construct_NASC_int_file(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    end
    dat2=out.dat;
else
    disp('---- No Hake Age 2+ processing !!')
end
toc

tic
if get(hdl.proc.adv_proc.radio_1yp,'value') == 1 & isfield(para.acoust.filename,'EV_export_input_filepath')
   if get(hdl.proc.adv_proc.radio_mixed,'value') == 1 & ~isempty(get(hdl.proc.adv_proc.text_mixed_file_path,'string'))
        disp('Hake Age 1+ & Hake Mix exported from EchoView simultanenously ...')
        species_name={'Age-1 Hake','Age-1 Hake Mix','Hake','Hake Mix'};  % also exclude Age-0
   else
       disp('Hake Age 1+ ...')
       species_name={'Age-1 Hake','Hake'};                               % also exclude Age-0
   end
    %% Transect_region_haul_filename (year 21)
    EV_export_input_filepath=para.acoust.filename.EV_export_input_filepath;
    para.proc.hake=1;     % age1+
    output_full_filename=para.acoust.filename.proc_EV_NASC;
    if platform_ID == 1      %FSV
        out=construct_NASC_int_file_new(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    else                  % SD
        out=SD_construct_NASC_int_file(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    end
    dat1=out.dat;
else
    disp('---- No Hake Age 1+ processing !!')
end
toc

tic
if get(hdl.proc.adv_proc.radio_1y_only,'value') == 1 & isfield(para.acoust.filename,'EV_export_input_filepath')
    disp('Hake Age 1 only ...')
    species_name={'Age-1 Hake', 'Age-1 Hake Mix'};
    %% Transect_region_haul_filename (year 1)
    EV_export_input_filepath=para.acoust.filename.EV_export_input_filepath;
    output_full_filename=para.acoust.filename.proc_EV_NASC;
    if platform_ID == 1   % FSV
        out=construct_NASC_int_file_new(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    else                  % SD
        out=SD_construct_NASC_int_file(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    end
    dat1=out.dat;
else
    disp('---- No Hake Age 1 ONLY processing !!')
end
toc

tic
if get(hdl.proc.adv_proc.radio_mixed,'value') == 1 & isfield(para.acoust.filename,'EV_export_mixed_filepath') & ~isempty(get(hdl.proc.adv_proc.text_mixed_file_path,'string')) ...
        &  get(hdl.proc.adv_proc.radio_2yp,'value') == 0 & get(hdl.proc.adv_proc.radio_1yp,'value') == 0
    disp('Hake mix only ...')
    species_name={'Hake Mix'};
    %% Transect_region_haul_filename (mix)
    EV_export_input_filepath=para.acoust.filename.EV_export_mixed_filepath;
    output_full_filename=para.acoust.filename.proc_EV_NASC;
    if platform_ID == 1    % FSV
        out=construct_NASC_int_file_new(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    else                   % SD
        out=SD_construct_NASC_int_file(species_name,haul_strata_filename,EV_export_input_filepath,transect_reg_haul_filename,para.proc.max_spacing);
    end
    dat_mix=out.dat;
else
    disp('---- No separate Hake Mix processing !!')
end
toc

dat=[];
if exist('dat2','var')
    dat=[dat; dat2];
end

if exist('dat1','var')
    dat=[dat; dat1];
end

toc
if exist('dat_mix','var')
    dat=[dat; dat_mix];
end

if isempty(dat)
    disp('Need to pick an input filepath !!!')
    return
end
%% excecute the save command
% cmd2=['xlswrite(' '''' output_full_filename ''',dat)'];
% disp(cmd2)
% eval(cmd)

var_name={'Transect','Region ID','VL start','VL end','Latitude','Longitude','Stratum','Spacing','Layer mean depth','Layer height','Bottom depth','NASC','Assigned haul'};
xlswrite(output_full_filename,var_name,1,'A1');
xlswrite(output_full_filename,dat,1,'A2');

para.acoust.filename.processed_data=output_full_filename;

return