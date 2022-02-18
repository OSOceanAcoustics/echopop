function      info_str=get_disp_str(disp_text_handle,para,show_indx)
% display variable values
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

switch show_indx
    case 1
        switch para.proc.source
            case 1
                DataSource_str='US';
            case 2
                DataSource_str='CAN';
            case 3
                DataSource_str='US & CAN';
        end
        switch para.proc.stratification_index
            case 1
               StatificationScheme_str='Hauls'; 
            case 2
               StatificationScheme_str='INPFS'; 
            case 3
               StatificationScheme_str='Kolmogorov-Smirnov'; 
            case 4
               StatificationScheme_str='Whole Survey Region'; 
            case 5
               StatificationScheme_str='Echo Energy Weighted'; 
        end
        switch para.proc.scale
            case 1
               EI_Region_str='Interval';
            case 2
               EI_Region_str='Judged Region';
            case 3
               EI_Region_str='Transect';
            case 4
               EI_Region_str='Leg';
            case 5
               EI_Region_str='Geographic Region';
            case 6
               EI_Region_str='Country';
        end
        info_str={['Raw Acoustcic Data Drive:  ' para.drive],'','Home directory:  ',['  ' para.home_dir], ...
            '',['Target Species Name:  ', para.target_species_name],'', ['Survey Year:  ' para.survey_year], ...
            '',['Data Source:  ' DataSource_str],'',['Statification Scheme:  ' StatificationScheme_str], ...
            '',['Start Transect #:  ' num2str(para.proc.start_transect)],'',['End Transect #:  ' num2str(para.proc.end_transect)], ...
            '',['Echo-Integration Scale :  ' EI_Region_str]};
    case 2
            info_str={['Species ID:  ' num2str(para.bio.species_code_ID)],['Number of Length Bins:  ' num2str(length(para.bio.hake_len_bin))], ...
            ['Number of Age Bins:  ' num2str(length(para.bio.hake_age_bin))],['Database:  ' para.bio.database_type], '', ...
            'Filenames:'};
            n=size(info_str,2);
            switch para.proc.source
            case   1  % US only
                info_str{n+1}=[' Trawl: ' para.bio.filename.trawl_US];
                info_str{n+2}=[' Gear: ' para.bio.filename.gear_US];
                info_str{n+3}=[' Catch: ' para.bio.filename.catch_US];
                info_str{n+4}=[' Length: ' para.bio.filename.length_US];
                info_str{n+5}=[' Specimem: ' para.bio.filename.specimen_US];               
            case   2  % CAN only
                info_str{n+1}=[' Trawl: ' para.bio.filename.trawl_CAN];
                info_str{n+2}=[' Gear: ' para.bio.filename.gear_CAN];
                info_str{n+3}=[' Catch: ' para.bio.filename.catch_CAN];
                info_str{n+4}=[' Length: ' para.bio.filename.length_CAN];
                info_str{n+5}=[' Specimem: ' para.bio.filename.specimen_CAN];
            case   3  % US & CAN
                info_str{n+1}=[' Trawl (US): ' para.bio.filename.trawl_US];
                info_str{n+2}=[' Trawl (CAN): ' para.bio.filename.trawl_CAN];
                info_str{n+3}='';
                info_str{n+4}=[' Gear (US): ' para.bio.filename.gear_US];
                info_str{n+5}=[' Gear (CAN: ' para.bio.filename.gear_CAN];
                info_str{n+6}='';
                info_str{n+7}=[' Catch (US): ' para.bio.filename.catch_US];
                info_str{n+8}=[' Catch (CAN): ' para.bio.filename.catch_CAN];
                info_str{n+9}='';
                info_str{n+10}=[' Length (US): ' para.bio.filename.length_US];
                info_str{n+11}=[' Length (CAN): ' para.bio.filename.length_CAN];
                info_str{n+12}='';
                info_str{n+13}=[' Specimem (US): ' para.bio.filename.specimen_US];
                info_str{n+14}=[' Specimem (CAN): ' para.bio.filename.specimen_CAN];
        end
    case 3
        info_str={['Available Frequencies:  ' num2str(para.acoust.freq0)  ' (kHz)'],'',['Processing Frequency:  ' num2str(para.acoust.freq)  ' (kHz)'], '', ...
            'Strata Filename:', para.acoust.filename.strata, '','NASC Filename:',para.acoust.filename.processed_data,'', ...
            ['Number of Biological Sampling Stations for TS:  ',num2str(para.acoust.TS_station_num)]};
end


set(disp_text_handle,'string',info_str,'HorizontalAlignment','left');
return