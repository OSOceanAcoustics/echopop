function   proc_raw_acoustic_data
%% process raw acoustic data used the region definition defined in 
%% the region def files exported from the EchoView
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

global para data


Final_out.Len_Age_Matrix_AcoustM=zeros(length(bio.para.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for male hake
Final_out.Len_Age_Matrix_AcoustF=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for female hake
Final_out.Len_Age_Matrix_AcoustALL=zeros(length(para.bio.hake_len_bin),length(para.bio.hake_age_bin)+1);  % Acoustic-weighted and stratum biological trawl normalized Abundance Matrix for all hake

para.grand_male_biomass=0;
para.grand_female_biomass=0;
para.grand_total_biomass=0;
para.leg_accum_male_biomass=0;
para.leg_accum_female_biomass=0;
para.leg_accum_total_biomass=0;
XnR_NASC=[];
VL_BiomassInt=[];
tmp=[];
start_region=1;
end_region=20;

for ii=1:length(para.proc.legs)
    bio_out=[];
    Leg=char(para.proc.legs(ii));
    fprintf('\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
    fprintf('*******************************  Processing Leg %s **************************************\n',Leg);
    fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
    t0=clock;
    for i=1:length(TX_list)                       % Transect loop: process transect data 
        clear Region
        para.cTX=['X',  num2str(TX_list(i))];
        para.cTx=['x',  num2str(TX_list(i))];
        fprintf('\n\n =============================  Processing %s ==============================  \n',para.cTX)
        para.VL=load_VL_data(para.Transect_VL_LatLon_dir,TX_list(i),para);
        fprintf('Get VL data = '); toc ;fprintf('\n')
        para.region_evr_filename=[para.region_root_evr_filename 'T' num2str(TX_list(i)) '.evr'];
        para.region_Transect_num=TX_list(i);
        Region=ReadRegionDef(para.acoust.filename.region_evr);
        if ~isstruct(Region)
            disp('Bad Region Filename !')
            break
        end
        XTb=sprintf('ST%d',TX_list(i));
        for j=1:length(Region)
            if strcmp(Region(j).name,XTb)
                break
            end
        end
        X(i).min_date=num2str2(Region(j).min_date,8);
        X(i).min_time=num2str2(Region(j).min_time,10);
        XTe=sprintf('ET%d',TX_list(i));
        for j=1:length(Region)
            if strcmp(Region(j).name,XTe)
                break
            end
        end     
        X(i).max_date=num2str2(Region(j).max_date,8);
        X(i).max_time=num2str2(Region(j).max_time,10);     
%   fprintf('Region Process:  %f\n',etime(clock,t0))
%% find correct hake regions on the current transect line
        for j=1:length(para.X_trawl_region)
            if strcmp(para.X_trawl_region(j).Xnum,para.cTX) | strcmp(para.X_trawl_region(j).Xnum,para.cTx)
                para.select_current_XR_indx=j;
                break
            end
        end
% process region data
        para.region=Region;
        hake_region=get_hake_bio_region_data(para);

%   fprintf('Bio Data:  %f\n',etime(clock,t0))
   % loop through selected hake region on the current Transect
        accum_len=[];
        accum_age=[];
        accum_age_len=[];
  %      clear trawl_no wgt_tot wgt_male wgt_female n_tot n_male n_female lsize slon slat
        trawl_no=[];
        wgt_tot=[];wgt_male=[];wgt_female=[];wgt_tot_ALL=[];
        nwgt_tot=[];nwgt_male=[];nwgt_female=[];nwgt_tot_ALL=[];
        n_tot=[];n_male=[];n_female=[];lsize=[];slon=[];slat=[];NASC=[];ABC=[];
        nn_tot=[];nn_male=[];nn_female=[];nn_int=[];n_int=[];
        n_tot_ping=[];wgt_tot_ping=[];
        N_tot=[];nN_tot=[];
        for j=start_region:min(end_region,length(hake_region))    % region loop
            NASC_int=[];ABC_int=[];VL_int_beg=[];VL_int_end=[];VL_dist=[];ntk=[];Lat_int=[];Lon_int=[];nntk=[];
            ABC_k=[];sd=[];nt=[];nd=[];dep=[];n_index=[];n_region=[];
            fprintf('\n\n -------------- Region = %d -----------------\n\n',Region(hake_region(j).reg_loc_no).n);
            region_out=get_selected_region_filename2009(para,Region,hake_region(j).reg_loc_no);
            filenames.sel_filenames=region_out.sel_filenames;
            filenames.filepath=region_out.filepath;
            data=get_selected_region_EK500_data(para,filenames,hake_region(j).reg_loc_no);   
            fprintf('Get EK60 & EK500 Raw Data = '); toc ;fprintf('\n')
 %     fprintf('Acoustic Data:  %f\n',etime(clock,t0))
            bad_data_region=get_bad_data_region(Region,Region(hake_region(j).reg_loc_no));          
            if ~isnan(bad_data_region.min_time)  % remove bad data region
                for k=1:length(bad_data_region.min_time)
                    min_dep=max(bad_data_region.min_dep(k),data.pings.dep(1));
                    max_dep=min(bad_data_region.max_dep(k), data.pings.dep(end));
                    dep_bad_region_ind=find( data.pings.dep >= min_dep & data.pings.dep <= max_dep);
                    min_time=max(bad_data_region.min_time(k),data.pings.time(1));
                    max_time=min(bad_data_region.max_time(k), data.pings.time(end));
                    time_bad_region_ind=find( data.pings.time >= min_time & data.pings.time <= max_time);  
                    data.pings.Sv(dep_bad_region_ind,time_bad_region_ind)=-999;
                end
            end
            dZ=data.calParms(1).sampleresolution;
            dep0=data.pings.range(1)+data.config.transducerdepth+dZ/2;
            [Ds,a12,a21] = m_idist(data.lon(2:end),data.lat(2:end),data.lon(1:end-1),data.lat(1:end-1));
            clear sm stot sv sa nt 
            ind_hake=[];
            for k=1:length(para.hake_length_sex)
                if para.hake_length_sex(k).trawl_no == para.X_trawl_region(para.select_current_XR_indx).region(j).trawl_no
                    ind_hake=para.X_trawl_region(para.select_current_XR_indx).region(j).trawl_no;
                    break
                end
            end
            if isempty(ind_hake)
                for k=1:length(para.hake_length_weight_sex_age)
                    if para.hake_length_weight_sex_age(k).trawl_no == para.X_trawl_region(para.select_current_XR_indx).region(j).trawl_no
                        ind_hake=para.X_trawl_region(para.select_current_XR_indx).region(j).trawl_no;
                        break
                    end
                end         
                if isempty(ind_hake)
                    fprintf('\n******************** Corresponding hake trawl not found!!! ********************\n\n');
                    data_flag=0;
                    break
                else
                    data_flag=2;
                end
            else
                data_flag=1;
            end
            haul_stratum_id = para.haul_strata(ind_hake);
            
    %% finc proportion coefficient for mix species
            mix_sa_ratio=1;
            if para.mix_region_flag == 1
                % find the proportion coefficient that matches all the IDs
                ind_mix=find(para.mix_trawls.TX_no == TX_list(i) & para.mix_trawls.region_ID == Region(hake_region(j).reg_loc_no).n  ...
                      & para.mix_trawls.stratum == haul_stratum_id);
                mix_sa_ratio=para.mix_trawls.species1_sa_proportion(ind_mix);
                if mix_sa_ratio < 0.01
                    mix_sa_ratio=100*mix_sa_ratio;   % if decimal number, change to percent
                end
            end
            np=size(data.pings(1).Sv,2);
            
     
           if para.VL.VL_flag ~= -1
 %%%%%%%%% Echo Integration for each interval %%%%%%%%
            para=find_VL_info(para,data);
  %        disp(para.VL.transect_spacing/1852)
            VL_ping_indx=[1 cumsum(para.VL.n_ping_interval)+1];
            npings=para.VL.npings;
            sig_b=para.strata(haul_stratum_id).sig_b;   % diff. backscattering cross section
            transect_spacing=para.VL.transect_spacing/1852;                       % transect spacing in nmi
            mean_transect_spacing=mean(transect_spacing);
            for k=1:length(npings)
                indx_k=VL_ping_indx(k):min(VL_ping_indx(k+1),np);
                Svk=data.pings.Sv(:,indx_k);
                Svk(Svk > -1000 & Svk < -69)= -999;
                ind=find(~isnan(Svk)== 1);                                   % thresholding at -70 dB
                sv=10.^(Svk(ind)/10);                                     % linear volume backscattering coef.
                ABC_int(k)=nansum(sv)*dZ/npings(k);
                NASC_int(k)=fac2*ABC_int(k);
                ntk(k)=round(mix_sa_ratio*transect_spacing(k)*para.interval*NASC_int(k)/sig_b);     % # hake/colume within the interval    
                Lat_int(k)=nanmean(data.lat(indx_k));
                Lon_int(k)=nanmean(data.lon(indx_k));
                VL_int_beg(k)=para.VL.VL(para.VL.start_VL_ind+k-2);
                VL_int_end(k)=para.VL.VL(para.VL.start_VL_ind+k-1);
                VL_dist(k)=VL_int_end(k)-VL_int_beg(k);
                nntk(k)=round(mix_sa_ratio*para.interval*NASC_int(k)/sig_b)/VL_dist(k);     % normalized # hake/colume within the interval and spacing #/area
            end
            nn_int(j)=round(sum(nntk));
            n_int(j)=round(sum(ntk));   
           else
               mean_transect_spacing=10;
           end  % end of VL file available loop
     % column sv integration
            for k=1:np     
                if k == 1
                    D=Ds(1);
                elseif k == size(data.pings(1).Sv,2)
                    D=Ds(end);
                else
                    D=0.5*(Ds(k-1)+Ds(k));
                end
                ind=find(~isnan(data.pings.Sv(:,k))== 1);                                   % thresholding at -70 dB
                sv=10.^(data.pings.Sv(ind,k)/10);                                     % linear volume backscattering coef.
                if data_flag == 1
                    sig_bs=para.strata(haul_stratum_id).sig_bs;   % diff. backscattering cross section
                else
                    sig_bs=para.strata(haul_stratum_id).sig_bs;   % diff. backscattering cross section
                end
                ABC_k(k)=nansum(sv)*dZ;                                                       % Areas Backscattering Coefficient (ABC)
                sd(k)=nansum(sv)*dZ*D;                                                % vertical integration of sv per unit distance
                                                                                      % (ping distance)
                nt(k)=round(mix_sa_ratio*mean_transect_spacing*1852*sd(k)/sig_bs);                                                   % # hake/colume per unit distance
                nd(k)=length(ind);                                                    % number of selected samples for the ping
                dep(k)=dep0+mean((ind-1/2)*dZ);
                if ~isempty(ind)
                  n_index(k,1:2)=[ind(1) ind(end)];
                  n_region(k,1:2)=([ind(1) ind(end)]-1)*dZ+region_out.min_dep;
                end
            end
            n_tot_ping(j)=nansum(nt);      % ping based acoustically estimated abundance


 %%%%%%%%%%%%% EchoPro - Echo Integration over fish aggregation region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum over pings to get quantities for each region
            Svj=data.pings.Sv(data.pings.Sv >= para.SvThreshold);
            Svj(Svj > para.Sv_min & Svj < para.SvThreshold)=para.Sv_min;
 %           Svj(Svj > para.Sv_min & Svj > -40)=nan;
            nsj=length(find(~isnan(data.pings.Sv) == 1));       % num of samples in the region
            Hmj=(nsj/np)*dZ;                                    % average layer thickness
            svj=10.^(Svj/10);                                   % dB -> linear
            sv_tot_mean_j=nansum(svj)/nsj;                      % mean sv: linear
            Sv_tot_mean_j=10*log10(sv_tot_mean_j);              % mean Sv (dB)
            ABC(j)=sv_tot_mean_j*Hmj;                           % Acoustic Backscattering Coefficient        m^2 / m^2
            NASC(j)=fac2*ABC(j);                                % Nautical Acoustic Scattering Coefficient   m^2 / mni^2
            
    %        sd_sum(j)=nansum(sd);                       % sum in linear space: sum = 10*log10( sum(10^(Sv/10)) )
    %        ABC_sum(j)=nansum(ABC)/size(data.pings(1).Sv,2);
            lsize(j)=nansum(Ds);               % fish school length-wise dimension in meters
            n_tot(j)=round(mix_sa_ratio*mean_transect_spacing*1852*ABC(j)/sig_bs*lsize(j));               % number of animals in the region estimated acoustically 
            
 %%%%% add values to Length-Age Abundance Matrix
            if para.integration_flag == 1 & ~isempty(ntk)
                N_tot(j)=n_int(j);
                nN_tot(j)=nn_int(j);                       % number of fish /nm^2
            else
                N_tot(j)=n_tot(j);
            end
 %% proportion the numbers
            N_len_M=N_tot(j)*para.strata(haul_stratum_id).Len_M_proportion;
            N_len_F=N_tot(j)*para.strata(haul_stratum_id).Len_F_proportion;
            N_len_age_M=N_tot(j)*para.strata(haul_stratum_id).Len_Age_M_proportion;
            N_len_age_F=N_tot(j)*para.strata(haul_stratum_id).Len_Age_F_proportion;
        
            Final_out.Len_Age_Matrix_AcoustM(:,end)=Final_out.Len_Age_Matrix_AcoustM(:,end)+round(N_len_M*para.strata(haul_stratum_id).Len_key_Mn(:));
            Final_out.Len_Age_Matrix_AcoustF(:,end)=Final_out.Len_Age_Matrix_AcoustF(:,end)+round(N_len_F*para.strata(haul_stratum_id).Len_key_Fn(:));                           
            Final_out.Len_Age_Matrix_AcoustALL(:,end)=Final_out.Len_Age_Matrix_AcoustALL(:,end) ...
                                                + round(N_len_M*para.strata(haul_stratum_id).Len_key_Mn(:) + N_len_F*para.strata(haul_stratum_id).Len_key_Fn(:));
            Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustM(:,1:end-1)+round(N_len_age_M*para.strata(haul_stratum_id).Len_Age_key_Mn(:,1:end));
            Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustF(:,1:end-1)+round(N_len_age_F*para.strata(haul_stratum_id).Len_Age_key_Fn(:,1:end));
            Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1)=Final_out.Len_Age_Matrix_AcoustALL(:,1:end-1) ...
                       +round(N_len_age_M*para.strata(haul_stratum_id).Len_Age_key_Mn(:,1:end)+N_len_age_F*para.strata(haul_stratum_id).Len_Age_key_Fn(:,1:end));

        %% number of hake included in the aggregation estimated based on acoustics  
            n_male(j)=N_len_M+N_len_age_M;
            n_female(j)=N_len_F+N_len_age_F;
            slon(j)=nanmean(data.lon);
            slat(j)=nanmean(data.lat);


fprintf('\n******************  EchoPro  *****************************\n')
fprintf('Sv_{mean} = %4.2f dB\n',Sv_tot_mean_j);
fprintf('Sv_{min} = %4.2f dB\n',min(Svj));
fprintf('Sv_{max} = %4.2f dB\n\n',max(Svj));

fprintf('NASC = %10.2f (m2/n. mi.2)\n',NASC(j));
fprintf(' ABC = %15.10f (m2/m2)\n\n',ABC(j));

fprintf('Height_{mean} = %5.2f m\n',Hmj);
%fprintf('Depth_{mean} = %4.2f m\n',nanmean(dep)-dZ);
fprintf('Depth_{mean} = %4.2f m\n',nansum(dep.*nd)/sum(nd));

fprintf('No. of Samples = %4d\n',nsj);
fprintf('No. of Pings = %4d\n',np);

%%%%%% normalized quantities
            nN_len_M=N_tot(j)*para.strata(haul_stratum_id).Len_M_proportion;
            nN_len_F=N_tot(j)*para.strata(haul_stratum_id).Len_F_proportion;
            nN_len_age_M=N_tot(j)*para.strata(haul_stratum_id).Len_Age_M_proportion;
            nN_len_age_F=N_tot(j)*para.strata(haul_stratum_id).Len_Age_F_proportion;
        
        %% number of hake included in the aggregation estimated based on acoustics  
            nn_male(j)=nN_len_M+nN_len_age_M;
            nn_female(j)=nN_len_F+nN_len_age_F;


%% region based estimates
%            idx_stratum=para.X_trawl_region(para.select_current_XR_indx).region(j).stratum;     % old program
            idx_stratum=haul_stratum_id;                                                        % new porgram 8-8-10
            if para.len_biomass == 1
                LaveM=para.strata(idx_stratum).LaveM;   % mean length of hake in assigned stratum corresponding to trawl kk
                len_wgt_M=para.len_wgt_all.reg_w0M*LaveM^para.len_wgt_all.reg_pM;
                LaveF=para.strata(idx_stratum).LaveF;   % mean length of hake in assigned stratum corresponding to trawl kk
                len_wgt_F=para.len_wgt_all.reg_w0F*LaveF^para.len_wgt_all.reg_pF;       
                len_wgt_ALL=para.len_wgt_all.reg_w0*para.strata(idx_stratum).Lave^para.len_wgt_all.reg_p;
            else
                len_wgt_M=nansum(para.strata(idx_stratum).Len_key_Mn.*para.len_wgt_all.reg_w0M.*para.hake_len_bin.^para.len_wgt_all.reg_pM);
                len_wgt_F=nansum(para.strata(idx_stratum).Len_key_Fn.*para.len_wgt_all.reg_w0F.*para.hake_len_bin.^para.len_wgt_all.reg_pF);       
                len_wgt_ALL=nansum(para.strata(idx_stratum).Len_key_Nn.*para.len_wgt_all.reg_w0.*para.hake_len_bin.^para.len_wgt_all.reg_p);
            end
            wgt_male(j)=n_male(j)*len_wgt_M*1e-6;                   % kg -> KMT
            wgt_female(j)=n_female(j)*len_wgt_F*1e-6;               % kg -> KMT
            wgt_tot(j)=wgt_male(j)+wgt_female(j);
            wgt_tot_ALL(j)=N_tot(j)*len_wgt_ALL*1e-6;
            trawl_no(j)=para.X_trawl_region(para.select_current_XR_indx).region(j).trawl_no;
            
       %%%%%% normalized quantities
            nwgt_male(j)=nn_male(j)*len_wgt_M*1e-6;                   % kg/nm^2 -> KMT/nm^2
            nwgt_female(j)=nn_female(j)*len_wgt_F*1e-6;               % kg/nm^2 -> KMT/nm^2
            nwgt_tot(j)=nwgt_male(j)+wgt_female(j);
            nwgt_tot_ALL(j)=nN_tot(j)*len_wgt_ALL*1e-6;

            
  %          XnR_NASC=[XnR_NASC ; TX_list(i) Region(hake_region(j).reg_loc_no).n NASC(j)];
            if  ~isempty(ntk)
            ntk_male=round(ntk*(para.strata(haul_stratum_id).Len_M_proportion+para.strata(haul_stratum_id).Len_Age_M_proportion));
            ntk_female=round(ntk*(para.strata(haul_stratum_id).Len_F_proportion+para.strata(haul_stratum_id).Len_Age_F_proportion));
            Wgt_male_int=ntk_male*len_wgt_M;           % in kg
            Wgt_female_int=ntk_female*len_wgt_F;           % in kg
%            Wgt_ALL_int=Wgt_female_int+Wgt_male_int;
            Wgt_ALL_int=ntk*len_wgt_ALL;  % in kg   % modified on 11/10/2010
          %%%%%% normalized quantities
            nntk_male=round(nntk*(para.strata(haul_stratum_id).Len_M_proportion+para.strata(haul_stratum_id).Len_Age_M_proportion));
            nntk_female=round(nntk*(para.strata(haul_stratum_id).Len_F_proportion+para.strata(haul_stratum_id).Len_Age_F_proportion));
            nWgt_male_int=nntk_male*len_wgt_M;           % in kg
            nWgt_female_int=nntk_female*len_wgt_F;           % in kg
  %          nWgt_ALL_int=nWgt_female_int+nWgt_male_int;  in kg
            nWgt_ALL_int=nntk*len_wgt_ALL;  % in kg   % modified on 10/7/2010
            
            VL_BiomassInt=[VL_BiomassInt; ...
              TX_list(i)*ones(length(VL_int_end),1) Region(hake_region(j).reg_loc_no).n*ones(length(VL_int_end),1) VL_int_beg(:) VL_int_end(:) Lat_int(:) Lon_int(:) NASC_int(:) ...
                                   ntk_male(:) ntk_female(:) ntk(:)  Wgt_male_int(:) Wgt_female_int(:) Wgt_ALL_int(:) ...
                                   nntk_male(:) nntk_female(:) nntk(:)  nWgt_male_int(:) nWgt_female_int(:) nWgt_ALL_int(:)];
            tmp=[tmp; TX_list(i)*ones(length(VL_int_end),1) Region(hake_region(j).reg_loc_no).n*ones(length(VL_int_end),1) VL_int_beg(:) VL_int_end(:) Lat_int(:) Lon_int(:)  ... 
                 haul_stratum_id*ones(length(VL_int_end),1) NASC_int(:) ntk(:) Wgt_ALL_int(:)];
            end
            wgt_tot_ping(j)=n_tot_ping(j)*len_wgt_ALL*1e-6;
            fprintf('\n Haul # = %d\t  Stratum # = %d\n',ind_hake, haul_stratum_id);
            fprintf('Total # of fish = %d\t, total Biomass = %6.3f (KMT)\n\n',N_tot(j), wgt_tot_ALL(j));
            XnR_NASC=[XnR_NASC ; TX_list(i) Region(hake_region(j).reg_loc_no).n NASC(j) N_tot(j) wgt_tot(j) wgt_tot_ALL(j)];
            disp(num2str([TX_list(i) Region(hake_region(j).reg_loc_no).n NASC(j) N_tot(j) wgt_tot(j) wgt_tot_ALL(j)]))
            disp(' ')

        end   % end hake region loop
   %     XnR_NASC=[XnR_NASC; NASC(:)];
        if data_flag ~= 0     % corresponding trawl found
            bio_out.trawl_no(i)=length(trawl_no);                          % number of trawls conducted within the transect
            bio_out.aggregation(i).trawl_no=trawl_no;
     %       bio_out.sd_sum(i)=nansum(sd_sum);
            bio_out.NASC(i)=fac2*nansum(ABC);
            bio_out.n_tot(i)=nansum(N_tot);
            bio_out.TX{i}=para.cTX;
            bio_out.n_male(i)=nansum(n_male);     
            bio_out.n_female(i)=nansum(n_female);     
            bio_out.TXaccum(i).length=accum_len;     
            bio_out.TXaccum(i).age=accum_age;
            bio_out.TXaccum(i).age_len=accum_age_len;     
            bio_out.wgt_tot(i)=nansum(wgt_tot_ALL);                     
            bio_out.wgt_male(i)=nansum(wgt_male);
            bio_out.wgt_female(i)=nansum(wgt_female);
            bio_out.aggregation(i).TX=para.cTX;
   %         bio_out.aggregation(i).sa_sum=sa_sum;
            bio_out.aggregation(i).lat=slat;
            bio_out.aggregation(i).lat=slon;
            bio_out.aggregation(i).lsize=lsize;
            bio_out.aggregation(i).wgt_male=wgt_male;
            bio_out.aggregation(i).wgt_female=wgt_female;
            bio_out.aggregation(i).wgt_tot=wgt_tot_ALL;
        end  
            
        fprintf('Current Transect Total Biomass (reg int) = %6.3f (KMT),\t  (ping int) = %6.3f\n',bio_out.wgt_tot(i),wgt_tot_ping(j))
        fprintf('Accumulative Total Biomass = %6.3f (KMT)\t Time Elapse = %f (sec)\n',sum(bio_out.wgt_tot),etime(clock,t0))

    end  % end Transect loop
    biomass(ii).male_biomass=sum(bio_out.wgt_male);
    biomass(ii).female_biomass=sum(bio_out.wgt_female);
    biomass(ii).tot_biomass=sum(bio_out.wgt_tot);
    para.leg_accum_male_biomass=para.leg_accum_male_biomass+biomass(ii).male_biomass;
    para.leg_accum_female_biomass=para.leg_accum_female_biomass+biomass(ii).female_biomass;
    para.leg_accum_total_biomass=para.leg_accum_total_biomass+biomass(ii).tot_biomass;
    if length(Leg0) > 1
        fprintf('---- Leg %s:  Estimated Biomass of Male: %5.2f (KMT),\t Female: %5.2f,\t Total: %5.2f (KMT) ----\n\n',...
            Leg,biomass(ii).male_biomass,biomass(ii).female_biomass,biomass(ii).tot_biomass)
%    else
%        fprintf('Estimated Biomass of Male: %5.2f (KMT),\t Female: %5.2f,\t Total: %5.2f (KMT) \n\n',male_biomass,female_biomass,tot_biomass)
    end
%     figure(3)
%     plot(TX_list,bio_out.wgt_tot,'s-','linewidth',1.5)
%     xlabel('Transect Number')
%     ylabel('Estimated Biomass (KMT)')
%     title(['Leg ' Leg])
%     grid
    toc
end  % end of Leg loop

para.grand_male_biomass=para.leg_accum_male_biomass;
para.grand_female_biomass=para.leg_accum_female_biomass;
para.grand_tot_biomass=para.leg_accum_total_biomass;

fprintf('\n ******************** All Legs:  Estimated Biomass *********** \n ******* Male: %5.2f (KMT)\n ******* Female: %5.2f (KMT)\n ******* Total: %5.2f (KMT) \n\n',...
    para.grand_male_biomass,para.grand_female_biomass,para.grand_tot_biomass)

% fid=fopen('US_hake2009_int_Biomass.csv','w');
% for i=1:size(VL_BiomassInt,1)
% %        VL_int_beg, VL_int_end, Lat_int,Lon_int, NASC_int,  ntk_male, ntk_female,  ntk,  Wgt_male_int, Wgt_female_int, Wgt_ALL_int 
%     fprintf(fid,'%8.2f,%8.2f,%10.4f,%10.4f,%8.2f,%12d,%12d,%12d,%12.2f,%12.2f,%12.2f\n', ...
%        VL_BiomassInt(i,1),VL_BiomassInt(i,2),VL_BiomassInt(i,3),VL_BiomassInt(i,4),VL_BiomassInt(i,5),VL_BiomassInt(i,6), ...
%        VL_BiomassInt(i,7),VL_BiomassInt(i,8),VL_BiomassInt(i,9),VL_BiomassInt(i,10),VL_BiomassInt(i,11));
% end
% fclose(fid);
toc