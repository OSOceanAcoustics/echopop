function plot_results
%% plotting the chosen variables 
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global hdl para data

for i=1:length(para.proc.disp_method)    
    switch para.proc.disp_method(i)
        case 1  % 2D curve
            figure(1)
            clf
            %% remove NaN's
            ind=find(~isnan(data.visual.var) == 1);
            data.visual.var=data.visual.var(ind);
            data.visual.func=data.visual.func(ind,:);
            var1=unique(round(data.visual.var));
            if get(hdl.vis.radio_visual_var_transect,'value') ~= 1 & get(hdl.vis.radio_visual_var_trawl,'value') ~= 1 & (length(var1) > para.proc.disp_bins | get(hdl.vis.radio_visual_var_lat,'value') == 1) 
                inc=(max(data.visual.var)-min(data.visual.var))/para.proc.disp_bins;
                min_var=max(0,min(data.visual.var)-inc);
                max_var=max(data.visual.var)+inc;
                var1=linspace(min_var,max_var,para.proc.disp_bins);
            end
            inc=[var1(2)-var1(1)];
            var1=var1(1):inc:var1(end);
            for j=1:length(var1)
                if j == 1
                    ind=find(data.visual.var <= var1(j)+inc/2);
                elseif j == length(var1)
                    ind=find(data.visual.var > var1(j)-inc/2);
                else
                    ind=find(data.visual.var > var1(j-1)+inc/2 & data.visual.var <= var1(j+1)-inc/2);
                end
                if length(ind) == 1
                    func(j,1:size(data.visual.func,2))=data.visual.func(ind,:);
                    sd(j,1:size(data.visual.func,2))=data.visual.func(ind,:);
                else
                    if get(hdl.vis.radio_visual_func_number,'value') == 1 | get(hdl.vis.radio_visual_func_biomass,'value') ==1
                        func(j,1:size(data.visual.func,2))=nansum(data.visual.func(ind,:));
                    else
                        func(j,1:size(data.visual.func,2))=nanmean(data.visual.func(ind,:));
                    end
                    sd(j,1:size(data.visual.func,2))=nanstd(data.visual.func(ind,:));
                end
                len(j)=length(ind);
            end
            data.visual.var1=var1(:);
            if size(data.visual.func,2) ==1
                data.visual.func1=func(:);
            else
                data.visual.func1=func;
            end
            data.visual.var1_str=data.visual.var_str;
            data.visual.func1_str=data.visual.func_str;
            if get(hdl.vis.radio_visual_var_weight,'value') == 1 
                if get(hdl.vis.radio_visual_func_length,'value') == 1
                   ind=find((~isnan(data.bio.len_wgt_all.len) & ~isnan(data.bio.len_wgt_all.wgt)) == 1);
                   log_wgt=log(data.bio.len_wgt_all.wgt(ind));
                   log_len=log(data.bio.len_wgt_all.len(ind));
                   p=polyfit(log_len,log_wgt,1);
                   yfit=exp(p(2))*para.bio.hake_len_bin.^p(1);
                   plot(data.bio.len_wgt_all.len,data.bio.len_wgt_all.wgt,'.',para.bio.hake_len_bin,data.bio.len_wgt_ALL,'-g',para.bio.hake_len_bin,yfit,'-r','linewidth',2);
                   legend('Data','Semi-fitted Curve', 'fitted Curve','location','northwest')
                   h=text(0.05,0.70,sprintf('W = %6.3e L(cm)^{%4.2f} (kg)',data.bio.len_wgt_all.reg_w0,data.bio.len_wgt_all.reg_p),'sc');
                   set(h,'FontWeight','bold')
                else
                   plot(data.bio.len_wgt_all.age,data.bio.len_wgt_all.wgt,'.',data.visual.var1,data.visual.func1,'-r','linewidth',2);                   
                end
            else
                plot(data.visual.var1,data.visual.func1,'.-','linewidth',2);
                axis([min(var1)-inc max(var1)+inc 0 1.2*max(max(func))])
            end
            xlabel(data.visual.var1_str,'fontsize',20,'fontweight','bold');
            ylabel(data.visual.func1_str,'fontsize',20,'fontweight','bold');
            if size(data.visual.func1,2) > 1
                legend('Male','Female','ALL')
            end
        case 2  % 2D histogram
            if i == 1
                %% remove NaN's
                if get(hdl.vis.radio_visual_var_age,'value') ~= 1
                    ind=find(~isnan(data.visual.var) == 1);
                else
                    ind=1:length(data.visual.var);
                    ind(end)=length(ind);
                end
                if get(hdl.vis.radio_visual_var_survey_region,'value') ~= 1 
                    data.visual.var=data.visual.var(ind);
                    data.visual.func=data.visual.func(ind,:);
                    var1=unique(round(data.visual.var));
                    if get(hdl.vis.radio_visual_var_age,'value') == 1
                       var1(end)=length(var1); 
                    end
                    if get(hdl.vis.radio_visual_var_transect,'value') ~= 1 & get(hdl.vis.radio_visual_var_trawl,'value') ~= 1 ...
                            & get(hdl.vis.radio_visual_var_length,'value') ~= 1 & get(hdl.vis.radio_visual_var_age,'value') ~= 1 & length(var1) > para.proc.disp_bins
                        inc=(max(data.visual.var)-min(data.visual.var))/para.proc.disp_bins;
                        min_var=max(0,min(data.visual.var)-inc);
                        max_var=max(data.visual.var)+inc;
                        var1=linspace(min_var,max_var,para.proc.disp_bins);
                    end
                    inc=var1(2)-var1(1);
                    for j=1:length(var1)
                        if j == 1
                            ind=find(data.visual.var <= var1(j));
                        elseif j == length(var1)
                            if get(hdl.vis.radio_visual_var_age,'value') ~= 1
                               ind=find(data.visual.var >= var1(j));
                            else
                               ind=length(data.visual.var);
                            end
                        else
                            ind=find(data.visual.var > var1(j-1) & data.visual.var < var1(j+1));
                        end                      
                        if length(ind) == 1
                            func(j,1:size(data.visual.func,2))=data.visual.func(ind,:);
                            sd(j,1:size(data.visual.func,2))=data.visual.func(ind,:);
                        else
                            if get(hdl.vis.radio_visual_func_number,'value') == 1 | get(hdl.vis.radio_visual_func_biomass,'value') ==1
                                func(j,1:size(data.visual.func,2))=nansum(data.visual.func(ind,:));
                            else
                                func(j,1:size(data.visual.func,2))=nanmean(data.visual.func(ind,:));
                            end
                            sd(j,1:size(data.visual.func,2))=nanstd(data.visual.func(ind,:));
                        end
                    end
                else   % survey region
                    ind=find(data.visual.func(:,end) > 0);
                    if get(hdl.vis.radio_visual_biological,'value') == 1
                        var0=data.visual.func(ind,:);
                    else
                        var0=data.visual.func(ind,:)*1e-3;
                    end
                    if get(hdl.vis.radio_visual_func_length,'value') == 1 & get(hdl.vis.radio_visual_var_survey_region,'value') == 1
                        [func,var1]=hist(var0,1:max(para.bio.hake_len_bin));
                    elseif get(hdl.vis.radio_visual_func_age,'value') == 1 & get(hdl.vis.radio_visual_var_survey_region,'value') == 1
                        [func,var1]=hist(var0,1:max(para.bio.hake_age_bin));
                    else
                        [func,var1]=hist(var0,para.proc.disp_bins);
                    end
                    inc=var1(2)-var1(1);
                    data.visual.func_str='Counts';
                end
                data.visual.var1=var1(:);
                if size(func,2) >  size(func,1)
                    data.visual.func1=func';
                else
                    data.visual.func1=func;
                end
                data.visual.var1_str=data.visual.var_str;
                data.visual.func1_str=data.visual.func_str;
            end
            figure(2)
            clf
            gh=bar(data.visual.var1,data.visual.func1);
            axis([min(var1)-inc max(data.visual.var1)+inc 0 1.2*max(max(data.visual.func1))])
            grid on
            xlabel(data.visual.var1_str,'fontsize',20,'fontweight','bold');
            ylabel(data.visual.func1_str,'fontsize',20,'fontweight','bold');
            if size(data.visual.func1,2) > 1
                legend('Male','Female','ALL')
            end
            set(gca,'fontsize',14)
            title(para.survey_year,'fontsize',20)
            set(gh,'FaceColor','r')
            if strcmp(data.visual.func_str,'Biomass(tons)') == 1
                gh_text=text(0.1,0.9,sprintf('Total Biomass = %5.3f (mmt)',sum(data.visual.func1(:,3))/1e6),'sc','fontsize',14,'fontweight', 'bold');
                set(gca,'ylim',[0 1.5e5])
            end
        
            if get(hdl.vis.radio_visual_var_age,'value') == 1 & get(hdl.vis.radio_visual_func_number,'value') == 1
                xtick0_lable=get(gca,'xticklabel');
                set(gca,'xtick',[0:2:20 21]);
                xtick0_lable=get(gca,'xticklabel');
                new_xtick_label=[xtick0_lable(1:end-1,:)];
                ntick=size(new_xtick_label,1);
                new_xtick_label(ntick+1,1:7)='Un-Aged';
                set(gca,'xticklabel',new_xtick_label);
                if size(data.visual.func,2) > 1
                    legend('Male','Female','ALL',2)
                end
            elseif get(hdl.vis.radio_visual_var_gender,'value') == 1
                new_xtick_label=[' M ';' F ';'ALL'];
                set(gca,'xticklabel',new_xtick_label);
            end
            if get(hdl.vis.radio_visual_var_survey_region,'value') == 1 & i == 1
                %% Mean & s.d.
                var_mean=sum(repmat(data.visual.var1,1,size(data.visual.func1,2)).*data.visual.func1)./sum(data.visual.func1);
                var_std=sqrt(sum(repmat(data.visual.var1.^2,1,size(data.visual.func1,2)).*data.visual.func1)./sum(data.visual.func1)-var_mean.^2);
                switch data.visual.var1_str
                    case 'Biomass(mt)'
                        hmn=text(0.1,0.9,sprintf('Mean (all) = %4.1f (cm)',var_mean(3)),'sc','fontweight','bold');
                        hsd=text(0.1,0.8,sprintf('s.d. (all) = %4.1f (cm)',var_std(3)),'sc','fontweight','bold');
                    case 'Length'
                        hmn=text(0.1,0.9,sprintf('Mean (all) = %4.1f (cm)',var_mean),'sc','fontweight','bold');
                        hsd=text(0.1,0.8,sprintf('s.d. (all) = %4.1f (cm)',var_std),'sc','fontweight','bold');
                    otherwise
                        disp('not available!')
                end
            end
        case 3  % 2D pcolor
            lat=data.visual.var(:,1);
            lon=data.visual.var(:,2);
            func=data.visual.func(:,end);    % no sex differentiation
            fac=para.visual.disp_clr_fac;   % display scale factor
            display_NASC_number_biomass(lat,lon,func,data.visual.func_str,fac);
            data.visual.var1=lat;
            data.visual.var2=lon;
            data.visual.var1_str='Latitude';
            data.visual.var2_str='Longitude';
            data.visual.func1=data.visual.func;
            data.visual.func1_str=data.visual.func_str;
        case 4  % 2D images
            n=size(data.visual.func,1);
            n1=n/3;
            if get(hdl.vis.radio_visual_func_number,'value') == 1
               title_str={'Abundance (Aged Male)','Abundance (Aged Female)','Abundance (Aged ALL)'};
                data.visual.func1_str='Abundance';
            else
               title_str={'Biomass (Male)','Biomass (Female)','Biomass (ALL)'};
                data.visual.func1_str='Biomass (kg)';
            end
            num_wgt_flag=get(hdl.vis.radio_visual_func_number,'value');
            for j=1:3   % M, F, ALL
                ind1=((j-1)*n1+1):(j*n1);
                func=data.visual.func(ind1,1:(end-num_wgt_flag));
                var2=data.visual.var2(1:(end-num_wgt_flag));
                var1=data.visual.var1;
                figure(3+j)
                clf
                imagesc(var2,var1,func);
                xlabel(data.visual.var1_str,'fontsize',20,'fontweight','bold')
                ylabel(data.visual.var2_str,'fontsize',20,'fontweight','bold')
                title(title_str(j),'fontsize',20,'fontweight','bold')
                colorbar
                fac=0.05;   % display scale factor
                caxis([0 fac*max(max(func))])
                func2(1:length(ind1),j)=data.visual.func(ind1,end);
            end     
            if num_wgt_flag == 1 % abundance or biomass
                figure(7)
                clf
                bar(var1,func2*1e-3)
                xlabel(data.visual.var1_str,'fontsize',20,'fontweight','bold')
                ylabel('Abundance (x1000)','fontsize',20,'fontweight','bold')
                title('Un-aged','fontsize',20,'fontweight','bold')
                legend('Male','Female','ALL')
            end
            data.visual.var1=var1;
            data.visual.var2=var2;
            data.visual.func1=data.visual.func(:,1:(end-num_wgt_flag));
            data.visual.func2=func2;
            data.visual.func2_str=data.visual.func_str;
        case 5  % 3D histogram (with gender specific)
            n=size(data.visual.func,1);
            n1=n/3;
            if get(hdl.vis.radio_visual_func_number,'value') == 1
               title_str={'Abundance (Aged Male)','Abundance (Aged Female)','Abundance (Aged ALL)'};
               zlabels='Abundance';
            else
               title_str={'Biomass (Male)','Biomass (Female)','Biomass (ALL)'};
               zlabels='Biomass (tons)';
            end
            if  i == 1
                fig0=6;
            else
                fig0=7-get(hdl.vis.radio_visual_func_biomass,'value');
            end         
            var1=data.visual.var2;       % age bins
            var2=data.visual.var1;       % length bins
            for j=1:3   % M, F, ALL
                ind1=((j-1)*n1+1):(j*n1);
                func=data.visual.func(ind1,1:end-1);
                figure(fig0+j)
                clf
                bar3(var2,func,0.75);
                hx=xlabel(data.visual.var2_str,'fontsize',20,'fontweight','bold');
                hy=ylabel(data.visual.var1_str,'fontsize',20,'fontweight','bold');
                hz=zlabel(zlabels,'fontsize',20,'fontweight','bold');
                ht=title(title_str(j),'fontsize',20,'fontweight','bold');
                set(hy,'Rotation',-35);
                set(hx,'Rotation',25);
            end    
            data.visual.var1=var1;
            data.visual.var2=var2;            
        case 6  % 3D-histogram (no gender specific)
            figure(11)
            bar3(data.visual.var,data.visual.func,0.75);
            hx=xlabel(data.visual.func_str,'fontsize',20,'fontweight','bold');
            hy=ylabel(data.visual.var_str,'fontsize',20,'fontweight','bold');
            if get(hdl.vis.radio_visual_func_biomass,'value') == 1
               hz=zlabel('Biomass (kg)','fontsize',20,'fontweight','bold');
            else
               hz=zlabel('Counts','fontsize',20,'fontweight','bold');
            end
            set(hy,'Rotation',-35);
            set(hx,'Rotation',25);
            data.visual.var1_str=data.visual.var_str;
            data.visual.func1_str=data.visual.func_str;
            data.visual.var2=data.visual.var';
            if get(hdl.vis.radio_visual_func_length,'value') == 1               
               data.visual.var1=para.bio.hake_len_bin(:);
               data.visual.func1=data.visual.func(:,2:2:end)';               
            else
               data.visual.var1=para.bio.hake_age_bin(:);
               data.visual.func1=data.visual.func';
            end
            if isfield(data.visual,'var0_str')
                data.visual.var_str=data.visual.var0_str;
            end
            if isfield(data.visual,'func0_str')
                data.visual.func_str=data.visual.func0_str;
            end
            if  get(hdl.vis.radio_visual_var_len_age,'value') == 1
                data.visual.var1_str='Length (cm)';
                data.visual.var2_str='Age';
                data.visual.func1_str=data.visual.func0_str;
            end
        case 7  % 2D bubble plot
            lat=data.visual.var(:,1);
            lon=data.visual.var(:,2);
            if lon(1) > 0
              lon=-lon;
            end
            func=data.visual.func;    % no sex differentiation
            fac=1;   % display scale factor
            if get(hdl.vis.radio_visual_biological,'value') == 1 & (get(hdl.vis.radio_visual_func_length,'value') == 1 | get(hdl.vis.radio_visual_func_age,'value') == 1 )
                func.func=func;
                func.biomass=data.visual.func_biomass;
            end
            bubble_plot_NASC_number_biomass(lat,lon,func,data.visual.func_str,fac);           
            data.visual.var1=lat;
            data.visual.var2=lon;
            data.visual.var1_str='Latitude';
            data.visual.var2_str='Longitude';
            data.visual.func1=data.visual.func;
            data.visual.func1_str=data.visual.func_str;
        case 8    % 2D scatter plot
            figure(1)
            clf
            %% remove NaN's
            ind=find(~isnan(data.visual.var) == 1);
            data.visual.var=data.visual.var(ind);
            data.visual.func=data.visual.func(ind,:);
            if get(hdl.vis.radio_visual_func_age,'value') == 1 & get(hdl.vis.radio_visual_var_length,'value') == 1  % age vs length
                func=data.visual.func3;
                var=data.visual.var3;
                var_str=data.visual.var3_str;
                func_str=data.visual.func3_str;
                ind_i=find(~isnan(func) ==1);
                %% add 2*ns samples to force the curve go through reasonable (x,y) values
                ns=400;
                var_i=[4+0.5*randn(ns,1); 14+0.5*randn(ns,1); var(ind_i)];
                func_i=[zeros(ns,1); 0.5*ones(ns,1); func(ind_i)];
                p=polyfit(var_i,func_i,4);
                Lmin=min(var_i);
                Lmax=max(var_i);
                len=linspace(Lmin,Lmax,30);
                yfit=polyval(p,len);
                plot(var_i,func_i,'.',len,yfit,'-r','linewidth',2)
                legend('Data','Fitted Curve','location','northwest')
                h=text(0.0,0.75,sprintf('Age = %5.2eL^4+%4.2eL^3+%4.3fL^2+%4.2fL+%4.2f',p),'sc','fontweight','bold');
                axis([Lmin Lmax 0 para.bio.hake_age_bin(end)])
           else
                var=data.visual.var;
                func=data.visual.func;
                var_str=data.visual.var_str;
                func_str=data.visual.func_str;
                plot(var,func,'.')
            end
            xlabel(var_str,'fontsize',20,'fontweight','bold');
            ylabel(func_str,'fontsize',20,'fontweight','bold');
        case 9    % 2D boxplot
            var=data.visual.var3;
            func=data.visual.func3;
            var_str=data.visual.var_str;
            func_str=data.visual.func_str;
            func(func == 0)=nan;
            figure(14)
            dat=boxplot(func,'labels',var);
            h=draw_boxplot_xtick_labels(var);
            xlabel(var_str,'fontsize',20,'fontweight','bold');
            ylabel(func_str,'fontsize',20,'fontweight','bold');
            if get(hdl.vis.radio_visual_var_gender,'value') == 1
                set(gca,'xticklabels',{'Male','Female','ALL'})
            end
            n=length(find(~isnan(func) & func ~= 0));
            title(sprintf('Total non-zero samples = %d',n));
            if get(hdl.vis.radio_visual_func_length,'value') == 1 & get(hdl.vis.radio_visual_var_age,'value') == 1  % length vs age
                len=nanmean(func)';
                age=var;
                p=polyfit(age,len,4);
                yfit=polyval(p,para.bio.hake_age_bin); 
                figure(14)
                hold on
                plot(para.bio.hake_age_bin,yfit,'-g','linewidth',2)
                legend('Fitted Curve','location','northwest')
                h=text(0.2,0.2,sprintf('L = %5.4fa^4+%4.3fa^3+%4.2fa^2+%4.2fa+%4.2f (cm)',p),'sc');
                hold off
                data.visual.len_wgt.p=p;
                data.visual.len_wgt.yfit=yfit;
                data.visual.len_wgt.age_bin=para.bio.hake_age_bin;
            end
    end
end
function   h=draw_boxplot_xtick_labels(var_in)
set(gca,'xticklabel',{' '});
int_var=unique(round(var_in));
for i=1:length(int_var)
    indx=find(var_in == int_var(i));
    if isempty(indx)
        [val,indx]=min(abs(var_in-int_var(i)));
        ind(i)=indx;
        cell_labels{i}=num2str(var_in(indx));
    else
        ind(i)=indx;
       cell_labels{i}=num2str(int_var(i));
    end
end
set(gca,'xtick',ind,'xticklabel',cell_labels)
h=gca;
return
    



