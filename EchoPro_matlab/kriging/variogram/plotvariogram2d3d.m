function 	plotvariogram2d3d(opt)
% function 	plotvariogram2d3d(opt) displays the 2D/3D semi-variogram/correlogram
% opt = 1	:   initial or re-plot with a new set of variogram/correlogram
%     = 2   :   replot selected variogram in azimuth (3D)
%     = 3   :   replot selected variogram in elevation  (3D)
%     = 4   :   adjust color scale display
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para hdl data color


if opt == 1 
  plotvariogram2d3dfig;
end
   
switch opt		
    case 1
        if para.vario.dim == 2
            n=data.out.vario.nlag2d3d;								% number of lags at each direction
            m=max(para.vario.nazm,para.vario.ndip);		% number of directions
            vario_corr=reshape(data.out.vario.gammah2d3d,n,m);
            x=reshape(data.out.vario.xx,n,m);
            y=reshape(data.out.vario.yy,n,m);
        else
            n=data.out.vario.nlag2d3d;								% number of lags at each direction
            m=para.vario.nazm;										% number of directions in azimuth
            l=para.vario.ndip;
            nidip=round(l/2);
            indx=n*m*(nidip-1)+1:n*m*nidip;
            vario_corr=reshape(data.out.vario.gammah2d3d(indx),n,m);
            x=reshape(data.out.vario.xx(indx),n,m);
            y=reshape(data.out.vario.yy(indx),n,m);
            set(hdl.dispvario2d3d.azm_val,'string','0-360 deg');
            mean_dip_ang=(para.vario.dip_beg+para.vario.dip_end)/2;
            set(hdl.dispvario2d3d.dip_slider,'value',0.5)
            set(hdl.dispvario2d3d.azm_slider,'value',0)
            set(hdl.dispvario2d3d.dip_val,'string',[num2str(mean_dip_ang)  '  deg'])
            set(hdl.dispvario2d3d.azm_val,'enable','off')
            set(hdl.dispvario2d3d.azm_slider,'enable','on')
        end
        xmax=max(max(x));
        ymax=max(max(y));
        x=[x x(:,end)];
        y=[y zeros(size(x,1),1)];
        vario_corr=[vario_corr vario_corr(:,end) ];
        xx=x;
        yy=-y;
        if get(hdl.vario.variogram,'Value') == 1
            hdl.dispvario2d3d.axes1=pcolor(y,x,vario_corr);
            hold on
            hdl.dispvario2d3d.axes1=pcolor(yy,xx,vario_corr);
            hold off
        else
            hdl.dispvario2d3d.axes1=pcolor(y,x,1-vario_corr);
            hold on
            hdl.dispvario2d3d.axes1=pcolor(yy,xx,1-vario_corr);
            hold off
        end
        if get(hdl.dispvario2d3d.shading_radio1,'value') == 1
            shading faceted;
        elseif get(hdl.dispvario2d3d.shading_radio2,'value') == 1
            shading flat;
        elseif get(hdl.dispvario2d3d.shading_radio3,'value') == 1
            shading interp;
        end
        hdl.dispvario2d3d.htxt2d(1)=text(xmax*1.15,0,'90^o');
        hdl.dispvario2d3d.htxt2d(2)=text(0,ymax*1.15,'0^o');
        hdl.dispvario2d3d.htxt2d(4)=text(0,-ymax*1.15,'180^o');
        hdl.dispvario2d3d.htxt2d(3)=text(-xmax*1.4,0,'270^o');
        set(hdl.dispvario2d3d.htxt2d(1:4),'fontweight','bold','fontsize',12);
        axis(1.1*[-xmax xmax -ymax ymax])
    case 2								% adjust azm angle
        n=data.out.vario.nlag2d3d;								% number of lags at each direction
        m=para.vario.nazm;										% number of directions in azimuth
        l=para.vario.ndip;										% number of directions in dip
        rate=(para.vario.azm_end-para.vario.azm_beg);
        azm_ang=get(hdl.dispvario2d3d.azm_slider,'value')*rate+para.vario.azm_beg;
        niazm=min(floor(round(get(hdl.dispvario2d3d.azm_slider,'value')*m)+1),m);
        set(hdl.dispvario2d3d.azm_val,'string',num2str(azm_ang))
        indx=(niazm-1)*n+1:n*niazm;
        for i=2:l
            indx1=n*m*(i-1)+(niazm-1)*n+1;
            indx=[indx indx1:indx1+n-1];
        end
        vario_corr=reshape(data.out.vario.gammah2d3d(indx),n,l);
        range=linspace(0,para.vario.range,n);
        dip_ang=linspace(para.vario.dip_beg,para.vario.dip_end,l);
        x=range(:)*cos(dip_ang*pi/180);
        z=-range(:)*sin(dip_ang*pi/180);
        xmax=max(max(abs(x)));
        zmax=max(max(abs(z)));
        if get(hdl.vario.variogram,'Value') == 1
            hdl.dispvario2d3d.axes1=pcolor(x,z,vario_corr);
        else
            hdl.dispvario2d3d.axes1=pcolor(x,z,1-vario_corr);
        end
        if get(hdl.dispvario2d3d.shading_radio1,'value') == 1
            shading faceted;
        elseif get(hdl.dispvario2d3d.shading_radio2,'value') == 1
            shading flat;
        elseif get(hdl.dispvario2d3d.shading_radio3,'value') == 1
            shading interp;
        end
        hdl.dispvario2d3d.htxt2d(1)=text(xmax*1.15,0,'0^o');
        hdl.dispvario2d3d.htxt2d(2)=text(0,zmax*1.15,'-90^o');
        hdl.dispvario2d3d.htxt2d(3)=text(0,-zmax*1.15,'90^o');
        set(hdl.dispvario2d3d.htxt2d(1:3),'fontweight','bold','fontsize',12);
        axis(1.1*[-xmax xmax -zmax zmax])
    case 3															% adjust dip angle
        n=data.out.vario.nlag2d3d;								% number of lags at each direction
        m=para.vario.nazm;										% number of directions in azimuth
        l=para.vario.ndip;										% number of directions in dip
        rate=(para.vario.dip_end-para.vario.dip_beg);
        dip_ang=get(hdl.dispvario2d3d.dip_slider,'value')*rate+para.vario.dip_beg;
        nidip=min(floor(round(get(hdl.dispvario2d3d.dip_slider,'value')*l)+1),l);
        set(hdl.dispvario2d3d.dip_val,'string',num2str(dip_ang))
        indx=n*m*(nidip-1)+1:n*m*nidip;
        vario_corr=reshape(data.out.vario.gammah2d3d(indx),n,m);
        range=linspace(0,para.vario.range,n);
        azm_ang=linspace(para.vario.azm_beg,para.vario.azm_end,m);
        x=range(:)*cos(azm_ang*pi/180);
        y=range(:)*sin(azm_ang*pi/180);
        xx=x;
        yy=-y;
        xmax=max(max(x));
        ymax=max(max(y));
        if get(hdl.vario.variogram,'Value') == 1
            hdl.dispvario2d3d.axes1=pcolor(y,x,vario_corr);
            hold on
            hdl.dispvario2d3d.axes1=pcolor(yy,xx,vario_corr);
            hold off
        else
            hdl.dispvario2d3d.axes1=pcolor(y,x,1-vario_corr);
            hold on
            hdl.dispvario2d3d.axes1=pcolor(yy,xx,1-vario_corr);
            hold off
        end
        if get(hdl.dispvario2d3d.shading_radio1,'value') == 1
            shading faceted;
        elseif get(hdl.dispvario2d3d.shading_radio2,'value') == 1
            shading flat;
        elseif get(hdl.dispvario2d3d.shading_radio3,'value') == 1
            shading interp;
        end
        hdl.dispvario2d3d.htxt2d(2)=text(0,ymax*1.15,'0^o');
        hdl.dispvario2d3d.htxt2d(1)=text(xmax*1.15,0,'90^o');
        hdl.dispvario2d3d.htxt2d(3)=text(-xmax*1.4,0,'270^o');
        hdl.dispvario2d3d.htxt2d(4)=text(0,-ymax*1.15,'180^o');
        set(hdl.dispvario2d3d.htxt2d,'fontweight','bold','fontsize',12);
        axis(1.1*[-xmax xmax -ymax ymax])
    case 4
        DispVar=max(get(hdl.dispvario2d3d.axes1,'Ydata'));
        var1=get(hdl.dispvario2d3d.cbar_slider_bot,'value');
        var2=get(hdl.dispvario2d3d.cbar_slider_top,'value');
        max_var=max(max(max(DispVar)));
        min_var=min(min(min(DispVar)));
        step=(max_var-min_var)/100;
        disp_var1=min_var+(max_var-min_var)*var1;
        disp_var2=max_var-(max_var-min_var)*(1-var2);
        caxis([disp_var1 disp_var2]);
end

axis off
axis square
hdl.dispvario2d3d.colorbar2=colorbar;
set(hdl.dispvario2d3d.colorbar2,'Position',[0.85 0.4 0.045 0.5]);
%set(hdl.dispvario2d3d.colorbar2,'visible','off');
%set(get(hdl.dispvario2d3d.colorbar2,'children'),'visible','off');
