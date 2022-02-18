function     [np, dis, gam, hm, tm, hv, tv, xx, yy, zz] = variogram3d(nd, x, y, z, vr, nlag, xlag, xltol, ...
          ndir, azm, atol, bandwh, dip, dtol, bandwd, nvarg, ...
          ivtail, ivhead, ivtype, nvar, isill);
% function variogram3d() computes normalized 2d/3d semi-variogram/correlogram model parameters
% 
%% INPUT
%    data.in.tv=f(x,y)		 (x,y) coordinates of transformed variable value
%    para.vario.res   = lag resolution
%    para.vario.range = range limit for semivariogram computation
%    para.vario
%            .angle = orientation angle from horizontal axis in degree
%            .ratio = aspect ratio of the scale in vertical (Y) to that of horizontal (X)
%% OUTPUT
%    data.out.vario.grd2d =  angle-lag grid
%    data.out.vario.gammah2d = gamma(h) 	or gamma(grd)   semi-variogram or correlogram
%
%%
%%  Kriging Software Package  version 3.0,   December 29, 2001
%%  Copyright (c) 1998, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para data


if nd < 1
    krig_message(1,'No usable data !!');
    return
end

%% parameters
res=xlag;
%range=sqrt(2);

azm_res=2*mean(atol);
dip_res=2*mean(dtol);
lag=0:res:res*(nlag-1);

msg1='Computing Semi-variogram/Correlogram, please be patient,... ';
msg2=sprintf('usable data points for each angle = %g',nd);
hdl0=krig_message(2,msg1,msg2);
tic

nnmax=ndir*nlag;
np=zeros(nnmax,1);
dis=nan*ones(nnmax,1);
gam=nan*ones(nnmax,1);
hm=nan*ones(nnmax,1);
tm=nan*ones(nnmax,1);
hv=nan*ones(nnmax,1);
tv=nan*ones(nnmax,1);

if nd*ndir < 1000
    nl_cnt=max(1,round(0.05*nd*ndir));
    dnorm=100;
elseif nd*ndir < 10000
    nl_cnt=max(1,round(0.01*nd*ndir));
    dnorm=100;
elseif nd*ndir < 100000
    nl_cnt=max(1,round(0.002*nd*ndir));
    dnorm=500;
else
    nl_cnt=max(1,round(0.001*nd*ndir));
    dnorm=1000;
end

%% update data based on sepecified anisotropy parameters
ii=sqrt(-1);

data.out.vario.len=nan*ones(ndir,nd);

for j=1:ndir						% orientation loop
    angz=azm(j);
    angd=dip(j);
    %% 2d-coordinates transformation
    % initialization
    cnt=zeros(nlag-1,1);
    gamma_sum=zeros(nlag-1,1);
    var_sum=zeros(nlag-1,1);
    var_mean=zeros(nlag-1,1);
    var2_sum=zeros(nlag-1,1);
    gammah=nan*zeros(nlag-1,1);
    c0=zeros(nlag-1,1);
    head_indx=zeros(nd,nlag-1);
    tail_indx=zeros(nd,nlag-1);
    ss=[];
    for i=1:nd					% sample loop
        ij=(j-1)*nd+i;
        if rem(ij,nl_cnt) == 0 & ~isempty(hdl0)
            percent=sprintf('%g percent done in %g sec.', floor(dnorm*ij/(ndir*nd))*100/dnorm,round(toc));
            krig_message(4,percent,' ',hdl0);
        end
        dx=x(i+1:nd)-x(i);
        dy=y(i+1:nd)-y(i);
        dz=z(i+1:nd)-z(i);
        dr=sqrt(dx.*dx+dy.*dy);
        %  dr3=sqrt(dx.*dx+dy.*dy+dz.*dz);
        ang_z=atan(dy./dx)*180/pi;           % change atan2(y,x) --> atan(y./x) on Nov 14, 2006
        ang_d=atan(dz./dr)*180/pi;           % change atan2(y,x) --> atan(y./x) on Nov 14, 2006
        neg_indx=find(ang_z < angz - 0.5*azm_res);
        if ~isempty(neg_indx)
            ang_z(neg_indx)=180+ang_z(neg_indx);
        end
        ang_indx=find(ang_z >= angz - 0.5*azm_res & ang_z < angz+0.5*azm_res  ...
            & ang_d >= angd - 0.5*dip_res & ang_d < angd + 0.5*dip_res);
        %  if j == 11
        %      disp([i x(i) y(i) length(ang_indx)])
        %  end
        if ~isempty(ang_indx)
            if nd-i ~= length(ang_indx)
                %       disp([i nd-i length(ang_indx) min(ang_z) max(ang_z)])
            end
            data.out.vario.len(j,i)=length(ang_indx);
            xi=x(ang_indx);yi=y(ang_indx);zi=z(ang_indx);var_i=vr(ang_indx);
            %plot(x(i+1:n),y(i+1:n),'.k',x(ang_indx),y(ang_indx),'or');pause
            ni=length(xi);
            dxi=x(i)-xi;			% dxi=dx(ang_indx-1)
            dyi=y(i)-yi;
            dzi=z(i)-zi;
            d=sqrt(dxi.*dxi+dyi.*dyi+dzi.*dzi);
            [dsort,indx_sort]=sort(d);
            var_sort=var_i(indx_sort);
            dindx=round(dsort/res)+1;
            k_acc=0;								% accumulated k index
            for k=1:nlag-1
                indxk=find(dindx(k_acc+1:ni) == k) + k_acc ;
                nk=length(indxk);
                if nk > 0
                    % tail portion
                    cnt(k)=cnt(k)+nk;
                    k_acc=k_acc+nk;
                    var_dif= vr(i) - var_sort(indxk);
                    var_sum(k)=var_sum(k)+sum(var_sort(indxk));
                    var2_sum(k)=var2_sum(k)+var_sort(indxk)'*var_sort(indxk);
                    gamma_sum(k)=gamma_sum(k)+ var_dif'*var_dif;
                    % head point index
                    head_indx(i,k)=nk;
                    jindx=indx_sort(indxk)+i;
                    tail_indx(jindx,k)=tail_indx(jindx,k)+ones(nk,1);
                    % 		      if k == 2 & i == 70000
                    % 					for kk=1:nk
                    %   						jindx=indx_sort(indxk(kk))+i;
                    %                	disp([i jindx cnt(k) var_sort(indxk(kk)) vr(i) dsort(indxk(kk)) gamma_sum(k)])
                    %                     end
                    %  			  end
                end
            end					% end of lag loop
        end					% not empty
    end					% end of selected observation data loop
    sigma_head=zeros(nlag-1,1);
    mean_head=zeros(nlag-1,1);
    sigma_tail=zeros(nlag-1,1);
    mean_tail=zeros(nlag-1,1);
    for k=1:nlag-1
        indx1=find(head_indx(:,k) >= 1);
        hl=length(indx1);
        if hl > 0
            npk=sum_nan(head_indx(:,k));
            pdf=head_indx(indx1,k)/npk;
            mean_head(k)=sum_nan(vr(indx1).*pdf);
            sigma_head(k)=sqrt(sum_nan(((vr(indx1)-mean_head(k)).^2).*pdf));
        end
        indx1=find(tail_indx(:,k) >= 1);
        hl=length(indx1);
        if hl > 0
            npk=sum_nan(tail_indx(:,k));
            pdf=tail_indx(indx1,k)/npk;
            mean_tail0(k)=sum_nan(vr(indx1).*pdf);
            sigma_tail0(k)=sqrt(sum_nan(((vr(indx1)-mean_tail0(k)).^2).*pdf));
        end
    end
    
    indx=find( cnt >= 1);
    if length(indx) > 0
        cnt_nz=cnt(indx);
        var_mean=var_sum(indx)./cnt_nz;
        var_std2=var2_sum(indx)./cnt_nz-var_mean.^2;
        mean_tail(indx)=var_mean;
        sigma_tail(indx)=sqrt(abs(var_std2));
        c0(indx)=sigma_tail(indx).*sigma_head(indx);
        % c0(indx)=0.5*(mean_tail(indx)+mean_head(indx));
        gammah(indx)=0.5*gamma_sum(indx)./(cnt_nz.*mean(c0(indx))+eps);
        %  c0=std(vr)^2;                    % normalization factor   % used in  original version and commented out on 11-06-2006
        gammah(indx)=0.5*gamma_sum(indx)./(cnt_nz.*c0(indx)+eps);  % used in the revised on the Nov. 06, 2006
    else
        gammah=nan.*ones(nlag-1,1);
    end
    indx=find(gammah > max(2.5,para.vario.max_value));
    gammah(indx)=nan*ones(size(indx));
    indx_ij=(j-1)*nlag+1:j*nlag;
    gam(indx_ij)=[0;gammah]';
    dis(indx_ij)=lag;
    np(indx_ij)=[nd-1;cnt];
    hm(indx_ij)=[0;mean_head];
    tm(indx_ij)=[0;mean_tail];
    hv(indx_ij)=[0;sigma_head];
    tv(indx_ij)=[0;sigma_tail];
    xx(indx_ij)  = lag*cos(angz*pi/180)*cos(angd*pi/180);
    yy(indx_ij)  = lag*sin(angz*pi/180)*cos(angd*pi/180);
    zz(indx_ij)  = lag*sin(angd*pi/180);
end									% angle of orientation loop
close(hdl0);