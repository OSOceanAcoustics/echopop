function		[Vp,Ep,Eps]=krig2d(xp,yp,xn,yn,var,model_para,hmsg,xb,yb)
% function		[Vp,Ep, Eps]=krig2d(xp,yp,xn,yn,var,model_para,hmsg,xb,yb) perform the 2-D kriging
%% INPUT:  
%%       xn, yn - coordinates of the input data
%%       xp, yp - coordinates of the kriged output
%%          var - input data value array
%%   model_para - model parameters for computing the semi-variogram/correlogram 
%%         hmsg - handle of the message window
%%           xb - x-value (lon) array of west-bound transect boundary from SCB to NWVI
%%           yb - y-value (lat) array of west-bound transect boundary from SCB to NWVI
%%
%%  OUTPUT:
%%            Vp - kriged values
%%            EP - Kriging variance
%%            Eps = sample variance
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.
%%  Last modification  June 24, 2013 

global para data 

model=para.vario.model;
model_type=2-para.vario.corr;
mod_opt=para.krig.model;
eplim=para.krig.elim;
EPS=para.krig.eps;
Ratio=para.krig.ratio;
blk_opt=para.krig.scheme;
nx=length(xp);
ny=length(yp);
nbx=para.krig.blk_nx;
nby=para.krig.blk_ny;
kmin=para.krig.kmin;
kmax=para.krig.kmax;
R=para.krig.srad;
Rtaper=para.vario.lscl;
% R=3*Rtaper;
% kmin=1;
% kmax=100;

dx=para.krig.dx;
dy=para.krig.dy;


%% process anisotrophic data
if para.proc.default_para == 1   
        cz_in_obs=[xn(:) yn(:)];
        cz_out_obs=[xn(:) yn(:)/para.vario.ytox_ratio];
        xn=cz_out_obs(:,1);yn=cz_out_obs(:,2);
        %% predictions
        cz_in_pred=[xp(:) yp(:)];
        cz_out_pred=[xp(:) yp(:)/para.vario.ytox_ratio];
        xp=cz_out_pred(:,1);yp=cz_out_pred(:,2);    
else
    if abs(para.vario.azm_rot) > eps |  abs(para.vario.dip_rot) > eps ...
            | abs(1-para.vario.ytox_ratio) > 1e-5 | abs(1-para.vario.ztox_ratio)
        %% observations
        cz_in_obs=[xn(:) yn(:)];
        cz_out_obs=coordtransform3d(cz_in_obs);
        xn=cz_out_obs(:,1);yn=cz_out_obs(:,2);
        %% predictions
        cz_in_pred=[xp(:) yp(:)];
        cz_out_pred=coordtransform3d(cz_in_pred);
        xp=cz_out_pred(:,1);yp=cz_out_pred(:,2);
    end
end

n=length(xn);
xn=reshape(xn,1,n);
yn=reshape(yn,1,n);

if model_type == 1
  model=-abs(model);
end
nmax_var=500;
rR=linspace(0,sqrt(2)+0.01,nmax_var);dr=rR(2)-rR(1);
vgram=variogrammodel3d(model,rR,model_para);
if blk_opt == 2 
   model_para1=[0 1 model_para(3:5)];
% The semi-variogram/covariance between a grid block and itself:
   CVV=blk2blk2d(0,0,0,0,nbx,nby,dx,dy,model,model_para1); 
end
pad=zeros(1,4);		% needed to pad K in the main routine

if 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if blk_opt == 2 & model < 0
        model_para1=[0 1 model_para(3:5)];
        % The covariance between a grid block and itself:
        CVV=blk2blk2d(0,0,0,0,nbx,nby,dx,dy,model,model_para1);
    end
    
end

n=length(xn);
val_m=max(var);
Anomaly_cnt1=0;
Anomaly_cnt2=0;
Anomaly_cnt3=0;
nv=nx;
nv_cnt=max(1,round(0.01*nv));
Vp=nan*ones(nv,1);
Ep=nan*ones(nv,1);
Eps=nan*ones(nv,1);

if mod_opt == 3 			% for universal kriging with linear drift 
%  delx=-(xp(ones(n,1),:)'-xn(ones(nv,1),:));
%  dely=-(yp(ones(n,1),:)'-yn(ones(nv,1),:));
end
pad=zeros(1,3);		% needed to pad K in the main routine

if nargin >= 7 tic;end
kkk=0;

for i=1:nv    % loop through grid cell
    %     if xp(i) < -0.4 & yp(i) >= 0.2 & yp(i)  <=0.4
    %         disp(i)
    %         disp([xp(i) yp(i)])
    %     end
%     if i == 5900
%         disp(i)
%     end
    Dx=(xp(i)-xn);
    Dx2=Dx.*Dx;
    Dy=(yp(i)-yn);
    Dy2=Dy.*Dy;
    if rem(i,nv_cnt) == 0 & nv > 100 & nargin >= 7
        if para.proc.bootstrap.limit > 1
            percent=sprintf('Bootstraping: %d of %d:   %g percent done in %g sec.', para.proc.bootstrap.cnt, para.proc.bootstrap.limit, floor(100*i/nv),round(toc));
            if ~isempty(hmsg) krig_message(4,percent,' ',hmsg); end
        else
            percent=sprintf('   %g percent done in %g sec.', floor(100*i/nv),round(toc));
            if ~isempty(hmsg) krig_message(4,percent,' ',hmsg); end
        end
    end
    r=sqrt(Dx2+Dy2);
    M2_unity=1;
    if kmin >= n
        indx_sort=1:n;
        indx1=1:n;nk=n;
        fprintf('** kmin > n *** \n');
    else
        [r_sort,indx_sort]=sort(r);
        indx=find(r_sort <= R);
        if length(indx) > kmax
            indx1=indx_sort(1:kmax);
        elseif length(indx) < kmin
            indx1=indx_sort(1:kmin);
    %%% provide a tapered function to handle extrapolation
%             x_sort=min(xn(indx_sort(indx1)));
%             y_sort=mean(yn(indx_sort(indx1)));     
            [var_dif,sort_ind]=min(abs(yp(i)-yb));
            if xp(i) < xb(sort_ind) - R   % kriging grid x-axis < west_bound_x - R search radius  
%                M2_unity(i)=exp(-abs(nanmean((r_sort(indx1)-Rtaper)))/Rtaper);
%                M2_unity(i)=exp(-max(abs(r(indx1)-Rtaper))/Rtaper);
               M2_unity=exp(-nanmean(r_sort(indx1))/R);
            end
        else
            indx1=indx_sort(indx);
        end
        nk=length(indx1);
%         smpl_n(i)=nk;
%         for j=1:30
%         var_ave(i,j)=mean(var(indx_sort(1:j)));
%         end
    end

    if blk_opt == 2
        M20=sta2blk2d(xn(indx1),yn(indx1),xp(i),yp(i),nbx,nby,dx,dy,model,model_para);
    else
        M20=variogrammodel3d(model,r(indx1),model_para);
    end
    switch mod_opt
        case 1		% Objective mapping, follow Journel and Huijbregts, p. 307
            M2=M20';
        case 2		% Ordinary Kriging, follow Journel and Huijbregts, p. 307
            M2=[M20 1 ]';
        case 3		% Universal Kriging w/ Linear drift, follow Journel and Huijbregts, p. 319
            M2=[M20 1 xp(i) yp(i) ]';		% x,y
    end
    if i == 1 | kmin < n				% for kmin > n only needs to compute SVD once
        % Construct the coefficient matrix and its inverse matrix
        % to avoid repeated computation in the loop
        x1=xn(indx1);
        y1=yn(indx1);
        var1=var(indx1);
        ind0=find(var1 > 1e-5);
        nvar_pts(i)=length(ind0);
        npts(i)=n;
        dx1=(x1(ones(nk,1),:)'-x1(ones(nk,1),:));
        dy1=(y1(ones(nk,1),:)'-y1(ones(nk,1),:));
        rs=sqrt(dx1.*dx1+dy1.*dy1);
        K0=variogrammodel3d(model,rs,model_para);	% semi-variagram
        if model_type == 1								% correlogram
            K0(1:nk+1:nk*nk)=ones(nk,1);
        else													% variogram
            K0(1:nk+1:nk*nk)=zeros(nk,1);
        end
        switch mod_opt
            case 1		% Objective mapping, follow Journel and Huijbregts, p. 307
                K=K0;
            case 2		% Ordinary Kriging, follow Journel and Huijbregts, p. 307
                K=[K0 ones(nk,1);ones(1,nk) 0];
            case 3		% Universal kriging w/ linear drift, follow Journel and Huijbregts, p. 319
                K=[K0 ones(nk,1) xn(indx1)' yn(indx1)' ; ...
                    ones(1,nk) pad; xn(indx1) pad; yn(indx1) pad];
        end
        [U,S,V]=svd(K);
        %   kindx=find(diag(S) > EPS);
        kindx=find(abs(diag(S)/S(1,1)) > Ratio);
        Sinv=diag(1./diag(S(kindx,kindx)));
        K_inv=V(:,kindx)*Sinv*U(:,kindx)';
    end																		% end of if loop for computing SVD
    lambda=K_inv*M2;
    lnan=find(isnan(var1));
    if length(lnan) > 0
        var1(lnan)=zeros(length(lnan),1);
    else
        var1(find(r(indx1) > R))=0;
    end
    Vp(i)=sum_nan(lambda(1:nk).*var1)*M2_unity;
    Vp_taper(i)=sum_nan(lambda(1:nk).*var1);
    Vpm(i)=nanmean(var1);
    ind_non0=find(abs(var1) > 2*eps);
    ind_non0=1:nk;
    CVs(i)=sqrt(nanvar(var1(ind_non0)))/abs(nanmean(var1(ind_non0)));              % localized sample CV within the search redius
    if abs(mean_nan(var(indx1)-Vp(i))) > 2
        % 		  disp(sprintf('|Vp-mean(var)|>2 i=%g',i));
        %       plot(1:length(indx1),var(indx1),'.-',length(indx1)/2,Vp(i),'o')
        Anomaly_cnt1=Anomaly_cnt1+1;
    end
    if blk_opt == 2
        if model_type == 1
            Ep(i)=CVV-sum_nan(lambda.*M2);	% Error estimate for tp(j)
        else
            Ep(i)=sum_nan(lambda.*M2)-CVV;	% CVV = GammaVV
        end
    else
        if model_type == 1
            Ep(i)=1-sum_nan(lambda.*M2);
        else
            Ep(i)=sum_nan(lambda.*M2);
        end
    end
    if abs(Vp(i)) > 3.0*val_m
        %			disp(sprintf('Var > 3*Var_max i=%g',i));
        Anomaly_cnt2=Anomaly_cnt2+1;
        %    		Ep(i)=3.0*val_m;
    end
    Eps(i)=sqrt(Ep(i)*nanvar(var1(ind_non0)))/abs(Vp(i));  % sample CV
    if Ep(i) < 0
        kkk=kkk+1;
        fprintf('%d\t Ep = %8.4g\t M2_unity = %5.3f\t  nt = %d\n',i,Ep(i),M2_unity,kkk);
        Ep(i)=abs(Ep(i));
    end
% if  Eps(i)> 10
% disp(num2str([i  Vp(i) nanstd(Ep(i)) Eps(i) ]))
% end
%     Vp_s(i)=nanmean(var1);
%     C0=9.0938e+09;    % variance
%     Var_k(i)=Ep(i)*C0;
%     Var_k0(i)=(var1-Vp(i))'*((var1-Vp(i)))/(length(var1)-1);    % kriging variance
%     Var_s(i)=nanstd(var1).^2;                                     % sample variance
%     Epp(i)=std(var1)^2;
    if abs(Ep(i)) > eplim
        %			disp(sprintf('Err > eplim i=%g',i));
        Vp(i)=NaN;
        Ep(i)=abs(Ep(i));
        Anomaly_cnt3=Anomaly_cnt3+1;
    elseif abs(Ep(i)) <= eps & nk > kmin
        fprintf('Possible bad kriging grid cell at: %d\t  Ep = %4.2f\n',i,Ep(i));
    end
end

if  Anomaly_cnt1 > 1
    %  disp(sprintf(' Anomaly_cnt for |<var(indx1)-Vp(i)>| > 2 = %g',Anomaly_cnt1));
end
if  Anomaly_cnt2 > 1
    %  disp(sprintf(' Anomaly_cnt for Vp > 3*Var_max =%g ',Anomaly_cnt2));
end
if  Anomaly_cnt3 > 1
    disp(sprintf(' Anomaly_cnt for |Ep| > Relative Error =%g ',Anomaly_cnt3));
end
% ind=find(M2_unity < 1);
% tot_var_taper=nansum(Vp_taper(ind));
% tot_var=nansum(Vp_taper);
% tot_Vp=nansum(Vp);
% fprintf('\n tot_var = %9.3e, \t tot_Vp = %9.3e\n',tot_var,tot_Vp);
% fprintf('No_taper = %d, \t tot_var_taper = %9.3e, \t Vp = %9.3e\n',length(ind),tot_var_taper,tot_var-tot_Vp);
% fprintf('--- Ratio = %10.7f --- \n\n',tot_var_taper/tot_Vp);
% figure(6); plot(xp,yp,'.',xp(ind),yp(ind),'.r')
return
