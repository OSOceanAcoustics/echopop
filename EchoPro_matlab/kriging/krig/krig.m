function [tp,ep]=krig(mod_opt,xc,yc,x,y,var,krig_para,model,model_para,hmsg)
% function [tp,ep]=krig(mod_opt,xc,yc,x,y,var,krig_para,model,model_para,hmsg)
% performs kriging and returns the kriged results at locations (xc,yc) using observed
% data at var(x,y). The kriging resolution is specified by res.
% Input Parameters:
% mod_opt:  determine what kind of kriging model is selected:   
%       mod_opt=1:     object analysis
%       mod_opt=2:     ordinary kriging -- second order stationary case
%       mod_opt=3;     kriging  -- universal kriging with linear drift w/o revision
% xc, yc:    coordinates vectors for specified kriging points  
%            if the kriging grid is N x M, the their length is N*M.
%    x,y:    vectors with the coordinates of the observation points
%    var:    vector with the measurements of the variable to be mapped
%
%% kriging parameters
% krig_para = [blk_opt,10,range,res,kmin,kmax,eplim,rot_ang,aspect];
% blk_opt:  option for whether using station to station kriging ...
%            or station - block kriging     
%       blk_opt=1:     station to station (excluding the observation itself)
%       blk_opt=2:     block to station
%    nbx:    no. of horizontal elements in a block,              
%    nby:    no. of vertical elements in a block, i.e, block size = nbx * nby             
%  range:    search radius for kriging
%   kmin:    minimum observed data points used for kriging, > 3
%   kmax:    maximum observed data points used for kriging
%            setting a reasonable kmax can speed up computations
%  eplim:    normalized kriging error limit
% anisotropy_para =[rot_ang aspect]
%  rot_ang:  rotate angle for anisotropy
%  aspect:   take into account anisotropy  scale(y)/scale(x)
%  model:    index of variogram or covariance model
%  model_para: variogram or correlogram model parameters
%  hmsg:		 message window handle
%
% Output parameters:
%
% tp:    the estimate of the measured variable at the grid-points
% ep:    the estimation variance for tp, divided by the total variance
% 
% By Dezhang Chu, Woods Hole Oceanographic Institution
% Last modification  31 August, 1999
% %
% Reference:
%  A. J. Journel and Ch. J. Huijbregts, Mining Geostatistics
%  Academic Press, New York, pp. 304-313, 1978. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.
%%
%%		modified on 2/4/00	to handle data files without enough data points

dbstop if error

if nargin < 10
  hmsg=[];
end
blk_opt=krig_para(1);
nbx=krig_para(2);
nby=krig_para(3);
range=krig_para(4);
kmin=krig_para(5);
kmax=krig_para(6);
eplim=krig_para(7);
dx=krig_para(8);
dy=krig_para(9);
EPS=krig_para(10);
if ~isempty(hmsg) tic;end

n=length(x);		% n = number of data points
ngx=length(xc);
ngy=length(yc);
if n < mod_opt+1	% added on 2/4/00
  tp=NaN*ones(ngx,1);
  ep=NaN*ones(ngx,1);
  return				
end
if ngx > 1e6			% 32 MB
   krig_message(1,sprintf('Total Grid points nx*ny  = %g', ngx));
end
   
kn=length(krig_para);
if  kn > 10
  anisotropy_para=krig_para(kn-1:kn);
% dealing with anisotropy
  ang=anisotropy_para(1);
  aspect=anisotropy_para(2);
  if abs(ang) > eps | abs(1-aspect) > 1e-5
%% observations
    cz=[x(:) y(:)];
    rx=1;
    ry=aspect;
    [cz,rot]=trans2d(cz,ang,rx,ry);
    x=cz(:,1);y=cz(:,2);
%% predictions
    cz=[xc(:) yc(:)];
    [cz,rot]=trans2d(cz,ang,rx,ry);
    xc=cz(:,1);yc=cz(:,2);
  end
end

if kmin < mod_opt 
   kmin=mod_opt;
end
   
x=reshape(x,1,n); 	% makes sure that x is row vector
y=reshape(y,1,n);
nv=length(xc);				% nv = number of grid points
xc=reshape(xc,1,nv);
yc=reshape(yc,1,nv);
   
var=reshape(var,n,1); 	% and var is column vector

if model > 0				% variogram
   model_type=2;
else
   model_type=1;			% covariance
end

if mod_opt == 3 			% for universal kriging with linear drift 
  delx=-(xc(ones(n,1),:)'-x(ones(nv,1),:));
  dely=-(yc(ones(n,1),:)'-y(ones(nv,1),:));
end
pad=zeros(1,3);		% needed to pad K in the main routine


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if blk_opt == 2 & model < 0
   model_para1=[0 1 model_para(3:5)];
% The covariance between a grid block and itself:
   CVV=blk2blk(0,0,0,0,nbx,nby,dx,dy,model,model_para1); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NaNcnt=0;
Anomaly_cnt=0;
ep=zeros(nv,1);
tp=zeros(nv,1);
val_m=max(abs(var));

% Ordinary Kriging routine
nv_cnt=max(1,round(0.01*nv));
for j=1:nv				% kriging through gridded points
  if rem(j,nv_cnt) == 0 
     percent=sprintf('   %g percent done in %g sec.', floor(100*j/nv),round(toc));
     if ~isempty(hmsg) krig_message(4,percent,' ',hmsg); end
  end
  difx=xc(j)-x;
  dify=yc(j)-y;
  r=sqrt( difx.*difx+dify.*dify);  % grid to station
  if kmin >= n
     indx_sort=1:n;
     k=1:n;nk=n;
  else
     [r_sort, indx_sort]=sort(r);
     ind=find(r_sort <= range);
     nd=length(ind);
     if isempty(ind)
        k=1:min(kmin,n);
     elseif ind(nd) >= kmin & ind(nd) < kmax 
        k = 1:ind(nd);
     elseif ind(nd) >= kmax
        k = 1:kmax;
     elseif ind(nd) < kmin
        k = 1:min(kmin,n);
     end
     nk=length(k);
  end
  %% station to grid 
  if blk_opt == 2 & model < 0        
     M20=sta2blk(x(indx_sort(k)),y(indx_sort(k)),xc(j),yc(j),nbx,nby,dx,dy,model,model_para);
  else
     M20=variogrammodel(model,r(indx_sort(k)),model_para);
  end
% plot(x(indx_sort(k)),y(indx_sort(k)),'.',xc(j),yc(j),'+r');pause
 %%%% model selection %%%%%
  switch mod_opt
      case 1		% Objective mapping, follow Journel and Huijbregts, p. 307
        M2=M20';
      case 2		% Ordinary Kriging, follow Journel and Huijbregts, p. 307
        M2=[M20 1 ]';
      case 3		% Universal Kriging w/ Linear drift, follow Journel and Huijbregts, p. 319
        M2=[M20 1 0 0]';
  end
   
% Construct the coefficient matrix and its inverse matrix
% to avoid repeated computation in the loop
  if  j == 1 | kmin < n 
    rs=sqrt(((x(ones(length(k),1),indx_sort(k))' ...
          -x(ones(length(k),1),indx_sort(k))).^2+...
          (y(ones(length(k),1),indx_sort(k))' ...
          -y(ones(length(k),1),indx_sort(k))).^2));
    K0=variogrammodel(model,rs,model_para);
    if model_type == 1
       K0(1:nk+1:nk^2)=ones(nk,1);	   		% diagonal = 1
    else
       K0(1:nk+1:nk^2)=zeros(nk,1);	   		% diagonal = 0
    end
    K=[K0 ones(nk,1);ones(1,nk) 0];

      switch mod_opt
      case 1		% Objective mapping, follow Journel and Huijbregts, p. 307
        K=K0;
      case 2		% Ordinary Kriging, follow Journel and Huijbregts, p. 307
        K=[K0 ones(nk,1);ones(1,nk) 0];
      case 3		% Universal kriging w/ linear drift, follow Journel and Huijbregts, p. 319
        K=[K0 ones(nk,1) x(indx_sort(k))' y(indx_sort(k))'; ...
               ones(1,nk) pad; delx(j,indx_sort(k)) pad; dely(j,indx_sort(k)) pad];
    end
	 [U,S,V]=svd(K);
	 kindx=find(diag(S) > EPS);
	 Sinv=diag(1./diag(S(kindx,kindx)));
	 K_inv=V(:,kindx)*Sinv*U(:,kindx)';
  end			
  if nk == 1 & mod_opt == 1	% special case for matrix  0/0
		lambda=1;
  else
  		lambda=K_inv*M2;
  end
  indx=find(isnan(M2));
  if length(indx) > 0 
     disp(sprintf('nans in lambda at j=%g ',j));
  end
  lnan=find(isnan(var(indx_sort(k))));
  if length(lnan) > 0
     var(k(lnan))=zeros(length(lnan),1);
  end
  tp(j)=sum_nan(lambda(1:nk).*var(indx_sort(k)));
  if blk_opt == 2 & model < 0
     ep(j)=CVV-sum_nan(lambda.*M2);	% Error estimate for tp(j)
  else
     if model_type == 1
        ep(j)=1-sum_nan(lambda.*M2);
     else
        ep(j)=sum_nan(lambda.*M2);
     end
  end
  if abs(tp(j)) > 3.0*val_m  
    Anomaly_cnt=Anomaly_cnt+1;
    tp(j)=val_m;
  end
  if abs(ep(j)) > eplim 
    tp(j)=NaN;
    ep(j)=abs(ep(j));
  end
end

if NaNcnt > 1 | Anomaly_cnt > 1
   disp(sprintf('NaNcnt = %g, Anomaly_cnt = %g',NaNcnt,Anomaly_cnt))
end
percent=sprintf('   %g percent done in %g sec.',100,round(toc));
if ~isempty(hmsg) krig_message(4,percent,' ',hmsg); end
