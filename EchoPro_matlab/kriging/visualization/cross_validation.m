function			out=cross_validation(opt1,opt2)
% function	cross_validation(opt) performs cross validation computation
%	opt1 =  1				Q1 cross validation
%			2				Q2 cross validation
% 		 	3				double kriging
%			4				Jackknife
%   opt2 =  1				Compute if it has not computed before
%			2				Re-compute use the current variogram parameter settings
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global	data para hdl

out=1;
if ~isfield(para,'vario')  
   krig_message(1,'Need to performing variogram first!!');
   out=-1;
   return
elseif ~isfield(para.vario,'model')
   krig_message(1,'Need to performing kriging first!!');
   out=-1;
   return
end


vmod=para.vario.model;
if para.vario.corr == 1
  vmod = -vmod;
end
model_para=[para.vario.nugt para.vario.sill para.vario.lscl para.vario.powr para.vario.hole];

%% kriging parameters
mod_opt=para.krig.model;
srad=para.krig.srad;
kmin=para.krig.kmin;
kmax=para.krig.kmax;
blk_opt=para.krig.scheme;
blknx=para.krig.blk_nx;
blkny=para.krig.blk_ny;
elim=para.krig.elim;
dx=para.krig.dx;
dy=para.krig.dy;
EPS=para.krig.eps;
var=data.in.tv;
x=data.in.x;
y=data.in.y;
if data.in.dim == 3
   dz=para.krig.dz;
   z=data.in.z;
   indx_nan=find(isnan(x)|isnan(y)|isnan(z));
   if ~isempty(indx_nan)
      z(indx_nan)=[];
   end
else
   indx_nan=find(isnan(x)|isnan(y));
end
if ~isempty(indx_nan)
   indx_orig=[1:length(x)]';
   data.out.krig.sta=setdiff(indx_orig,indx_nan);
   x(indx_nan)=[];
   y(indx_nan)=[];
   var(indx_nan)=[];
else
   data.out.krig.sta=1:length(x);
end

msg1='Cross-Validation in progress, ...';
msg2=' ';
nmin=8;						% standard for displaying time elapse

if ~isfield(para.krig,'model')
   krig_mod=2;							% default kriging method: ordinary Kriging
else
   krig_mod=para.krig.model;
end
n=length(x);

tic
if opt1 <= 2
	%			% Q-1
	 if para.dispkrig.Qcheck == 0 | opt2 == 1
		sigma0=sqrt(para.vario.c0);
		n=length(x);
		if n >= nmin hmsg=krig_message(2,msg1,' ');drawnow;end
		nv=n-mod_opt;
		nv_cnt=max(1,round(0.01*nv));
		for i=mod_opt+1:n
			j=i-mod_opt;
  			if rem(j,nv_cnt) == 0 
     			percent=sprintf('   %g percent done in %g sec.', floor(100*j/nv),round(toc));
     		   if n >= nmin krig_message(4,percent,' ',hmsg);drawnow; end
        end
        if data.in.dim == 2
			xc=x(i);yc=y(i);
			xi=x(1:i-1);yi=y(1:i-1);var_in=var(1:i-1);
			if krig_mod == 1			% remove the mean for simple kriging
				var_mean=mean(var_in);
				var_in=var_in-var_mean;
   			end
   %         para.krig.kmin=mod_opt;
   %         para.krig.kmax=length(xi);
   %         para.krig.srad=0.3;
			[var_out,err_krig]=krig2d(xc,yc,xi,yi,var_in,model_para);
  			if para.krig.model == 1			% add in the mean for simple kriging
				var_out=var_out+var_mean;
			end
         else                       % 3-D kriging
			xc=x(i);yc=y(i);zc=z(i);
			xi=x(1:i-1);yi=y(1:i-1);zi=z(1:i-1);var_in=var(1:i-1);
			if krig_mod == 1			% remove the mean for simple kriging
				var_mean=mean(var_in);
				var_in=var_in-var_mean;
   			end
            [var_out,err_krig]=krig3d(xc,yc,zc,xi,yi,zi,var_in,model_para);
  			if para.krig.model == 1			% add in the mean for simple kriging
				var_out=var_out+var_mean;
			end
         end    
		 if err_krig < 0 err_krig=nan;end
		 ek(i-mod_opt)=(var(i)-var_out)/(sigma0*sqrt(err_krig));
         errk(i-mod_opt)=err_krig;
        end
        data.out.krig.ek=ek;
	    data.out.krig.q1=mean_nan(ek);
		data.out.krig.q2=mean_nan(ek.^2);
		data.out.krig.ek=ek;
		if n > nmin close(hmsg);end
	 else
        data.out.krig.q1=mean_nan(data.out.krig.ek);
		data.out.krig.q2=mean_nan(data.out.krig.ek.^2);
	 end
    para.dispkrig.Qcheck=1;
    set(hdl.validation.recompute,'enable','on');
elseif opt1 ==  3			% double kriging - default
    if para.dispkrig.Dblkrig == 0 | opt2 == 1
	   n=length(x);
	   if n >= nmin hmsg=krig_message(2,msg1,' ');end
       if data.in.dim == 2
		    xc=x;yc=y;
            xi=data.out.krig.xp;
            yi=data.out.krig.yp;
            var_in=data.out.krig.var;
			if krig_mod == 1			% remove the mean for simple kriging
				var_mean=mean(var_in);
				var_in=var_in-var_mean;
   		    end
			[var_out,err_krig]=krig2d(xc,yc,xi,yi,var_in,model_para,hmsg);
			if para.krig.model == 1			% add in the mean for simple kriging
				var_out=var_out+var_mean;
			end
       else
			xc=x;yc=y;zc=z;
         xi=data.out.krig.xp;
         yi=data.out.krig.yp;
         zi=data.out.krig.zp;
         var_in=data.out.krig.var;
			if krig_mod == 1			% remove the mean for simple kriging
				var_mean=mean(var_in);
				var_in=var_in-var_mean;
   		end
         [var_out,err_krig]=krig3d(xc,yc,zc,xi,yi,zi,var_in,model_para,hmsg);
			if para.krig.model == 1			% add in the mean for simple kriging
				var_out=var_out+var_mean;
			end
       end          
       data.out.krig.Is=var_out(:);  
	   if n > nmin close(hmsg);end
    end
    para.dispkrig.Dblkrig=1;
    set(hdl.validation.recompute,'enable','on');
elseif opt1 == 4			% jackknife
	 if para.dispkrig.JKcheck == 0 | opt2 == 1
		var=data.in.tv;
		n=length(x);
		if n >= nmin hmsg=krig_message(2,msg1,' ');end
        Ijk(1)=var(1);
		nv=n-mod_opt;
		nv_cnt=max(1,round(0.01*nv));
		for i=2:n
            j=i-mod_opt;
  			if rem(j,nv_cnt) == 0 
     			percent=sprintf('   %g percent done in %g sec.', floor(100*j/nv),round(toc));
     		   if n >= nmin krig_message(4,percent,' ',hmsg);drawnow; end
  			end
            if data.in.dim == 2
				xc=x(i);yc=y(i);
				xi=[x(1:i-1); x(i+1:n)];yi=[y(1:i-1); y(i+1:n)];var_in=[var(1:i-1); var(i+1:n)];
				if krig_mod == 1			% remove the mean for simple kriging
					var_mean=mean(var_in);
					var_in=var_in-var_mean;
   			    end
                [var_out,err_krig]=krig2d(xc,yc,xi,yi,var_in,model_para);
   				if para.krig.model == 1			% add in the mean for simple kriging
					var_out=var_out+var_mean;
				end
            else
				xc=x(i);yc=y(i);zc=z(i);
				xi=[x(1:i-1); x(i+1:n)];yi=[y(1:i-1); y(i+1:n)];zi=[z(1:i-1); z(i+1:n)];var_in=[var(1:i-1); var(i+1:n)];
 				if krig_mod == 1			% remove the mean for simple kriging
					var_mean=mean(var_in);
					var_in=var_in-var_mean;
   			    end
                [var_out,err_krig]=krig3d(xc,yc,zc,xi,yi,zi,var_in,model_para);
  				if para.krig.model == 1			% add in the mean for simple kriging
					var_out=var_out+var_mean;
				end
            end
			Ijk(i)=var_out;
	    end
		data.out.krig.Ijk=Ijk(:);
		if n > nmin close(hmsg);end
	 end
	 para.dispkrig.JKcheck=1;
    set(hdl.validation.recompute,'enable','on');
end