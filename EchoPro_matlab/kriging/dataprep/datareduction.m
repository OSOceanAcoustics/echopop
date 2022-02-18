function   datareduction(opt)
%% data reduction or subsampling
%% opt = 1					get filter parameters from the window panel
%%		 2					get filter parameters from variable struct
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global  hdl data para

if opt == 1
  para.dataprep.reduct_fac=str2num(get(hdl.dataprep.reduct_fac,'string'));
  para.dataprep.filter_supt=str2num(get(hdl.dataprep.filter_supt,'string'));
  para.dataprep.filter_type=get(hdl.dataprep.filter,'value');				% filter type: 1-simple reduction,
																								%	2-mean (box car)
																								% 	3-median
end

ninc=para.dataprep.reduct_fac;
filter_type=para.dataprep.filter_type;
supt=para.dataprep.filter_supt;

x=data.in.x;
y=data.in.y;
z=data.in.z;
var=data.in.v;

if ninc > 1 | supt > 1
 n=length(x);
 if filter_type >= 2			% not simple data reduction
   k=1;
   if rem(ninc,2) ~= 0  & ninc ~= 1 % odd ninc
      ninc=ninc+1;
   end   
   if rem(supt,2) ~= 0   % odd support
      supt=supt-1;
   end   
   for i=1:ninc:n
      if i <= supt/2
         indx=1:i+supt/2;
      elseif i > n-supt
         indx=i-supt/2:n;
      else
         indx=i-supt/2:i+supt/2;
      end
      switch filter_type 
        case 2 % box car (mean)
           var_new(k)=mean_nan(var(indx));
        case 3	% median filter to avoid outliers
          var_new(k)=median_nan(var(indx));
      end
   % data position should be that of one of the data involved points
      [val, indxm]=min(abs(var_new(k)-var(indx)));
      xnew(k)=x(indx(indxm));		
      ynew(k)=y(indx(indxm));
      if data.in.dim == 3
      	znew(k)=z(indx(indxm));
      end        
      k=k+1;
   end
   x=xnew(:);
   y=ynew(:);
   if data.in.dim == 3
      z=znew(:);
   end
   var=var_new(:);
 else							% simple data reduction or subsampling
   x=x(1:ninc:n);
   y=y(1:ninc:n);
   var=var(1:ninc:n);
   if data.in.dim == 3
      z=z(1:ninc:n);
   end
 end
end

data.in.x=x;
data.in.y=y;
data.in.z=z;
data.in.v=var;