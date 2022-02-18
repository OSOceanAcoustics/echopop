function y = variogrammodel3d(type, r, model_para)
%  function y = variogrammodel3d(type, r, model_para) computes the
%  semi-variogram/correlogram based on the given model_parameters
%   y = variogrammodel ( type, r, model_para)
%   |type| = model index for semivariogram/correlogram
%   type > 0: semi-variogram
%        < 0: correlogram  
%      r = vector lag distances
%  model_para = [ p1 p2 p3 p4 p5];
%      
%      p1: Nugt = nugget effect
%      p2: Sill = sill
%      p3: L = length scale for the main lobe
%      p4: p = power for the expenential 
%      p5: b = length scale for hole effect     
%
%      model type:
%          01 = spherical
%          02 = exponential
%          03 = gaussian
%          04 = linear

%      models including hole effects
%          05 = C * [ 1 - (sin b*r) / r ] + Nugt
%          06 = C * [ 1 - (exp(-r/L))     * cos(br) ] + Nugt
%          07 = C * [ 1 + (exp(-r/L))     * cos(br) ] + Nugt
%          08 = C * [ 1 - (exp(-(r/L)^2)) * cos(br) ] + Nugt
%          09 = C * [ 1 -                   Jo (br) ] + Nugt
%          10 = C * [ 1 -  exp(-r/L)      * Jo (br) ] + Nugt
%          11 = C * [ 1 -  exp(-(r/L)^2)  * Jo (br) ] + Nugt
%	        12 = C * [ 1 -  exp(-(r/L)^2)  * (1 - br)] + Nugt
%          13 = C * [ 1 -  exp(-(r/L)^p)  * Jo (br) ] + Nugt	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Marcotte
% Version 2.1  97/aug/18
% Revised by Dezhang Chu,   10-29-98
%%
%%  Kriging Software Package  version 3.0,   December 29, 2001
%%  Copyright (c) 1998, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


Nugt=model_para(1);
Sill=model_para(2);
L=model_para(3);
p=model_para(4);
b=model_para(5);	
C=Sill-Nugt;
n = length(r);
rL = r ./ L;

switch abs(type) 
case 1
     indx1=find(rL < 1);
     indx2=find(rL >= 1);
     rL1=rL(indx1);
     rL2=rL(indx2);
     y1 = C * ( 1.5 .* rL1  -  0.5 .* rL1.^3)+Nugt;
     y2 = Sill * ones(size(rL2));
     y(indx1)=y1;
     y(indx2)=y2;
     y=reshape(y,size(r,1),size(r,2));
 case 2
     y = C * ( 1 - exp(-(r/L)))+Nugt;
 case 3
   y = C * ( 1 - exp(-(r/L).^2))+Nugt;
 case 4
   y = C .* r+Nugt;
 case 5
    y = C .* ( 1 - sin(b.*(r+eps))./(r+eps) )+Nugt;
 case 6
   y = C * ( 1 - exp(-r/L) .* cos(b*r) )+Nugt;
 case 7
   y = C * ( 1 + exp(-r/L) .* cos(b*r) )+Nugt;
 case 8
   y = C * (1 - exp(-(r/L).^2) .* cos(b*r) )+Nugt;
 case 9
   y = C * ( 1 - besselj(0, b*r) )+Nugt;
 case 10
   y = C * ( 1 - besselj(0, b*r) .* exp(-r/L) )+Nugt;
 case 11
   y = C * ( 1 - exp(-(r/L).^2) .* besselj(0, b*r) )+Nugt;
 case 12
    y = C * ( 1 - exp(-(r/L).^2) .* (1 - b*r.^2) )+Nugt;
 case 13							% generalized exponetial-Bessel 
    y= C* (1 - exp(- (r/L).^p).* besselj(0, b*r))+Nugt;		
end
if type < 0
  if type == -1 | type == -4
	y=1-y;
  else
   y=(1-(y-Nugt)/C)*C;
  end
end