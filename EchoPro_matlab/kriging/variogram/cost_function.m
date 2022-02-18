function f = cost_function(model_para,data_in)
%% function f = cost_function(model_para,data_in) constructs a cost function
% 	f = fun(para,data)
%
% Input: 
%     model_para  = [nugt sill lscl powr hole]
%	         data_in = [distance variance weighting function]
%
% Output:  
%%    f =cost function
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global para

  x = data_in(:,1);							% lag
  y = data_in(:,2);							% variogram computed form observed data
  
% compute selected model with parameter 'para'
 yr = variogrammodel3d(para.vario.model,x,model_para);
 w = data_in(:,3);   % weigthing function or proportionality coef. in least-square sense
 f = (yr(:) - y).*w;	% cost function

