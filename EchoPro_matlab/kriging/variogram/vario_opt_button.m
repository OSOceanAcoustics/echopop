function  vario_opt_button()
% function  vario_opt_button() enables and disables the
% the parameter option buttons depending on the selected model
%%
%%  Kriging Software Package  version 3.0,   December 29, 2001
%%  Copyright (c) 1998, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl para

model_str={'1.  sphere','2.  exponential','3.  gaussian','4.  linear','5. sinc','6. exponential-cosine (type I)', ...
    '7. exponential-cosine (type II)','8. gaussian-cosine','9. bessel (1-Jo)','10. exponential-bessel', ...
    '11. gaussian-bessel','12. gaussian-linear','13. general exponential-bessel'};

%%% read in parameters
model=get(hdl.vario.model,'Value');
switch model 
	case 1			% spherical
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='off';
   case 2			% exponential
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='off';
   case 3			% gaussian
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='off';
   case 4			% linear
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='off';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='off';
   case 5			% sinc
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='off';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 6			% exponential-cosine (type I)
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	 POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 7			% exponential-cosine (type II)
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 8			% gaussian-cosine
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 9			% bessel (1-Jo)
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='off';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 10			% exponential-bessel
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 11			% gaussian-bessel
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 12			% gaussian-linear
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='off';
     HOLE_ON_OFF='on';
   case 13			% generalized exponetial-Bessel 
     NUGT_ON_OFF='on';
     SILL_ON_OFF='on';
     LSCL_ON_OFF='on';
	  POWR_ON_OFF='on';
     HOLE_ON_OFF='on';
end
%% value editing window
set(hdl.vario.nugt_edit,'enable',NUGT_ON_OFF);
set(hdl.vario.sill_edit,'enable',SILL_ON_OFF);
set(hdl.vario.lscl_edit,'enable',LSCL_ON_OFF);
set(hdl.vario.powr_edit,'enable',POWR_ON_OFF);
set(hdl.vario.hole_edit,'enable',HOLE_ON_OFF);

%% slider
set(hdl.vario.nugt_slider,'enable',NUGT_ON_OFF);
set(hdl.vario.sill_slider,'enable',SILL_ON_OFF);
set(hdl.vario.lscl_slider,'enable',LSCL_ON_OFF);
set(hdl.vario.powr_slider,'enable',POWR_ON_OFF);
set(hdl.vario.hole_slider,'enable',HOLE_ON_OFF);
