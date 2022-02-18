function		[x,pdf,Qj,fQ]=pdf_func(opt,n)
%% function		[x,pdf,Qj,fQ]=pdf_func(opt,n) computes the PDF based on the parameter opt
%% INPUT:
%%    opt = task option
%%        1 --  Q-1 cross-validation
%%        2 --  Q-2 cross validation
%%      n = number of opints in PDF output
%% OUPUT:
%%    x, pdf - variable and its PDF
%%    Qj, fQ - x values corresponding to 0.0275 and 0.0975 of the
%%             probability function  int_(-inf)^(x) PDF(t) dt
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


switch opt
	case 1		% Q-1 model (Gaussian)
		m=50;
		xb=5/sqrt(n);				% 5 times standard deviation
		x=linspace(-xb,xb,m);
		pdf=exp(-0.5*n*x.^2)*sqrt(0.5*n/pi);
		Qj=2/sqrt(n);
		fQ=exp(-2)*sqrt(0.5*n/pi);
	case 2		% Q-2 model (kappa 2)
		m=800;
		xb=max(2,10/sqrt(n));				% > 10 times standard deviation
		x=linspace(0,xb,m);
		dx=x(2)-x(1);
		if n < 100     % chi-square with n degree of freedom
	  	  pdf=(0.5*n)^(n/2)*exp(-n*x/2).*x.^(n/2-1)/(gamma(n/2));
		else           % Gaussian with mean = 1, and std = 2/sqrt(n)
		  pdf=exp(-0.25*n*(x-1).^2)*sqrt(0.25*n/pi);
		end
		S=cumsum(pdf)*dx;
		ft=[0.025 0.975];
		[var indx1]=min(abs(S-ft(1)));
		[var indx2]=min(abs(S-ft(2)));
		Qj(1)=x(indx1);
		Qj(2)=x(indx2);
		fQ=pdf([indx1 indx2]);
end