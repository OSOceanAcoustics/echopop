function validation_proc(opt)
%%  function validation_proc(opt) performs validation computation
%%   opt   =  1   Kriging map
%%            2   Kriging variance map
%%            3   Cross validation
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global data para hdl


para.dispkrig.validation_model=get(hdl.validation.model,'value');


% computing   
flag=cross_validation(para.dispkrig.validation_model,opt);
if flag < 0
   return
end


%% plotting results of cross-validation computations
if data.in.dim == 2
   posvec=[0.15 0.4 0.6 0.5];
else
   posvec=[0.15 0.4 0.45 0.5];
end

%% clear graphic window
	figure(hdl.dispkrig3d.h0)
   delete(hdl.dispkrig3d.axes1);
	hdl.dispkrig3d.axes1 = axes('Parent',hdl.dispkrig3d.h0, ...
	'Color',[1 1 1], ...
   'Tag','dispkrig3dAxes1', ...
	'Position',posvec,'visible','on');

switch	para.dispkrig.validation_model 
		case 1 		% Q-1
			nq=2;
			n=length(data.in.x)-para.krig.model;
	   	    [x,fq,Qrej,fr]=pdf_func(para.dispkrig.validation_model,n);		  
			plot(x,fq,'m','linewidth',2.0)
			hold on
			plot(-Qrej*ones(1,2),[0 fr],'k',Qrej*ones(1,2),[0 fr],'k','linewidth',2.0);
			ht=text(-Qrej+0.02,0.5*fr,'Accept Region');
			if data.out.krig.q1 <= max(x) & data.out.krig.q1 >= min(x)
		   	plot(data.out.krig.q1*ones(1,nq),linspace(0,max(fq),nq),'-r','linewidth',1.5);
				hdt=text(data.out.krig.q1+0.02,0.8*max(fq),'Computed Q_1');
				set(hdt,'fontsize',8,'fontweight','bold');
			end
			title(sprintf('Q1 = %g',data.out.krig.q1),'fontweight','bold')
			hold off
			hx=xlabel('Q_1');
			hy=ylabel('PDF(Q_1)');
			set([hx hy ht],'fontweight','bold')
		case 2 		% Q-2
			nq=2;
			n=length(data.in.x)-para.krig.model;
	        [x,fq,Qrej,fr]=pdf_func(para.dispkrig.validation_model,n);		  
			plot(x,fq,'m','linewidth',2.0)
			hold on
			plot(Qrej(1)*ones(1,2),[0 fr(1)],'k',Qrej(2)*ones(1,2),[0 fr(2)],'k','linewidth',2.0);
			ht=text(Qrej(1)+0.1,0.5*fr(1),'Accept Region');
			if data.out.krig.q2 <= max(x) & data.out.krig.q2 >= min(x)
		   	plot(data.out.krig.q2*ones(1,nq),linspace(0,max(fq),nq),'-r','linewidth',1.5);
				hdt=text(data.out.krig.q2+0.02,0.8*max(fq),'Computed Q_2');
				set(hdt,'fontsize',8,'fontweight','bold');
			end
			title(sprintf('Q2 = %g',data.out.krig.q2),'fontweight','bold')
			hold off
			hx=xlabel('Q_2');
			hy=ylabel('PDF(Q_2)');
			set([hx hy ht],'fontweight','bold')
      case 3					% double kriging
         Is=data.out.krig.Is;
			sta=1:length(data.in.tv);
			plot(sta,data.in.tv,'-o',data.out.krig.sta,Is,'+-r');
			hh1=legend('Observation','Prediction',0);
			set(hh1,'fontsize',10)
			hx=xlabel('SAMPLE NUMBER');
            hy=ylabel('OBSERVATION AND PREDICTION');
            val_mean=mean_nan(data.in.tv(data.out.krig.sta)-Is);
            val_std=std_nan(data.in.tv(data.out.krig.sta)-Is);
            ht=title(['Double Kriging: <z - z*> = ' sprintf('%g',val_mean) '  \sigma = ' sprintf('%g',val_std)]);
			set([hx hy ht],'fontweight','bold')
		case 4		% JackKnife
			Ijk=data.out.krig.Ijk;
			sta=1:length(data.in.tv);
			plot(sta,data.in.tv,'-o',data.out.krig.sta,Ijk,'+-r');
			hh1=legend('Observation','Prediction',0);
			set(hh1,'fontsize',10)
			hx=xlabel('SAMPLE NUMBER');
			hy=ylabel('OBSERVATION AND PREDICTION');
         val_mean=mean_nan(data.in.tv(data.out.krig.sta)-Ijk);
         val_std=std_nan(data.in.tv(data.out.krig.sta)-Ijk);
         ht=title(['Leave One Out: <z - z*> = ' sprintf('%g',val_mean) '  \sigma = ' sprintf('%g',val_std)]);
			set([hx hy ht],'fontweight','bold')			
	end
   figure(hdl.validation.h0)