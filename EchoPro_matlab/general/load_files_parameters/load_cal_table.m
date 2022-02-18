function out=load_cal_table(filename,freq_ind)
% load calibration parameters from a file
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global para

%% pulse length table 18, 38, 70, 120, 200 kHz
%%                              18 kHz                  38 kHz                      70 kHz              120 kHz             200 kHz
pulselengthtable=[512 1024 2048 4096 8192; 256 512 1024 2048 4096; 128 256 512 1024 2048; 64 128 256 512 1024; 64 128 256 512 1024];
fid=fopen(filename);
if fid < 0
    disp('Calibration Table File %s not found !!');
    out=-1;
    return
end
fclose(fid);

dat=load(filename);
gain=dat(1,:);
phi_2way=dat(2,:);
alpha=dat(3,:);
tx_power=dat(4,:);
cw=dat(5,:);
D_offset=dat(6,:);
if size(dat,1) > 6 
   pulse_duration=dat(7,:)*1e-6;
   sacorrectiontable=dat(8,:);
end

for i=1:length(freq_ind);
   ind=freq_ind(i);
   out(i).gain=gain(ind);
   out(i).equivalentbeamangle=phi_2way(ind);
   out(i).absorptioncoefficient=alpha(ind);  
   out(i).transmitpower=tx_power(ind);
   out(i).soundvelocity=cw(ind);
   out(i).transducerdepth=D_offset(ind);
   if size(dat,1) > 6 
      out(i).pulselength=single(pulse_duration(ind));
      pulse_indx=find(abs(pulselengthtable(ind,:) - out(i).pulselength) < 1e-6);
      out(i).sacorrectiontable(pulse_indx)=sacorrectiontable(ind);
   end
   if para.proc.source == 2 & para.acoust.file_sys == 2 &  strcmp(para.survey_year,'2003')
      out(i).sampleresolution=0.1*out(i).soundvelocity/1500;                      % 2003 ek500 sampling interval with correction
 %     out(i).sampleresolution=0.1;                      % 2003 ek500 sampling interval w/o correction 
      out(i).sampleinterval=2*out(i).sampleresolution/out(i).soundvelocity;   % sampling distance of 10 cm (1480 m/s)~ 7400 Hz sampling rate
   else
      out(i).sampleinterval=out(i).pulselength/4;
      out(i).sampleresolution=cw(ind)*out(i).pulselength/8;
   end
end

return