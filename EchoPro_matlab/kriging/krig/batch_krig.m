function		batch_krig()
%% function		batch_krig()
%% Batch kriging on slected data files using current variogram/correlogram and kriging parameters
%%
%%  Modified on 2/4/00		to be able to include the last file on the filename list
%%									without CR at the end
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.


global para hdl data
	
if ~isfield(para.vario,'nugt')
   krig_message(1,'Need to set variogram/correlogram parameters for batch processing information');
   return;
end

if 0
if para.krig.xmin == -0.5 & para.krig.xmax == 0.5 ...
   & para.krig.ymin == -0.5 & para.krig.ymax == 0.5
   krig_message(1,'You may need to set kriging parameters for batch processing information');
   return;
end
end

if ~isfield(para.krig,'batch_data_file')
   krig_message(1,'Need to provide a file that contains a list of filenames');
   return;
end
if ~isfield(para.krig,'batch_log_file')
   krig_message(1,'Need to specify a log filename that record the batch processing information');
   return;
end

fid=fopen(para.krig.batch_data_file,'r');
if isempty(fid)
  return
end
para.krig.batch_file_proc=1;



%%% get input filenames and set output filenames   3/2/00
%file_len;
k=0;
while 1
  line = fgetl(fid);
  if ~isstr(line) | isempty(line), break, end   %%%%%%xxxxx modified 15June2004
  k=k+1;
  fname_len(k)=length(line);
  para.krig.batch_filename_in(k,1:fname_len(k))=getfilename(line,para); 
  if length(para.krig.batch_filename_in(k,:)) > fname_len(k)                %%%%%%xxxxx modified 15June2004
      dflen=length(para.krig.batch_filename_in(k,:))-fname_len(k);
      len2=length(para.krig.batch_filename_in(k,:));
      len1=fname_len(k)+1;
      para.krig.batch_filename_in(k,len1:len2)=' '*ones(1,dflen);
  end
  fname=filenametruncation(line);
  out_file=[para.file_dir.batch_log fname '.mat']; 
  file_len(k)=length(out_file);
  para.krig.batch_filename_out(k,1:length(out_file))=out_file;
end
fclose(fid);              

para.krig.bat_proc_time=0;
fid=fopen(para.krig.batch_log_file,'w');
for i=1:k
   para.krig.bat_proc_cnt=i;
   line1=sprintf('file #%3d',i);
   fprintf(fid,'%s\n',line1);
   fid1=fopen(para.krig.batch_filename_in(i,1:fname_len(i)),'r');              %%%%%%xxxxx modified 15June2004
   if fid1 < 0
     line1=sprintf('input filename %s is not found',para.krig.batch_filename_in(i,:));
     krig_message(1,line1);
   else
	  fclose(fid1);						% added on 2/5/00
     line2=sprintf('    input filename = %s',para.krig.batch_filename_in(i,:));    
     fprintf(fid,'%s\n',line2);
%% start data processing
     para.krig.data_file=para.krig.batch_filename_in(i,1:fname_len(i));
     indx=find(para.krig.data_file == '/' | para.krig.data_file == '\' | para.krig.data_file ==':');
     if ~isempty(indx)
        pindx=max(indx);
        filename=para.krig.data_file(pindx+1:length(para.krig.data_file));
     else
        filename=para.krig.data_file;
     end        
%     truncated_filename=filenametruncation(filename);
     [pathi truncated_filename ext ver]=fileparts(para.krig.data_file);
	  set(hdl.krig.fileID,'string',truncated_filename);
	  para.dataprep.fileID=truncated_filename;
	  loaddatfile(2);			% load data file from kriging window
	  datachk(2);					% calling datachk.m from kriging window
     if ~isfield(para.dataprep,'transform_index') para.dataprep.transform_index=1;end
     data.in.tv=datatransform(1,data.in.v,para.dataprep.transform_index);	% Forward Data Transformation
     krig3dmanager;
  %   if ~isfield(para.dispkrig,'trackline')   %% added 3/2/00
			para.dispkrig.trackline.type_indx=1;
			para.dispkrig.trackline.line_color=1;
			para.dispkrig.colormap_indx=11;
			para.dispkrig.num_of_contour=4;
			para.dispkrig.digits_of_contour=3;
			para.dispkrig.trackline.color_indx=8;
			para.dispkrig.trackline.size_indx=8;
			para.dispkrig.validation_model=1;
			para.dispkrig.Qcheck=0;
			para.dispkrig.JKcheck=0;
  %   end
     line3=sprintf('    output filename = %s',para.krig.batch_filename_out(i,:));
     cmd=['save ''' para.krig.batch_filename_out(i,1:file_len(i)) ''' para data'];
     eval(cmd);
     fprintf(fid,'%s\n',line3);
  end
end
fclose(fid);
