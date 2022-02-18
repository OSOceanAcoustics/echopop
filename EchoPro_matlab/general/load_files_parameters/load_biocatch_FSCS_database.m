function   out=load_biocatch_FSCS_database(filename)
% load biocatch data in FSCS data format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013


fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Catch file: %s !!',filename));
    return
end
%out.filename=filename;
line=fgetl(fid);            % information line

i=0;
j=0;
line=fgetl(fid);
while line > 0
% process the first line of any trawl
   ind=find(line == ',');
   trawl_num=str2num(line(ind(1)+1:ind(2)-1));
   if i ~= 0 & out(i).trawl_no == str2num(line(ind(1)+1:ind(2)-1))
           j=j+1;
   else
      j=1;
      i=i+1;
      out(i).trawl_no=trawl_num;
   end
   out(i).species(j).ID=line(ind(2)+1:ind(3)-1);
   num1a=str2num(line(ind(4)+1:ind(5)-1));
   num2a=str2num(line(ind(5)+1:ind(6)-1));
   out(i).species(j).wgt=str2num(line(ind(6)+1:ind(7)-1));
   out(i).species(j).exp_wgt=str2num(line(ind(7)+1:ind(8)-1));   
   out(i).species(j).name=line(ind(8)+1:ind(9)-1);
   out(i).species(j).sub_sample=line(ind(9)+1:ind(10)-1);
   num1b=str2num(line(ind(11)+1:ind(12)-1));
   num2b=str2num(line(ind(12)+1:ind(13)-1));
   if isempty(num1a) & isempty(num1b)
       out(i).species(j).cnt=[];
   elseif isempty(num1b)
       out(i).species(j).cnt=num1a;
   elseif isempty(num1a)
       out(i).species(j).cnt=num1b;
   else
       out(i).species(j).cnt=max(num1a,num1b);
   end
   if isempty(num2a) & isempty(num2b)
       out(i).species(j).exp_cnt=[];
   elseif isempty(num2b)
       out(i).species(j).exp_cnt=num2a;
   elseif isempty(num2a)
       out(i).species(j).exp_cnt=num2b;
   else
       out(i).species(j).exp_cnt=max(num2a,num2b);
   end
   line=fgetl(fid);
end

fclose(fid);