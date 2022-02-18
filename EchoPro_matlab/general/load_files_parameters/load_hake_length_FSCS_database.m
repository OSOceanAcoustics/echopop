function     out=read_hake_length_FSCS_database(filename,species_code_id,haul_num_offset)
% read hake length data in FSCS data format
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       3/24/2013

j=0;
i=0;
fid=fopen(filename,'r');
line=fgetl(fid);                % information line
line=fgetl(fid);
if nargin < 3
    haul_num_offset=0;
end
while line > 0
% process the first line of any trawl
   ind=find(line == ',');
   if (str2num(line(ind(2)+1:ind(3)-1)) == species_code_id) | (str2num(line(ind(2)+1:ind(3)-1)) == 100000+species_code_id)
      trawl_num=str2num(line(ind(1)+1:ind(2)-1))+haul_num_offset;
      if i~= 0 & out(i).trawl_no == trawl_num
         j=j+1;
      else
        j=1;
        i=i+1;
        out(i).trawl_no=trawl_num;
      end
      out(i).length(j)=str2num(line(ind(7)+1:ind(8)-1));
      out(i).sex(j)=line(ind(11)+1:ind(12)-1);   
   end
   line=fgetl(fid);
end
for j=1:length(out);
    len=out(j).length;
    out(j).meanLen=mean(len);
    out(j).stdLen=std(len);
    TS0=20*log10(len)-68;
    TS=10*log10(mean(10.^(TS0/10)));
    out(j).TS_lin=TS;                   % 10*log10(<sv>
    out(j).TS_sd=std(TS0);
    out(j).TS_log=mean(TS0);            % <10*log10(sv)>   
    out(j).n=length(len);
    out(j).Gender(1:out(j).n)=nan;
    out(j).Male_ind=find(out(j).sex == 'M');
    out(j).Gender(out(j).Male_ind)=1;
    out(j).Female_ind=find(out(j).sex == 'F');
    out(j).Gender(out(j).Female_ind)=2;
    out(j).nM=length(out(j).Male_ind);
    out(j).nF=length(out(j).Female_ind);
end
fclose(fid);
return