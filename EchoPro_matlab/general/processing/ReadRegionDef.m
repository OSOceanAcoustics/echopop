%% read region definition file exported from EchoView for processing raw data instead of the exported
%% NASC data from EchoView
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

function Region=ReadRegionDef(filename)

fid=fopen(filename,'r');
if fid < 0
   Region=-1;
   error('Region definition file cannot be found!!')
   return
end

file_ver=fgetl(fid);
n_reg=str2num(fgetl(fid));

for i=1:n_reg
    if i == 4
%        disp(i)
    end
    Region(i).seq_num=i;
    D=fscanf(fid,'%d%d%d%d%d%d%d%d%f%f%d%f%f');
    Region(i).ver=D(1);                             % region defination version
    Region(i).npts=D(2);                            % number of data points defining the region
    Region(i).n=D(3);                               % region sequence number
    Region(i).creation_type=D(5);                % region creation type (rec, polygons, etc.)
    Region(i).data_valid=D(7);                      % flag of data validation for the next 4 fields
    if Region(i).data_valid == 1
        Region(i).min_date=D(8);
        Region(i).min_time=round(D(9));
        Region(i).min_dep=D(10);
        Region(i).max_date=D(11);
        Region(i).max_time=round(D(12));
        Region(i).max_dep=D(13);
    end
    tmp1=fgetl(fid);
    id=strfind(tmp1,' log ');
    if ~isempty(id)
       Region(i).notes={tmp1(1:id-2)};
       Region(i).log=str2num(tmp1(id+4:end));
       Region(i).class={fgetl(fid)};
       d=fscanf(fid,'%f');
    else
       tmp2=fgetl(fid);
       d=str2num(tmp2);
       Region(i).log='none';
       if length(d) > 0         
          Region(i).notes='';
          Region(i).class={tmp1};
       else
          tmp3=fgetl(fid);
          d=str2num(tmp3);          
          Region(i).notes={tmp1};
          Region(i).class={tmp2};
       end
       
    end
    
    for j=1:Region(i).npts
        Region(i).date(j)=d(3*(j-1)+1);
        Region(i).time(j)=d(3*(j-1)+2);
        Region(i).dep(j)=d(3*(j-1)+3);
    end
% remove dulplicate points
    diff_t=abs(diff(Region(i).time));
    diff_d=abs(diff(Region(i).dep));
    ind=find(diff_t == 0 & diff_d < 1e-3);
    Region(i).time(ind)=[];
    Region(i).dep(ind)=[];
    Region(i).npts=Region(i).npts-length(ind);
    Region(i).type=d(end);
    Region(i).name=fgetl(fid);
    if feof(fid)
        fclose(fid);
        return
    end
end
fclose(fid);
return
    