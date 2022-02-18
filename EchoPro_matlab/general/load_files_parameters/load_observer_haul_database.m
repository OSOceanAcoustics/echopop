function out=load_observer_haul_database(filename)
% load observer haul data
% 
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       10/23/2018

global para data

fid=fopen(filename,'r');
if fid < 0
    out=-1;
    disp(sprintf('\n Could not open Bio-Haul file: %s !!',filename));
    return
end
fclose(fid);

[dat0,dat_str]=xlsread(filename);

if str2num(para.survey_year) > 2007
    dat = dat0;
    %% generate ship_ID according to the ascending numbers
    ship_id = cell2mat(dat_str(2:end,1));
    ship_haul_id0 = ship_id(:, [1:11 17:20]);
    ship_haul_id_num0 = int64(str2num(ship_haul_id0));
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(ship_haul_id_num0);
    ship_haul_id = ship_haul_id0(sort_ind,:);
    dat(:, 1) = ship_haul_id_num;
    dat(:, 2:end) =  dat0(sort_ind, 2:end);
    ship_haul_id_num = str2num(ship_haul_id);
    dat_str = dat_str(sort_ind+1,:);
else
    %% sort to ascending order
    [ship_haul_id_num, sort_ind] = sort(dat0(:,1));
    dat = dat0(sort_ind,[1 3:end]);
    ship_haul_id = num2str(ship_haul_id_num);
end


if para.proc.exclude_age1 == 1
    %% exclude age-1 hauls
    [intersect_hauls,IA,IB]= intersect(dat(:,1),para.proc.age1_haul);
    if ~isempty(intersect_hauls)
        [selected_hauls,IA,IB]= setxor(dat(:,1),intersect_hauls);
        ind=[];
        for i=1:length(IA)
            ind0=find(dat(:,1) == selected_hauls(i));
            ind=[ind; ind0];
        end
        dat=dat(ind,:);
    end
end

out.trawl_no=dat(:,1);

out.haultype=[];
out.performance=dat(:,8);
out.duration=dat(:,11);
out.distance=[];
out.stratum=[];
out.ship_haul_id = ship_haul_id;
if para.proc.hemisphere(1) == 'N'
   out.EQlat=abs(mean(dat(:,[3 5]),2));
   out.EQlon=abs(mean(dat(:,[4 6]),2));
elseif para.proc.hemisphere(1) == 'S'
   out.EQlat=-abs(mean(dat(:,[3 5]),2));
   out.EQlon=-abs(mean(dat(:,[4 6]),2));
else
   fprintf('Wrong N/S Hemisphere definition!!\n');
end
if para.proc.hemisphere(2) == 'W'
   out.EQlon=-abs(mean(dat(:,[4 6]),2));
elseif para.proc.hemisphere(2) == 'E'
   out.EQlon=abs(mean(dat(:,[4 6]),2));
else
   fprintf('Wrong E/W Hemisphere definition!!\n');
end
out.ave_dep=dat(:,9)*6*0.0254*12;   % fathoms to meter
out.ave_bot_dep=dat(:,10)*6*0.0254*12;   % fathoms to meter

return