%% sort exported csv files from EV files by transect numbere
function   file_indx=sort_files_by_transect_num(files)
% INPUTS: 
%     files: interval file struct
%
% OUTPUT:
%     file_indx:  file sorting index 


n=length(files);

for i=1:n
    ind1=strfind(files(i).name,'-T')+2;
    ind2=strfind(files(i).name,'-Z')-1;
    transect(i)=str2num(files(i).name(ind1:ind2));
end
[ordered_file,file_indx]=sort(transect);
% figure(1)
% plot(ordered_file,'.-')
