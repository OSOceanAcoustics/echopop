% get highest plane for better visualization colormap with TrueColar 
function   getXvisual

global para

!xdpyinfo > info.m

para.Xvisual=[];
fid=fopen('color_info.dat');
m=1;
if fid >= 0
%%% find xVisual
  while 1
    line = fgetl(fid);
    if ~isstr(line) | (isempty(line) & m > 1), break, end   
    i=strfind(line,'visual:');
%    disp(sprintf('%d, %s',k,line))
    if ~isempty(i)
          line = fgetl(fid);
          ii=strfind(line,'0x'); 
          Hex=line(ii:ii+3);
          line = fgetl(fid);
          iii=strfind(line,'TrueColor');  
          if ~isempty(iii)
            XvisualHex(m,1:4)=Hex;
            line = fgetl(fid);
            indx=strfind(line,'planes');    
            planes(m)=str2num(line(indx-3:indx-1));
            m=m+1;
          end
    end
  end
  fclose(fid);   
  if m > 1
    [var,indx]=max(planes);
    para.Xvisual=XvisualHex(indx,:);
    para.Xvisual_planes=planes(indx);
  end
end
disp(m)
if isempty(para.Xvisual)
    disp('not found!!')
end