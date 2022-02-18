function  s=merge_two_struct(s1,s2)
% merge two structures to s, i.e. combine the fields in s1 and s2 into s
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

s=s1;
fld1=fieldnames(s1);
fld2=fieldnames(s2);

for i=1:length(fld1)
    for j=1:length(fld2)
        if strcmp(fld1(i),fld2(j))
            cmd=['s.' char(fld1(i)) '=[s.' char(fld1(i)) '; s2.' char(fld2(j)) '];'];
            eval(cmd)
            break
        end
    end
end

for i=1:length(fld2)
    n=0;
    for j=1:length(fld1)
        if ~strcmp(fld2(i),fld1(j))
            n=n+1;
            continue
        end
    end
    if n == length(fld1)
        cmd=['s.' char(fld2(i)) '=s2.' char(fld2(i)) ';'];
        eval(cmd)
    end
end

return