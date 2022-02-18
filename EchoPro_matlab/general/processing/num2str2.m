function out=num2str2(num_in,str_len)
% convert number to string having a specified string length with added '0's in front
% 

val=num2str(num_in(:));
[n,m]=size(val);


for i=1:length(num_in) 
  str_zeros='';
  val=num2str(num_in(i));
  for j=1:str_len-length(val)
     str_zeros=[str_zeros '0']; 
  end
  out(i,1:str_len)=[str_zeros val];
end
return