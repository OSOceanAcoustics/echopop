function  combine_two_biological_data_output(out,out1,out2)
% combine two sets of biological data (US & CAN)
%% out = samples from station #1
%% out1 = samples from station #2
%% out2 = length, weight, sex, and trawl # for all samples in arrays
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global hdl para data

%% length_sex data
n=length(data.bio.hake_length_sex);  % number of US trawl length data at station #1
n1=length(out);                      % number of CAN trawl length data at station #1
data.bio.hake_length_sex(n+1:n+n1)=out;  % combine US and CAN stations #1 length data

%%% length-weight_sex_age data
n=length(data.bio.hake_length_weight_sex_age); % number of US trawl sample data at station #2
n1=length(out1);                               % number of CAN trawl sample data at station #2
data.bio.hake_length_weight_sex_age(n+1:n+n1)=out1; % combine US and CAN stations #2 data

%% all length data from station #1 and 2 of US and CAN trawl data 
field_names=fieldnames(data.bio.len_wgt_all);
n2=length(field_names);
for i=1:n2
    cmd=['data.bio.len_wgt_all.' char(field_names(i)) '=[data.bio.len_wgt_all.' char(field_names(i)) ' out2.' char(field_names(i)) '];'];
    eval(cmd)
end

