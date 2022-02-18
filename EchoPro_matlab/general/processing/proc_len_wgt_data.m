function   proc_len_wgt_data
%% process length weight data (all hake trawls) to obtain 
%% (1) length-weight regression or length-weight-key
%% (2) length-age keys
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013

global  data para

ni=length(data.bio.hake_length_weight_sex_age);     % number of trawls
n=length(data.bio.len_wgt_all.len);                 % number of length-weight measurements
len=6:90;                                           % length bins to get the theoretical regression curve


%% length-weight regression for all trawls (male & female)
ind=find(~isnan(data.bio.len_wgt_all.len) & ~isnan(data.bio.len_wgt_all.wgt));
x=log10(data.bio.len_wgt_all.len(ind));
y=log10(data.bio.len_wgt_all.wgt(ind));
p=polyfit(x,y,1);               % linear regression
w0=10^p(2);
w=w0*len.^p(1);
data.bio.len_wgt_all.reg_w0=w0;      % w=w0*len^p;
data.bio.len_wgt_all.reg_p=p(1);

%% %% length-weight regression for all trawls for male
indM=find(data.bio.len_wgt_all.sex(ind) == 1);
xM=log10(data.bio.len_wgt_all.len(ind(indM)));
yM=log10(data.bio.len_wgt_all.wgt(ind(indM)));
pM=polyfit(xM,yM,1);               % linear regression
w0M=10^pM(2);
wM=w0M*len.^pM(1);
data.bio.len_wgt_all.reg_w0M=w0M;      % w=w0*len^p;
data.bio.len_wgt_all.reg_pM=pM(1);

%% %% length-weight regression for all trawls for female
indF=setxor(1:length(ind),indM);
xF=log10(data.bio.len_wgt_all.len(ind(indF)));
yF=log10(data.bio.len_wgt_all.wgt(ind(indF)));
pF=polyfit(xF,yF,1);               % linear regression
w0F=10^pF(2);
wF=w0F*len.^pF(1);
data.bio.len_wgt_all.reg_w0F=w0F;      % w=w0*len^p;
data.bio.len_wgt_all.reg_pF=pF(1);

%% total number of fish individuals at length
data.bio.len_nM=hist(data.bio.len_wgt_all.len(ind(indM)),para.bio.hake_len_bin);
data.bio.len_nF=hist(data.bio.len_wgt_all.len(ind(indF)),para.bio.hake_len_bin);
data.bio.len_nALL=hist(data.bio.len_wgt_all.len,para.bio.hake_len_bin);

%% length-key
data.bio.len_key_M=data.bio.len_nM/sum(data.bio.len_nM);
data.bio.len_key_F=data.bio.len_nF/sum(data.bio.len_nF);
data.bio.len_key_ALL=data.bio.len_nALL/sum(data.bio.len_nALL);


%% weight at length or length-weight-key
len_wgt_M=data.bio.len_wgt_all.reg_w0M*para.bio.hake_len_bin.^data.bio.len_wgt_all.reg_pM;
len_wgt_F=data.bio.len_wgt_all.reg_w0F*para.bio.hake_len_bin.^data.bio.len_wgt_all.reg_pF;
len_wgt_ALL=data.bio.len_wgt_all.reg_w0*para.bio.hake_len_bin.^data.bio.len_wgt_all.reg_p;

%% all valid length-weight measurements (no nan's)
L=data.bio.len_wgt_all.len(ind);
W=data.bio.len_wgt_all.wgt(ind);
Lm=data.bio.len_wgt_all.len(ind(indM));
Wm=data.bio.len_wgt_all.wgt(ind(indM));
Lf=data.bio.len_wgt_all.len(ind(indF));
Wf=data.bio.len_wgt_all.wgt(ind(indF));

%% create length-weight sex structured relations
for i=1:length(para.bio.hake_len_bin)
    %% Male
    if data.bio.len_nM(i) < 5    % bins with less than 5 sample will be replaced by the regression curve
        len_wgt_M(i)=data.bio.len_wgt_all.reg_w0M*para.bio.hake_len_bin(i)^data.bio.len_wgt_all.reg_pM;
    else
        indm=find(para.bio.hake_len_bin(i)-1 < Lm & para.bio.hake_len_bin(i)+1 >= Lm); 
        len_wgt_M(i)=mean(Wm(indm));
    end
    %% Female
    if data.bio.len_nF(i) < 5    % bins with less than 5 sample will be replaced by the regression curve
        len_wgt_F(i)=data.bio.len_wgt_all.reg_w0F*para.bio.hake_len_bin(i)^data.bio.len_wgt_all.reg_pF;
    else
        indf=find(para.bio.hake_len_bin(i)-1 < Lf & para.bio.hake_len_bin(i)+1 >= Lf); 
        len_wgt_F(i)=mean(Wf(indf));
    end
    %% combination of Male & Female
    if data.bio.len_nALL(i) < 5  % bins with less than 5 sample will be replaced by the regression curve
        len_wgt_ALL(i)=data.bio.len_wgt_all.reg_w0*para.bio.hake_len_bin(i)^data.bio.len_wgt_all.reg_p;
    else
        ind1=find(para.bio.hake_len_bin(i)-1 < L & para.bio.hake_len_bin(i)+1 >= L); 
        len_wgt_ALL(i)=mean(W(ind1));
    end
end

%% average length-weight-key per fish over entire survey region (a scalar)
data.bio.ave_len_wgt_M=len_wgt_M*data.bio.len_key_M(:);
data.bio.ave_len_wgt_F=len_wgt_F*data.bio.len_key_F(:);
data.bio.ave_len_wgt_ALL=len_wgt_ALL*data.bio.len_key_ALL(:);

%% length-weight-key per fish over entire survey region (an array)
data.bio.len_wgt_M=len_wgt_M;
data.bio.len_wgt_F=len_wgt_F;
data.bio.len_wgt_ALL=len_wgt_ALL;

return
    