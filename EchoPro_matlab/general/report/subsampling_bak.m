function out=subsampling(dat,int_for_zeros)

xt=unique(dat(:,1));

out=[];
for i=1:length(xt)
  %% find all samples on the same transect
    ind=find(dat(:,1) == xt(i));
    vlm=nanmean(dat(ind,3:4)')';
    na=dat(ind,12);
    lat=dat(ind,5);
    lon=dat(ind,6);
    w=dat(ind,21);
 %% merge data samples age-2+ and hake mix, and possibly age-1
    [vlm_sort, indx_sort]=sort(vlm);
    lat_sort=lat(indx_sort);
    lon_sort=lon(indx_sort);
    if i == 67
        disp(i)
    end
    w_sort=w(indx_sort);
    na_sort=na(indx_sort);
    ind_int=find(diff(vlm_sort) < 0.1);
    if ~isempty(ind_int)
        na_sort(ind_int)=na_sort(ind_int)+na_sort(ind_int+1);
        w_sort(ind_int)=w_sort(ind_int)+w_sort(ind_int+1);
        na_sort(ind_int+1)=[];
        w_sort(ind_int+1)=[];
        lat_sort(ind_int+1)=[];
        lon_sort(ind_int+1)=[];
    end
    %% find intervels with zeros & non_zeros
    n=length(lat_sort);
    ind_zero=find(na_sort(ind:n) == 0);
    ind_non_zero=find(na_sort(ind:n) ~= 0);
    dind_zero=find(diff(ind_zero) > 1);
    if isempty(dind_zero)
        ind=1:int_for_zeros:n;
        out=[out; lat_sort(ind) lon_sort(ind) w_sort(ind)];
    else
        dind_non_zero=find(diff(ind_non_zero) > 1);
        for j=1:length(dind_non_zero)
            if j == 1
                in_zeros=1:(ind_non_zero(1)-1);
                ind=in_zeros(1):int_for_zeros:in_zeros(end);
                out=[lat_sort(ind) lon_sort(ind) w_sort(ind)];
                in_non_zero=ind_non_zero(1):ind_zero(dind_zero(i))-1;
                out=[out; lat_sort(in_non_zero) lon_sort(in_non_zero) w_sort(in_non_zero)];
            elseif j == length(dind_non_zero)
                in_zeros=ind_non_zero(end)+1:n;
                ind=in_zeros(1):int_for_zeros:in_zeros(end);
                out=[out; lat_sort(ind) lon_sort(ind) w_sort(ind)];
            else
                in_zeros=ind_zero(dind_zero(j-1)):(ind_non_zero(dind_non_zero(j))-1);
                ind=in_zeros(1):int_for_zeros:in_zeros(end);
                out=[out; lat_sort(ind) lon_sort(ind) w_sort(ind)];
                in_non_zero=ind(end)+1:ind_zero(dind_zero(j))-1;
                out=[out; lat_sort(in_non_zero) lon_sort(in_non_zero) w_sort(in_non_zero)];
            end
        end
    end
end

return