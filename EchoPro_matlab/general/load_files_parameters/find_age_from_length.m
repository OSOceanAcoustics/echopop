function age = find_age_from_length(data, para, len)
%% 

len_all = data.final.table.trawl(:,9);
age_all = data.final.table.trawl(:,11);
ind = find(isnan(age_all) == 1 | isnan(len_all) == 1);
len_all(ind) =[];
age_all(ind) = [];

if isfield(para.bio, 'length_age_reg_order')
    polyfit_order = para.bio.length_age_reg_order;
else
    polyfit_order = 2;
end

p = polyfit(len_all, age_all, polyfit_order);
age = round(max(1, polyval(p, len)));

% figure(1)
% plot(len_all, age_all, '.', para.bio.hake_len_bin, polyval(p, para.bio.hake_len_bin), '-r')
% grid

return