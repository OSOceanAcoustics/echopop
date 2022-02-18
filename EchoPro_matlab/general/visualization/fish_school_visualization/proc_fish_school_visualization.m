function proc_pb_fish_school_visualization
% visualize the 2D fish school results
%
%% Written by Dezhang Chu    Written by Dezhang Chu, NOAA Fisheries, NWFSC
%% Last Modification:       4/4/2013
global para hdl data

figure(13)
if ~isfield(data,'visual') | ~isfield(data.visual,'school')
    data.visual.school.func=data.final.table.school(:,8)*1e-6;
    data.visual.school.func_str='Area (km^2)';
    data.visual.school.var=data.final.table.school(:,1);
    data.visual.school.var_str='Transect';
else
    if ~isfield(data.visual.school,'func')
        data.visual.school.func=data.final.table.school(:,8)*1e-6;
        data.visual.school.func_str='Area (km^2)';
    elseif ~isfield(data.visual.school,'var')
        data.visual.school.var=data.final.table.school(:,1);
        data.visual.school.var_str='Transect';
    end
end

if get(hdl.vis.school.radio_school_hist,'value') == 1
    %% plot histogram
    [y,x]=hist(data.visual.school.func);
    bar(x,y);
    data.visual.school.var1=x;
    data.visual.school.func1=y;
    data.visual.school.var1_str=data.visual.school.func_str;
    data.visual.school.func1_str='Counts';
    xlabel(data.visual.school.func_str,'fontsize',20,'fontweight','bold');
    ylabel(data.visual.school.var_str,'fontsize',20,'fontweight','bold');
    return
end

if strcmp(data.visual.school.var_str,'Transect')
    var=unique(data.visual.school.var);
    for i=1:length(var)
        ind=find(data.visual.school.var == var(i));
        func(i)=nanmean(data.visual.school.func(ind));
    end
    [sort_var,sort_ind]=sort(var,'ascend');
    plot(sort_var,func(sort_ind),'.-')
    data.visual.school.var1=sort_var;
    data.visual.school.func1=func(sort_ind);
    data.visual.school.var1_str=data.visual.school.func_str;
    data.visual.school.func1_str='Counts';
elseif strcmp(data.visual.school.var_str,'Region Sequence No.')
    var=1:length(data.visual.school.var);
    plot(var,data.visual.school.func,'.-')
    data.visual.school.var1=var;
    data.visual.school.func1=data.visual.school.func;
else
    [sort_var,sort_ind]=sort(data.visual.school.var,'ascend');
    plot(sort_var,data.visual.school.func(sort_ind),'.-')
    data.visual.school.var1=sort_var;
    data.visual.school.func1=data.visual.school.func(sort_ind);
end
ylabel(data.visual.school.func_str,'fontsize',20,'fontweight','bold');
xlabel(data.visual.school.var_str,'fontsize',20,'fontweight','bold');
data.visual.school.var1_str=data.visual.school.var_str;
data.visual.school.func1_str=data.visual.school.func_str;

return