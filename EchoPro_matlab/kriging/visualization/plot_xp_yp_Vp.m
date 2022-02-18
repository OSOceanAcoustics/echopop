function plot_xp_yp_Vp(xp,yp,DispVar)

  
clr_max=64;
cmap=colormap;
markersize=12;
Vmax=max(max(max(DispVar)));
Vmin=min(min(min(DispVar)));
clr_gd=min(max(floor(clr_max*(DispVar-Vmin)/(Vmax-Vmin)),0)+1,clr_max);

figure(2)
for i=1:length(xp)
  if DispVar(i) ~= 0
     plot(xp(i),yp(i),'.','color',cmap(clr_gd(i),:),'markersize',markersize);
  else
     disp([xp(i) yp(i)])
  end
  if i == 1
      hold on
  end
end   
hold off
grid
return