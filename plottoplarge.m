
load('topologyexp1.mat','z_h','z_u');
x_h = real(z_h);
y_h = imag(z_h);
x_u = real(z_u);
y_u = imag(z_u);
r = 5;
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')

%  for i = 1:num_helpers
%      xxx = x_h(i) - r;
%      yyy = y_h(i) - r;
%      rr = 2*r;
%      rectangle('Position',[xxx,yyy,rr,rr],'Curvature',[1,1],'Linestyle',':','LineWidth',1);
%      hold  all
%  end
 
scatter(x_u,y_u,'b*');
hold all
scatter(x_h,y_h,'*r');
%scatter(0,0,'^g');
set(gca,'XTick',-22.5:5:17.5);
set(gca, 'YTick',-22.5:5:17.5);
grid(gca)