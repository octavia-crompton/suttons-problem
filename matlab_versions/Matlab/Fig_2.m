%---------- Plot wv concentration and vertical flux
figure(2)
clf
QQ=sum(Q'-Qa)*dz;
subplot (2,1,1)
plot (x,QQ,'-','linewidth',3)
hold on
plot (x,max(QQ)*(x/Lx).^0.7,'r--','linewidth',3)
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{\int [Q(x,z)-Q_a]dz}','fontweight','bold','fontsize',10)
title ('Excess water vapor (gm m^{-2})')
legend ('Numerical','Sutton \delta(x)~x^{0.7}','Location','northwest')

subplot (2,1,2)
plot(x,FluxQ(:,1),'-','linewidth',3)
hold on
plot (x,ones(1,nx+1)*FluxQ(nx+1,1),'r--','linewidth',2)

xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{ET} (gm m^{-2} s^{-1})','fontweight','bold','fontSize',10)
title ('Surface Flux (gm m^{-2} s^{-1})')

