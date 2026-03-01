%---------- Plot wv concentration and vertical flux
figure(3)
clf
Qd=(Q-Qa)/(Qs-Qa); %Normalized concentration
subplot(2,1,1)
pcolor (x,z,Qd')
%colormap ('cool')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontsize',10)
colorbar
caxis([0 1])
title ('y=(Q-Q_a)/(Q_s-Q_a)')

subplot(2,1,2)
pcolor (x(1:nx),z(1:nz),Q_sutton(1:nx,1:nz)')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontsize',10)
colorbar
caxis([0 1])
title ('Sutton')

