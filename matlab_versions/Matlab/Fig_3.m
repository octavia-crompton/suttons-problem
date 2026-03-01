%---------- Plot wv concentration and vertical flux
figure(3)
clf
Qd=(Q-Qa)/(Qs-Qa); %Normalized concentration
pcolor (x,z,Qd')
%colormap ('cool')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontsize',10)
colorbar
title ('y=(Q-Q_a)/(Q_s-Q_a)')

