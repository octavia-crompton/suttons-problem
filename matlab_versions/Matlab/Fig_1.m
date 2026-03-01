%---------- Plot wv concentration and vertical flux
figure(1)
clf
%Qd=(Q-Qa)/(Qs-Qa); Normalized concentration
subplot (2,1,1)
pcolor (x,z,Q')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontsize',10)
colorbar
title ('Water vapor concentration (gm m^{-3})')

subplot (2,1,2)
pcolor (x,z,FluxQ')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontSize',10)
colorbar
title ('Water vapor Flux (gm m^{-2} s^{-1})')

