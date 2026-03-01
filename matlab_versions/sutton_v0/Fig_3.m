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

%--------- Compare to Higgins et al. (2013)
aa1=0.8;nn1=0.25;
mm1=nn1/(2-nn1);
yy2=mm1/(2+mm1);
yy1a=(aa1*mm1)/(ustar*(2+mm1-nn1)^2);
yy1b=yy1a*z.^(2+mm1-nn1);
Q_sutton=[];
Q_sutton(1,1:nz)=0*z;
for i=2:nx+1
    yy1=yy1b/x(i);
    Q_sutton(i,:)=gammainc(yy1,yy2,'upper');
end

subplot(2,1,2)
pcolor (x(1:nx),z(1:nz),Q_sutton(1:nx,1:nz)')
shading ('interp')
xlabel ('\it{x} (m)','fontweight','bold','fontsize',10)
ylabel ('\it{z} (m)','fontweight','bold','fontsize',10)
colorbar
caxis([0 1])
title ('Sutton')

