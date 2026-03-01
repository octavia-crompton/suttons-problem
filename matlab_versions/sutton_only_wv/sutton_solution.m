function [Q_sutton]=sutton_solution (nx,nz,x,z,ustar,aa1,nn1)
%--------- Compare to Higgins et al. (2013)
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
