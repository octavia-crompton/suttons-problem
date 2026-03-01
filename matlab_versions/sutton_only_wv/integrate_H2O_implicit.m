%--------------------------------------------------------
%Note: This function computes the coefficients
%      of the second-order system of variable s(t+dt)
%      formed by the (known) profile values at s(t).
%      The system is AA1 d2s/dz2 + AA2 ds/dz + AA3 s=AA4
%      The coefficients in AA4 are used from s(t)
%      The coefficients in AA1 are mainly the diffusion
%      The coefficients in AA2 include advection
%      The coefficients in AA3 include advection + s(t+dt) 
%--------------------------------------------------------
function [Q1, Fq]=integrate_H2O_implicit (n,m,dx,dz,A,B,C,Qup,Qs,Qa)
Q1=[]; Flux=[];
%---- setup the tridiagonal solver for mean H2O concentration
AA1=-A.*B;
AA2=-C.*B;
AA3=1/dx;
AA4=Qup/dx;

upd=(AA1/(dz^2)+AA2/(2*dz));
dia=(-2*AA1/(dz^2)+AA3);
lod=(AA1/(dz^2)-AA2/(2*dz));
co=AA4;

%--- These are set to ensure the bc are state not flux
lod(1)=0;
lod(m)=0;  % change to -1 for zero gradient BC at bottom
dia(1)=1;
dia(m)=1;
upd(1)=0;
upd(m)=0;

% --- surface and upper H2O concentration enforced 
co(1)=Qs;
co(m)=Qa;  % multiply by 0 for zero gradient BC at bottom

aa=lod;
bb=dia;
cc=upd;
dd=co;

%--------Call the tridiagonal solver
Q1=[];Fq=[];
Q1=Thomas(aa,bb,cc,dd);
[dQdz]=our_central_difference(Q1,dz);

Fq=-A.*dQdz;