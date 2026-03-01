%-------------------------------------------------------
%Main_Sutton_Advection.m
%
%Authors: Gaby Katul and Octavia Crompton
%
%Date: May 25/2024
%
%Purpose: Models advection of water vapor (wv) from a dry 
%         to a wet surface.
%         
%Reference:Brutsaert, W.B., 1982, Evaporation into 
%          the atmosphere: Theory, history, and 
%          applications, Kluwer Academic Publishers,
%          pp.168-173
%-------------------------------------------------------
clear all
clc
%------- Constants
k=0.4;       % von Karman

%------- flow conditions
ustar=0.25;   % Friction velocity (m/s)
zom=0.1;      % Momentum roughness length (m)

%Water vapor boundary conditions (gm/m3):
Qs=11;       % Surface wv concentration
Qa=8;        % Upwind background atmospheric wv concentration

%------- Domain specification
Lx=200; Hmax=50;
xmin=0;
xmax=Lx;
zmin=zom;
zmax=Hmax;
%-------- Generate grid 
nx=200;
nz=2000;
dx=(xmax-xmin)/nx;
dz=(zmax-zmin)/nz;
z=[zmin+dz:dz:zmax];
x=[xmin:dx:xmax];
%--------- Generate the mean velocity from log-law
U=(ustar/k)*log(z/zom);
%--------- Specify upwind wv concentration (as background)
Qup=ones(1,nz)*Qa;

%--------- Setup coefficients for implicit scheme
A=[]; B=[]; C=[];Q=[];Q1=[]; Q2=[];
A=k*z*ustar;
B=1./U;
[C]=our_central_difference(A,dz);

%--------- Upwind wv concentrations and fluxes
Q1=Qup;
Q(1,1:nz)=Q1;
FluxQ(1,1:nz)=0;

%---------- Begin downwind calculations by marching along x
for i=1:nx
	Q2=[]; Fq=[];
	[Q2,Fq]=integrate_H2O_implicit (nx,nz,dx,dz,A,B,C,Q1,Qs,Qa);
	Q(i+1,1:nz)=Q2; 
	FluxQ(i+1,1:nz)=Fq;
Q1=Q2;
end

%---------- Plot results
Fig_1
Fig_2
Fig_3


