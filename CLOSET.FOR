      Program Closure
      
C     This program is designed for solving the Internal Boundary Layer
C       by a Second-Order Closure model by Rao et al. (1974).
C     Cheng-I Hsieh, Nov. 5, 1996 (draft);
C     1/14/97 (Revised which Wi is added); 
C     4/8/97  (Revised with boundary condition determination)      
C     -------------Definition of Parameters----------------------------
C     x=horizontal distance (m)
C     z=vertical distance (m)  
C     dx, dz=step size of x and z coordinate, respectively 
C     dsx, dsz=coefficient for determing dx and dz, respectively
C     zo=surface roughness (m).
C     us=u*=friction velocity (m/s).  
C     ts=t*=temperature scale (K).  
C     qs=q*=specific humidity scale (g/g).                          
C     U=<U>=mean velocity of U, and so on. <> denotes averaging .
C     uu=<uu>=variance of u, and so on. 
C     tw=<tw>=covariance of t and w, and so on. 
C     e=dissipation rate of turbulent kinetic energy (m^2/s^3).
C     e2=uu+vv+ww  
C     nx, nz=size of dimension    
C     Constants kuu,kvv,kww...... are calculated by Hsieh, 11/13,96
C     Constants au,av,aw...... are based on Kader and Yaglom (1990)
C     ----------------------------------------------------------------- 
      REAL k,ke,kuu,kvv,kww,kuw,ktu,ktw,kqu,kqw,ktt,kqq,kqt,Le                            
      Parameter (g=9.8, Cp=1005.0, Le=2440000.0,RH=0.6) 
C     *****g(m/s^2), Cp(J/kg/K), Le(J/kg), RH=relative humidity********                                      
      Parameter (nx=1000, nz=55, dsx=0.04, dsz=0.35) 
C     *****Values for dsx and dsz should be specified to get stable
C          and converged solutions*************************************   
      Parameter (zo1=0.1,Uo1=0.0,Wo1=0.0,To1=298.0,Qo1=6.0)
      Parameter (uso1=0.4,two1=0.0,qwo1=0.1)  
      Parameter (zo2=0.1)
C      Parameter (zo2=0.1,Uo2=0.0,Wo2=0.0,To2=291.0,Qo2=8.9)      
C     *****The unit for Q is (g/kg), the value should be consistent
C          with the value used for <qw>********************************
C     *****To2 and Qo2 are determined from dT/dQ=-Le/Cp**************** 
C     *****The following constants are calculated by Cheng-I***********
c      Parameter (au=7.29,av=6.25,aw=1.56,at=8.41,aq=8.41,am=1.0)
c      Parameter (atu=4.0,aqu=4.0,aqt=4.0, ae=1.0)
c      Parameter (k=0.4,ke=0.017,kuu=8.93,kvv=-8.28,kww=2.90,kuw=23.60)
c      Parameter (ktu=7.55,ktw=24.25,kqu=7.55,kqw=24.25)
c      Parameter (ktt=1.80,kqq=1.80,kqt=1.80)                           
C     *****The following constants are from Rao et al., 1974***********
      Parameter (au=4.0,av=1.75,aw=1.75,at=4.0,aq=4.0)
      Parameter (atu=3.0,aqu=3.0,aqt=3.0)     
      Parameter (k=0.4,ke=0.0487,kuu=6.67,kvv=6.67,kww=6.67,kuw=13.125)
      Parameter (ktu=4.35,ktw=9.87,kqu=4.35,kqw=9.87)
      Parameter (ktt=1.39,kqq=1.39,kqt=1.39)
      Real U(nx,nz),W(nx,nz),T(nx,nz),Q(nx,nz),uu(nx,nz),vv(nx,nz)
      Real ww(nx,nz),uw(nx,nz),tu(nx,nz),tw(nx,nz),qu(nx,nz),qw(nx,nz)
      Real tt(nx,nz),qq(nx,nz),qt(nx,nz),e2(nx,nz),e(nx,nz),x(nx),z(nz)
      Real us(nx,nz),ts(nx,nz),qs(nx,nz),dz(nz),dx(nx)
C     -----------------------------------------------------------------
C     --------------Construct the Log(z) and x Coordinates-------------
      IF (zo1 .GE. zo2) THEN
        z(1)=zo1
       ELSE                 
        z(1)=zo2
      ENDIF  
c      DO n=2, nz 
      nzu=14
      DO n=2, nzu
       z(n)=Exp(Log(z(1))+(n-1)*dsz)
      dz(n-1)=z(n)-z(n-1)
      ENDDO 
C     *****After n=nzu, dz is uniform***********************************          
      DO n=nzu+1, nz
       z(n)=z(nzu)+(n-nzu)*dz(nzu-1)
      dz(n-1)=dz(nzu-1)
      ENDDO 
      
      x(1)=0.0
      DO j=2, nx
       x(j)=x(1)+(j-1)*dsx
      dx(j-1)=x(j)-x(j-1)
      ENDDO
C     -----------------------------------------------------------------
      To=To1
C     --------------Upstream Initial Conditions------------------------ 
      WRITE (*,*) 'Calculating the Initial Conditions'
      us(1,1)=uso1        
      tw(1,1)=two1
      qw(1,1)=qwo1
      ts(1,1)=-1.*tw(1,1)/us(1,1)
      qs(1,1)=-1.*qw(1,1)/us(1,1)
       U(1,1)=us(1,1)/k*Log(z(1)/zo1)
       W(1,1)=Wo1 
       T(1,1)=To1+0.74*ts(1,1)/k*Log(z(1)/zo1)
       Q(1,1)=Qo1+0.74*qs(1,1)/k*Log(z(1)/zo1)
      uu(1,1)=au*(us(1,1))**2.
      vv(1,1)=av*(us(1,1))**2.
      ww(1,1)=aw*(us(1,1))**2.
      uw(1,1)=-1.*us(1,1)**2. 
      tt(1,1)=at*(ts(1,1))**2.
      qq(1,1)=aq*(qs(1,1))**2.
      tu(1,1)=atu*ts(1,1)*us(1,1)
      qu(1,1)=aqu*qs(1,1)*us(1,1)
      qt(1,1)=aqt*qs(1,1)*ts(1,1) 
      e2(1,1)=uu(1,1)+vv(1,1)+ww(1,1) 
       e(1,1)=ke*(e2(1,1)**(3./2.))/(k*z(1))
C     *****Set the initial profiles by MOST with neutral condition*****            
      DO n=2, nz
       us(1,n)=us(1,1)
       ts(1,n)=ts(1,1)
       qs(1,n)=qs(1,1)
        U(1,n)=us(1,n)/k*Log(z(n)/zo1)
        W(1,n)=W(1,1)
        T(1,n)=To1+0.74*ts(1,n)/k*Log(z(n)/zo1)
        Q(1,n)=Qo1+0.74*qs(1,n)/k*Log(z(n)/zo1)
       uu(1,n)=uu(1,1)
       vv(1,n)=vv(1,1)
       ww(1,n)=ww(1,1)
       uw(1,n)=uw(1,1)
       tt(1,n)=tt(1,1)
       qq(1,n)=qq(1,1)
       qt(1,n)=qt(1,1)
       tu(1,n)=tu(1,1)
       tw(1,n)=tw(1,1)
       qu(1,n)=qu(1,1)
       qw(1,n)=qw(1,1)
       e2(1,n)=uu(1,n)+vv(1,n)+ww(1,n)  
       eL=k*z(n)                           
      IF (z(n) .GT. 10.0) THEN
      eL=k*z(n)/(1+k*(z(n)-10.0)/40.0)
      ENDIF
        e(1,n)=ke*(e2(1,n)**(3./2.))/(eL)
C     *****Assign the boundary conditions for x < 0********************       
      DO j=2, nx 
       us(j,n)=us(1,n)
       ts(j,n)=ts(1,n)
       qs(j,n)=qs(1,n)
        U(j,n)=U(1,n)
        W(j,n)=W(1,n)
        T(j,n)=T(1,n)
        Q(j,n)=Q(1,n)
       uu(j,n)=uu(1,n)
       vv(j,n)=vv(1,n)
       ww(j,n)=ww(1,n)
       uw(j,n)=uw(1,n)
       tt(j,n)=tt(1,n)
       qq(j,n)=qq(1,n)
       qt(j,n)=qt(1,n)
       tu(j,n)=tu(1,n)
       tw(j,n)=tw(1,n)
       qu(j,n)=qu(1,n)
       qw(j,n)=qw(1,n)
       e2(j,n)=e2(1,n)
        e(j,n)=e(1,n)
      ENDDO          
      ENDDO     
C     *****Calculate the Equilibrium Initial Profiles******************
      DO n=5, nz-1
      duu=1.0
      de2=1.0     
      DO j=1, nx-1 
      write (*,*) 'n,j',n,j
C     (1) U-->uw
      uw(j+1,n)=uw(1,1)
C     (2) T-->tw
      tw(j+1,n)=tw(1,1)
C     (3) Q-->qw
      qw(j+1,n)=qw(1,1)
C     (4) W
      W(j+1,n)=W(1,1)
C     (8) uw-->U
c      U(j,n+1)=((-kuw*e(j,n)/e2(j,n)*uw(j,n))+(g/To*tu(j,n)))*
c     + (z(n+1)-z(n-1))/ww(j,n)+U(j,n-1) 
      U(j,n)=((-kuw*e(j,n)/e2(j,n)*uw(j,n))+(g/To*tu(j,n)))*
     + (z(n)-z(n-1))/ww(j,n)+U(j,n-1)
      Ue=U(j,n)
      write (*,*) 'U(j,n),U(1,n),Uj,n-1',U(j,n),U(1,n),U(j,n-1)
C     (5)-uu                  
      uu(j+1,n)=((-kuu*e(j,n)/e2(j,n)*(uu(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))-2*(0*uu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+
     + uw(j,n)*(U(j,n)-U(j,n-1))/(z(n)-z(n-1)))-(W(j,n)*
     + (uu(j,n)-uu(j,n-1))/(z(n)-z(n-1))))*(dx(j))/U(j,n)+uu(j,n) 
      uue=uu(j,n)
C     (6)-vv
      vv(j+1,n)=((-kvv*e(j,n)/e2(j,n)*(vv(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))-(W(j,n)*(vv(j,n)-vv(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+vv(j,n)            
      vve=vv(j,n)
C     (7)-ww
      ww(j+1,n)=((-kww*e(j,n)/e2(j,n)*(ww(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))+2*(g/To*tw(j,n))-2*(ww(j,n)*
c     + (-1)*(U(j+1,n)-U(j,n))/(dx(j))+0*uw(j,n)*(W(j+1,n)-W(j,n))/ 

c     + (-1)*(U(j+1,n)-U(j,n))/(dx(j))+uw(j,n)*(W(j+1,n)-W(j,n))/

     + (W(j,n)-W(j,n-1))/(z(n)-z(n-1))+0*uw(j,n)*(W(j+1,n)-W(j,n))/ 
     + (dx(j)))-(W(j,n)*(ww(j,n)-ww(j,n-1))/(z(n)-z(n-1))))*(dx(j))/
     + U(j,n)+ww(j,n) 
      wwe=ww(j,n)

      write (*,*) 'uu(j,n),uu(1,n)',uu(j,n),uu(1,n)     
      write (*,*) 'ww(j,n),ww(1,n)',ww(j,n),ww(1,n)  

C     (10) tw-->T
c      T(j,n+1)=((-ktw*e(j,n)/e2(j,n)*tw(j,n))+(g/To*tt(j,n)))*
c     + (z(n+1)-z(n-1))/ww(j,n)+T(j,n-1)   
      T(j,n)=((-ktw*e(j,n)/e2(j,n)*tw(j,n))+(g/To*tt(j,n)))*
     + (z(n)-z(n-1))/ww(j,n)+T(j,n-1) 
      Te=T(j,n)
C     (9)-tu      
      tu(j+1,n)=((-ktu*e(j,n)/e2(j,n)*tu(j,n))-(0*uu(j,n)*(T(j+1,n)-
     + T(j,n))/(dx(j))+uw(j,n)*(T(j,n)-T(j,n-1))/(z(n)-z(n-1)))-
     + (0*tu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+tw(j,n)*(U(j,n)-U(j,n-1))/
     + (z(n)-z(n-1)))-(W(j,n)*(tu(j,n)-tu(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+tu(j,n)  
      tue=tu(j,n)
C     (12) qw-->Q
c      Q(j,n+1)=((-kqw*e(j,n)/e2(j,n)*qw(j,n))+(g/To*qt(j,n)))*
c     + (z(n+1)-z(n-1))/ww(j,n)+Q(j,n-1)    
      Q(j,n)=((-kqw*e(j,n)/e2(j,n)*qw(j,n))+(g/To*qt(j,n)))*
     + (z(n)-z(n-1))/ww(j,n)+Q(j,n-1)   
      Qe=Q(j,n)
C     (11)-qu      
      qu(j+1,n)=((-kqu*e(j,n)/e2(j,n)*qu(j,n))-(0*uu(j,n)*(Q(j+1,n)-
     + Q(j,n))/(dx(j))+uw(j,n)*(Q(j,n)-Q(j,n-1))/(z(n)-z(n-1)))-
     + (0*qu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+qw(j,n)*(U(j,n)-U(j,n-1))/
     + (z(n)-z(n-1)))-(W(j,n)*(qu(j,n)-qu(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+qu(j,n) 
      que=qu(j,n)
C     (13)-tt
      tt(j+1,n)=((-2*ktt*e(j,n)/e2(j,n)*tt(j,n))-2*(0*tu(j,n)*
     + (T(j+1,n)-T(j,n))/(dx(j))+tw(j,n)*(T(j,n)-T(j,n-1))/
     + (z(n)-z(n-1)))-(W(j,n)*(tt(j,n)-tt(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+tt(j,n) 
      tte=tt(j,n)
C     (14)-qq
      qq(j+1,n)=((-2*kqq*e(j,n)/e2(j,n)*qq(j,n))-2*(0*qu(j,n)*
     + (Q(j+1,n)-Q(j,n))/(dx(j))+qw(j,n)*(Q(j,n)-Q(j,n-1))/
     + (z(n)-z(n-1)))-(W(j,n)*(qq(j,n)-qq(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+qq(j,n) 
      qqe=qq(j,n)
C     (15)-qt
      qt(j+1,n)=((-2*kqt*e(j,n)/e2(j,n)*qt(j,n))-(0*tu(j,n)*(Q(j+1,n)-
     + Q(j,n))/(dx(j))+tw(j,n)*(Q(j,n)-Q(j,n-1))/(z(n)-z(n-1)))-
     + (0*qu(j,n)*(T(j+1,n)-T(j,n))/(dx(j))+qw(j,n)*(T(j,n)-T(j,n-1))/
     + (z(n)-z(n-1)))-(W(j,n)*(qt(j,n)-qt(j,n-1))/(z(n)-z(n-1))))*
     + (dx(j))/U(j,n)+qt(j,n)        
      qte=qt(j,n)
C     e & e2
      e2(j+1,n)=uu(j+1,n)+vv(j+1,n)+ww(j+1,n)  
      e2e=e2(j,n)
      eL=k*z(n)                           
      IF (z(n) .GT. 10.0) THEN
      eL=k*z(n)/(1+k*(z(n)-10.0)/40.0)
      ENDIF
      e(j+1,n)=ke*(e2(j+1,n)**(3./2.))/(eL)  
      ee=e(j,n)
                                 
      duu=ABS(uu(j+1,n)-uu(j,n))                           
      de2=ABS(e2(j+1,n)-e2(j,n))
      IF (duu .LT. 0.001) THEN
      je=j
      GO TO 5  
      ENDIF
      ENDDO
5     DO L=1, nx   
c        write (*,*) 'jf+1',jf+1
        U(L,n)=Ue
c        W(L+1,n)=W(je,n)
        T(L,n)=Te
        Q(L,n)=Qe
       uu(L+1,n)=uue
       vv(L+1,n)=vve
       ww(L+1,n)=wwe
c       uw(L+1,n)=uw(je,n)
       tt(L+1,n)=tte
       qq(L+1,n)=qqe
       qt(L+1,n)=qte
       tu(L+1,n)=tue
c       tw(L+1,n)=tw(je,n)
       qu(L+1,n)=que
c       qw(L+1,n)=qw(je,n)
       e2(L+1,n)=e2e
        e(L+1,n)=ee  
c        write (*,*) 'uu(jf+1,n)',uu(jf+1,n)  
       
      ENDDO  
      write (*,*) 'je', je
      pause              
      ENDDO
C     -----------------------------------------------------------------      
C     --------------Calculating the Profiles as f(x, z)----------------
C     -----Finite Difference Method - Forward x, Centred z Scheme------ 
      WRITE (*,*) '=====> Calculating the Profiles'
      DO j=1, nx-1  
      write (*,*) 'x(j)',j,x(j)
      DO n=2, nz-1
C     (4)-W
c      W(j+1,n)=((-1*(ww(j,n+1)-ww(j,n-1))/(z(n+1)-z(n-1)))-
c     + (W(j,n)*(W(j,n+1)-W(j,n-1))/(z(n+1)-z(n-1))))*dx(j)/U(j,n)+W(j,n)
C     This Eq. for Wi is an extra to Rao et. al.'s (1974) model.
C     (1)-U 
      U(j+1,n)=((-1*(uw(j,n+1)-uw(j,n-1))/(z(n+1)-z(n-1)))-
     + (W(j,n)*(U(j,n+1)-U(j,n-1))/(z(n+1)-z(n-1))))*dx(j)/U(j,n)+U(j,n)
c      write (*,*) 'W(j+1,n)-W(j,n)',W(j+1,n)-W(j,n),ww(j,3)
c      write (*,*) 'U(j,n),U(j+1,n)',j,n,U(j,n),U(j+1,n) 
c      pause  
C     (4)-W  
c      W(j,2)=0.0
c      W(j,n+1)=(-1.*(U(j+1,n)-U(j,n))/(dx(j)))*(z(n+1)-z(n-1))+W(j,n-1) 
      W(j+1,2)=0.0
      W(j+1,n+1)=(-1.*(U(j+1,n)-U(j,n))/(dx(j)))*(z(n+1)-z(n-1))+
     + W(j+1,n-1)
c      write (*,*) 'W(j,n),W(j,n+1)',j,n,W(j,n),W(j,n+1)  
c      pause
C     (5)-uu
      uu(j+1,n)=((-kuu*e(j,n)/e2(j,n)*(uu(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))-2*(uu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+
     + uw(j,n)*(U(j,n+1)-U(j,n-1))/(z(n+1)-z(n-1)))-(W(j,n)*
     + (uu(j,n+1)-uu(j,n-1))/(z(n+1)-z(n-1))))*(dx(j))/U(j,n)+uu(j,n)
C     (6)-vv
      vv(j+1,n)=((-kvv*e(j,n)/e2(j,n)*(vv(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))-(W(j,n)*(vv(j,n+1)-vv(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+vv(j,n)
C     (7)-ww
      ww(j+1,n)=((-kww*e(j,n)/e2(j,n)*(ww(j,n)-e2(j,n)/3))-
     + (2./3.*e(j,n))+2*(g/To*tw(j,n))-2*(ww(j,n)*
c     + (-1)*(U(j+1,n)-U(j,n))/(dx(j))+0*uw(j,n)*(W(j+1,n)-W(j,n))/ 

c     + (-1)*(U(j+1,n)-U(j,n))/(dx(j))+uw(j,n)*(W(j+1,n)-W(j,n))/

     + (W(j,n+1)-W(j,n-1))/(z(n+1)-z(n-1))+0*uw(j,n)*(W(j+1,n)-W(j,n))/ 
     + (dx(j)))-(W(j,n)*(ww(j,n+1)-ww(j,n-1))/(z(n+1)-z(n-1))))*(dx(j))/
     + U(j,n)+ww(j,n)

c      write (*,*) 'uu(j,n),uu(j+1,n)',uu(j,n),uu(j+1,n)     
c      write (*,*) 'ww(j,n),ww(j,n+1)',ww(j,n),ww(j,n+1)  
c      pause
C     (8)-uw
      uw(j+1,n)=((-kuw*e(j,n)/e2(j,n)*uw(j,n))+(g/To*tu(j,n))-
     + (ww(j,n)*(U(j,n+1)-U(j,n-1))/(z(n+1)-z(n-1))+0*uu(j,n)*   
     
c     + (ww(j,n)*(U(j,n+1)-U(j,n-1))/(z(n+1)-z(n-1))+uu(j,n)*   
  
     + (W(j+1,n)-W(j,n))/(dx(j)))-(W(j,n)*(uw(j,n+1)-uw(j,n-1))/
     + (z(n+1)-z(n-1))))*(dx(j))/U(j,n)+uw(j,n)
C     (2)-T
      T(j+1,n)=((-1*(tw(j,n+1)-tw(j,n-1))/(z(n+1)-z(n-1)))-
     + (W(j,n)*(T(j,n+1)-T(j,n-1))/(z(n+1)-z(n-1))))*dx(j)/U(j,n)+T(j,n) 
C     (9)-tu      
      tu(j+1,n)=((-ktu*e(j,n)/e2(j,n)*tu(j,n))-(uu(j,n)*(T(j+1,n)-
     + T(j,n))/(dx(j))+uw(j,n)*(T(j,n+1)-T(j,n-1))/(z(n+1)-z(n-1)))-
     + (tu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+tw(j,n)*(U(j,n+1)-U(j,n-1))/
     + (z(n+1)-z(n-1)))-(W(j,n)*(tu(j,n+1)-tu(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+tu(j,n) 
C     (10)-tw
      tw(j+1,n)=((-ktw*e(j,n)/e2(j,n)*tw(j,n))+(g/To*tt(j,n))-
     + (uw(j,n)*(T(j+1,n)-T(j,n))/(dx(j))+ww(j,n)*(T(j,n+1)-T(j,n-1))/
c     + (z(n+1)-z(n-1)))-(tw(j,n)*(-1)*(U(j+1,n)-U(j,n))/(dx(j))+0*tu(j,n)*  
     
c     + (z(n+1)-z(n-1)))-(tw(j,n)*(-1)*(U(j+1,n)-U(j,n))/(dx(j))+tu(j,n)*  

     + (z(n+1)-z(n-1)))-(tw(j,n)*(W(j,n+1)-W(j,n-1))/(z(n+1)-z(n-1))+
     + 0*tu(j,n)* 
     + (W(j+1,n)-W(j,n))/(dx(j)))-(W(j,n)*(tw(j,n+1)-tw(j,n-1))/
     + (z(n+1)-z(n-1))))*(dx(j))/U(j,n)+tw(j,n)
C     (3)-Q
      Q(j+1,n)=((-1*(qw(j,n+1)-qw(j,n-1))/(z(n+1)-z(n-1)))-
     + (W(j,n)*(Q(j,n+1)-Q(j,n-1))/(z(n+1)-z(n-1))))*dx(j)/U(j,n)+Q(j,n) 
C     (11)-qu      
      qu(j+1,n)=((-kqu*e(j,n)/e2(j,n)*qu(j,n))-(uu(j,n)*(Q(j+1,n)-
     + Q(j,n))/(dx(j))+uw(j,n)*(Q(j,n+1)-Q(j,n-1))/(z(n+1)-z(n-1)))-
     + (qu(j,n)*(U(j+1,n)-U(j,n))/(dx(j))+qw(j,n)*(U(j,n+1)-U(j,n-1))/
     + (z(n+1)-z(n-1)))-(W(j,n)*(qu(j,n+1)-qu(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+qu(j,n) 
C     (12)-qw
      qw(j+1,n)=((-kqw*e(j,n)/e2(j,n)*qw(j,n))+(g/To*qt(j,n))-
     + (uw(j,n)*(Q(j+1,n)-Q(j,n))/(dx(j))+ww(j,n)*(Q(j,n+1)-Q(j,n-1))/
c     + (z(n+1)-z(n-1)))-(qw(j,n)*(-1)*(U(j+1,n)-U(j,n))/(dx(j))+0*qu(j,n)*  
     
c     + (z(n+1)-z(n-1)))-(qw(j,n)*(-1)*(U(j+1,n)-U(j,n))/(dx(j))+qu(j,n)*   

     + (z(n+1)-z(n-1)))-(qw(j,n)*(W(j,n+1)-W(j,n-1))/(z(n+1)-z(n-1))+
     + 0*qu(j,n)* 
     + (W(j+1,n)-W(j,n))/(dx(j)))-(W(j,n)*(qw(j,n+1)-qw(j,n-1))/
     + (z(n+1)-z(n-1))))*(dx(j))/U(j,n)+qw(j,n)  
     
c      write (*,*) 'qw(j+1,n-1),qw(j+1,n)',qw(j+1,n-1),qw(j+1,n) 
c      write (*,*) '-kqw*e(j,n)/e2(j,n)*qw',-kqw*e(j,n)/e2(j,n)*qw(j,n)
c      write (*,*) 'g/To*qt(j,n)',g/To*qt(j,n)
c      write (*,*) 'uw(j,n)*dQ/dx',uw(j,n)*(Q(j+1,n)-Q(j,n))/(dx(j))
c      write (*,*) 'ww*dQ/(2*dz)',ww(j,n)*(Q(j,n+1)-Q(j,n-1))/(2*dz(n))
c      write (*,*) 'qw*dW/(2*dz)',qw(j,n)*(W(j,n+1)-W(j,n-1))/(2*dz(n))
c      pause
C     (13)-tt
      tt(j+1,n)=((-2*ktt*e(j,n)/e2(j,n)*tt(j,n))-2*(tu(j,n)*
     + (T(j+1,n)-T(j,n))/(dx(j))+tw(j,n)*(T(j,n+1)-T(j,n-1))/
     + (z(n+1)-z(n-1)))-(W(j,n)*(tt(j,n+1)-tt(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+tt(j,n)
C     (14)-qq
      qq(j+1,n)=((-2*kqq*e(j,n)/e2(j,n)*qq(j,n))-2*(qu(j,n)*
     + (Q(j+1,n)-Q(j,n))/(dx(j))+qw(j,n)*(Q(j,n+1)-Q(j,n-1))/
     + (z(n+1)-z(n-1)))-(W(j,n)*(qq(j,n+1)-qq(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+qq(j,n)
C     (15)-qt
      qt(j+1,n)=((-2*kqt*e(j,n)/e2(j,n)*qt(j,n))-(tu(j,n)*(Q(j+1,n)-
     + Q(j,n))/(dx(j))+tw(j,n)*(Q(j,n+1)-Q(j,n-1))/(z(n+1)-z(n-1)))-
     + (qu(j,n)*(T(j+1,n)-T(j,n))/(dx(j))+qw(j,n)*(T(j,n+1)-T(j,n-1))/
     + (z(n+1)-z(n-1)))-(W(j,n)*(qt(j,n+1)-qt(j,n-1))/(z(n+1)-z(n-1))))*
     + (dx(j))/U(j,n)+qt(j,n) 
C     e & e2
      e2(j+1,n)=uu(j+1,n)+vv(j+1,n)+ww(j+1,n)  
      eL=k*z(n)                           
      IF (z(n) .GT. 10.0) THEN
      eL=k*z(n)/(1+k*(z(n)-10.0)/40.0)
      ENDIF
      e(j+1,n)=ke*(e2(j+1,n)**(3./2.))/(eL) 
C     (4)-W
c      W(j,n)=(-1.*(U(j+1,n-1)-U(j,n-1))/(dx(j)))*(z(n+1)-z(n-1))+W(j,n-2)
C      W(j,n) is the best one so far 1/21/97 10:35 am  
      ENDDO    
C     *****Assign the Boundary Value for nx and nz*********************
        U(j+1,nz)=U(j+1,nz-1)
c        W(j,nz)=W(j,nz-1) 
c        W(j+1,nz)=W(j+1,nz-1)
        T(j+1,nz)=T(j+1,nz-1)
        Q(j+1,nz)=Q(j+1,nz-1)
       uu(j+1,nz)=uu(j+1,nz-1)
       vv(j+1,nz)=vv(j+1,nz-1)
       ww(j+1,nz)=ww(j+1,nz-1)
       uw(j+1,nz)=uw(j+1,nz-1)
       tt(j+1,nz)=tt(j+1,nz-1)
       qq(j+1,nz)=qq(j+1,nz-1)
       qt(j+1,nz)=qt(j+1,nz-1)
       tu(j+1,nz)=tu(j+1,nz-1)
       tw(j+1,nz)=tw(j+1,nz-1)
       qu(j+1,nz)=qu(j+1,nz-1)
       qw(j+1,nz)=qw(j+1,nz-1)
       e2(j+1,nz)=e2(j+1,nz-1)
        e(j+1,nz)=e(j+1,nz-1)      
C     -----Calculate the downstream boundary conditions (j,1), j>=2----
C     *****Calculate us, T, Q by Rao et al., 1974 and Stull, 1988******
      us(j+1,1)=U(j+1,2)*k/Log(z(2)/zo2)
c      us(j+1,1)=k*(z(1))*((U(j+1,2)-0.0)/(z(2)-zo2)) 
      Pa=100.0 
      Q(j+1,2)=Q(j+1,2)/1000. 
      qw(1,1)=qwo1/1000.
C     *****Transform unit for Q and qw*********************************                       
      dT=2.0    
      dsT=0.001
C     *****Make a guess for dT*****************************************      
      T(j+1,1)=T(j+1,2)-dT 
      Esat=0.6112*Exp(17.67*(T(j+1,1)-273.16)/(T(j+1,1)-29.66))
      Qsat=0.622*Esat/Pa
      Q(j+1,1)=RH*Qsat
      Diff0=ABS(((T(j+1,2)-T(j+1,1))/dz(1))+Le/Cp*((Q(j+1,2)-Q(j+1,1))/
     +    dz(1))-0.74/(k*z(1)*us(j+1,1))*(-1.)*(tw(1,1)+Le/Cp*qw(1,1))) 
      
      dT=dT+dsT     
      T(j+1,1)=T(j+1,2)-dT 
      Esat=0.6112*Exp(17.67*(T(j+1,1)-273.16)/(T(j+1,1)-29.66))
      Qsat=0.622*Esat/Pa
      Q(j+1,1)=RH*Qsat
      Diff=ABS(((T(j+1,2)-T(j+1,1))/dz(1))+Le/Cp*((Q(j+1,2)-Q(j+1,1))/
     +    dz(1))-0.74/(k*z(1)*us(j+1,1))*(-1.)*(tw(1,1)+Le/Cp*qw(1,1)))
     
      IF (Diff .LT. Diff0) THEN 
        Diff1=Diff
10      IF (Diff1 .GE. 0.05) THEN
        dT=dT+dsT     
        T(j+1,1)=T(j+1,2)-dT 
        Esat=0.6112*Exp(17.67*(T(j+1,1)-273.16)/(T(j+1,1)-29.66))
        Qsat=0.622*Esat/Pa
        Q(j+1,1)=RH*Qsat
        Diff1=ABS(((T(j+1,2)-T(j+1,1))/dz(1))+Le/Cp*((Q(j+1,2)-Q(j+1,1))
     +    /dz(1))-0.74/(k*z(1)*us(j+1,1))*(-1.)*(tw(1,1)+Le/Cp*qw(1,1)))
        GO TO 10
        ENDIF
       ELSE      
        Diff2=Diff0
20      IF (Diff2 .GE. 0.05) THEN
        dT=dT-dsT      
        T(j+1,1)=T(j+1,2)-dT 
        Esat=0.6112*Exp(17.67*(T(j+1,1)-273.16)/(T(j+1,1)-29.66))
        Qsat=0.622*Esat/Pa
        Q(j+1,1)=RH*Qsat
        Diff2=ABS(((T(j+1,2)-T(j+1,1))/dz(1))+Le/Cp*((Q(j+1,2)-Q(j+1,1))
     +    /dz(1))-0.74/(k*z(1)*us(j+1,1))*(-1.)*(tw(1,1)+Le/Cp*qw(1,1))) 
        GO TO 20
        ENDIF  
      ENDIF  
      
      Q(j+1,2)=Q(j+1,2)*1000.
      Q(j+1,1)=Q(j+1,1)*1000.  
      qw(1,1)=qw(1,1)*1000.
C     *****Transform unit for Q and qw*********************************
C     *****Calculate the other parameters by MOST**********************
c      ts(j+1,1)=(T(j+1,2)-T(j+1,1))*k/0.74/(Log(z(2)/zo2))
c      qs(j+1,1)=(Q(j+1,2)-Q(j+1,1))*k/0.74/(Log(z(2)/zo2))
      ts(j+1,1)=k*z(1)/0.74*((T(j+1,2)-T(j+1,1))/(dz(1)))
      qs(j+1,1)=k*z(1)/0.74*((Q(j+1,2)-Q(j+1,1))/(dz(1)))
      tw(j+1,1)=-1.*ts(j+1,1)*us(j+1,1)
      qw(j+1,1)=-1.*qs(j+1,1)*us(j+1,1)
       U(j+1,1)=us(j+1,1)/k*Log(z(1)/zo2)
       W(j+1,1)=0.0 
      uu(j+1,1)=au*(us(j+1,1))**2.
      vv(j+1,1)=av*(us(j+1,1))**2.
      ww(j+1,1)=aw*(us(j+1,1))**2.
      uw(j+1,1)=-1.*us(j+1,1)**2. 
      tt(j+1,1)=at*(ts(j+1,1))**2.
      qq(j+1,1)=aq*(qs(j+1,1))**2.
      tu(j+1,1)=atu*ts(j+1,1)*us(j+1,1)
      qu(j+1,1)=aqu*qs(j+1,1)*us(j+1,1)
      qt(j+1,1)=aqt*qs(j+1,1)*ts(j+1,1) 
      e2(j+1,1)=uu(j+1,1)+vv(j+1,1)+ww(j+1,1) 
       e(j+1,1)=ke*(e2(j+1,1)**(3./2.))/(k*z(1)) 
       
      us(j+1,2)=(-1.*uw(j+1,2))**0.5
c      write (*,*) 'j,T(j+1,1),Q(j+1,1)',j,T(j+1,1),Q(j+1,1)  
c      write (*,*) 'U(j+1,2),U(j+1,1)',U(j+1,2),U(j+1,1)
c      write (*,*) 'us,tw,qw(j+1,2)',us(j+1,2),tw(j+1,2),qw(j+1,2)
c      write (*,*) 'us,tw,qw(j+1,1)',us(j+1,1),tw(j+1,1),qw(j+1,1)
c      pause  
C     -----------------------------------------------------------------
      ENDDO                             
C     -----------------------------------------------------------------      
C     *****Writing the Outputs***************************************** 
      OPEN (1,FILE='testfa.01')
c      DO j=1,nx,1000
      DO n=1,nz
c      WRITE (1,50) x(j),z(n),U(j,n),W(j,n),T(j,n),Q(j,n),
c     +              uw(j,n),ww(j,n),tw(j,n),qw(j,n) 
c50    FORMAT (2(F11.4),1X,8(F12.4))  
      WRITE (1,60) x(1),z(n),U(1,n),W(1,n),T(1,n),Q(1,n),
     +             uw(1,n),ww(1,n),tw(1,n),qw(1,n),
     +             x(10),z(n),U(10,n),W(10,n),T(10,n),Q(10,n),
     +             uw(10,n),ww(10,n),tw(10,n),qw(10,n), 
     +             x(100),z(n),U(100,n),W(100,n),T(100,n),Q(100,n),
     +             uw(100,n),ww(100,n),tw(100,n),qw(100,n), 
     +             x(1000),z(n),U(1000,n),W(1000,n),T(1000,n),Q(1000,n),
     +             uw(1000,n),ww(1000,n),tw(1000,n),qw(1000,n) 
c     +             x(2000),z(n),U(2000,n),W(2000,n),T(2000,n),Q(2000,n),
c     +             uw(2000,n),ww(2000,n),tw(2000,n),qw(2000,n) 
c     +             x(4000),z(n),U(4000,n),W(4000,n),T(4000,n),Q(4000,n),
c     +             uw(4000,n),ww(4000,n),tw(4000,n),qw(4000,n) 
60    FORMAT (40(F12.4))      
      ENDDO                       
c      WRITE (1,*) '---------------------------------------------------'
c      WRITE (1,*) '---------------------------------------------------' 
c      ENDDO
      WRITE (*,*) '=====> Program Finished'   
         
      END 
      
               
