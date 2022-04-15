!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine UMAT_MODEL
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine UMAT_MODEL(DDSDDE,STRESS,STATEV,DSTRAN,PROPS,KINC,
     +                      NTENS,NSTATV,NPROPS)
      implicit none
!-----------------------------------------------------------------------
!-----Declaration variables
!-----------------------------------------------------------------------
      real*8 STRESS(NTENS),STATEV(NSTATV),DSTRAN(NTENS),PROPS(NPROPS)
      real*8 DDSDDE(NTENS,NTENS)
      integer NTENS,NSTATV,NPROPS,KINC
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i,j
!DECLARE -- Strain increments
      real*8 d1,d2,d3,d4,d5,d6
!DECLARE -- Stress tensor components
      real*8 s1,s2,s3,s4,s5,s6
!DECLARE -- Old stress tensor components
      real*8 s1old,s2old,s3old,s4old,s5old,s6old
!DECLARE -- Trial stress tensor components
      real*8 t1,t2,t3,t4,t5,t6
!DECLARE -- Equivalent plastic strain
      real*8 p,pold
!DECLARE -- Elasticity constants
      real*8 E0,NU
!DECLARE -- Elasticity matrix
      real*8 Ce(6,6)
!DECLARE -- Yield surface parameters
      real*8 SIGMAY,dSIGMAYdp
!DECLARE -- Yield function
      real*8 YF
!DECLARE -- Equivalent stress
      real*8 PHI
!DECLARE -- Yield function derivatives
      real*8 dfds1,dfds2,dfds3,dfds4,dfds5,dfds6,con
      real*8 dyfds,dyfdL,dfdSY
!DECLARE -- Work-hardening parameters
      real*8 SIGMA0,TR1,QR1,TR2,QR2,TR3,QR3
      real*8 T1oQ1,T2oQ2,T3oQ3
!DECLARE -- Fracture parameters
      real*8 WC,DCRIT
!DECLARE -- Damage variables
      real*8 tr
      real*8 e1,e2,e3,e4,e5,e6
      real*8 f1,f2,f3,f4
      real*8 SIG1P
      real*8 damage,damageold
!DECLARE -- RMAP variables
      real*8 TOL
      parameter(TOL=1e-8)
      real*8 de1,de2,de3,de4,de5,de6
      real*8 dL,ddL
      real*8 temp1,temp2,temp3,temp4,temp5,temp6
      integer inum
!DECLARE -- Data check
      integer iflUMAT
      data iflUMAT/0/
!-----------------------------------------------------------------------
!     Read material parameters
!-----------------------------------------------------------------------
c     Elasticity parameters
      E0 = PROPS(1)
      NU = PROPS(2)
c     Work-hardening parameters
      SIGMA0 = PROPS(9)
      TR1    = PROPS(10)
      QR1    = PROPS(11)
      TR2    = PROPS(12)
      QR2    = PROPS(13)
      TR3    = PROPS(14)
      QR3    = PROPS(15)
c     Failure parameters
      WC    = PROPS(17)
      DCRIT = PROPS(18)
!-----------------------------------------------------------------------
!     Define elasticity matrix
!-----------------------------------------------------------------------
      Ce = 0.0
      Ce(1,1) = E0*(1.0-NU)/((1.0+NU)*(1.0-2.0*NU))
      Ce(1,2) = NU*Ce(1,1)/(1.0-NU)
      Ce(1,3) = Ce(1,2)
      Ce(2,1) = Ce(1,2)
      Ce(2,2) = E0*(1.0-NU)/((1.0+NU)*(1.0-2.0*NU))
      Ce(2,3) = Ce(1,2)
      Ce(3,1) = Ce(1,2)
      Ce(3,2) = Ce(1,2)
      Ce(3,3) = E0*(1.0-NU)/((1.0+NU)*(1.0-2.0*NU))
      Ce(4,4) = 0.5*E0/(1.0+NU)
      Ce(5,5) = 0.5*E0/(1.0+NU)
      Ce(6,6) = 0.5*E0/(1.0+NU)
!-----------------------------------------------------------------------
!     Define additional hardening parameters
!-----------------------------------------------------------------------
      IF(abs(QR1).gt.0.0)THEN
         T1oQ1 = TR1/QR1
      ELSE
         T1oQ1 = 0.0
      ENDIF
c
      IF(abs(QR2).gt.0.0)THEN
         T2oQ2 = TR2/QR2
      ELSE
         T2oQ2 = 0.0
      ENDIF
c
      IF(abs(QR3).gt.0.0)THEN
         T3oQ3 = TR3/QR3
      ELSE
         T3oQ3 = 0.0
      ENDIF
!-----------------------------------------------------------------------
!     Print parameters when time = 0.0
!-----------------------------------------------------------------------
      if(iflUMAT.eq.0)then
         write(*,'(a30)') '______________________________'
         write(*,'(a30)') '|                            |'
         write(*,'(a30)') '|        SIMLab Example      |'
         write(*,'(a30)') '|            UMAT            |'
         write(*,'(a30)') '|____________________________|'
         write(*,*)
         write(*,'(a10,es14.6)') 'E0       =', E0
         write(*,'(a10,es14.6)') 'NU       =', NU
         write(*,'(a10,es14.6)') 'SIGMA0   =', SIGMA0
         write(*,'(a10,es14.6)') 'THETAR1  =', TR1
         write(*,'(a10,es14.6)') 'QR1      =', QR1
         write(*,'(a10,es14.6)') 'THETAR2  =', TR2
         write(*,'(a10,es14.6)') 'QR2      =', QR2
         write(*,'(a10,es14.6)') 'THETAR3  =', TR3
         write(*,'(a10,es14.6)') 'QR3      =', QR3
         write(*,'(a10,es14.6)') 'WC       =', WC
         write(*,'(a10,es14.6)') 'DCRIT    =', DCRIT
         write(*,*)
         iflUMAT = 1
      endif
!-----------------------------------------------------------------------
!     Packaging strain
!-----------------------------------------------------------------------
      d1 = DSTRAN(1)
      d2 = DSTRAN(2)
      d4 = DSTRAN(4)
c
      IF(NTENS.eq.6)THEN
         d3 = DSTRAN(3)
         d5 = DSTRAN(6)
         d6 = DSTRAN(5)
      ELSE
         d3 = DSTRAN(3)
         d5 = 0.0
         d6 = 0.0
      ENDIF
!-----------------------------------------------------------------------
!     Grab old stresses
!-----------------------------------------------------------------------
      s1old = STRESS(1)
      s2old = STRESS(2)
      s3old = STRESS(3)
      s4old = STRESS(4)
      IF(NTENS.eq.6)THEN
         s5old = STRESS(6)
         s6old = STRESS(5)
      ELSE
         s5old = 0.0
         s6old = 0.0
      ENDIF
!-----------------------------------------------------------------------
!     Extract history variables
!-----------------------------------------------------------------------
      pold      = STATEV(1)
      damageold = STATEV(2) 
!-----------------------------------------------------------------------
!     Compute elastic prediction
!-----------------------------------------------------------------------
      t1 = s1old+(Ce(1,1)*d1+Ce(1,2)*d2+Ce(1,3)*d3)
      t2 = s2old+(Ce(2,1)*d1+Ce(2,2)*d2+Ce(2,3)*d3)
      t3 = s3old+(Ce(3,1)*d1+Ce(3,2)*d2+Ce(3,3)*d3)
      t4 = s4old+(Ce(4,4)*d4)
      t5 = s5old+(Ce(5,5)*d5)
      t6 = s6old+(Ce(6,6)*d6)
!-----------------------------------------------------------------------
!     Compute equivalent stress
!-----------------------------------------------------------------------
      PHI = (t1*t1+t2*t2+t3*t3-t1*t2-t2*t3-t3*t1
     +     +3.0*(t4*t4+t5*t5+t6*t6))
      if(PHI.gt.0.0)then
         PHI = sqrt(PHI)
      else
         PHI = 0.0
      endif
!-----------------------------------------------------------------------
!     Compute Yield stress
!-----------------------------------------------------------------------
      SIGMAY = SIGMA0+QR1*(1.0-exp(-T1oQ1*pold))
     +               +QR2*(1.0-exp(-T2oQ2*pold))
     +               +QR3*(1.0-exp(-T3oQ3*pold))
!-----------------------------------------------------------------------
!     Compute Yield function
!-----------------------------------------------------------------------
      YF = PHI-SIGMAY
!-----------------------------------------------------------------------
!     CHECK FOR PLASTICITY
!-----------------------------------------------------------------------
      IF(YF.lt.0.0)THEN! ELASTICITY
         s1 = t1
         s2 = t2
         s3 = t3
         s4 = t4
         s5 = t5
         s6 = t6
c
         p   = pold
         dL  = 0.0
      ELSE ! PLASTICITY
         de1 = 0.0
         de2 = 0.0
         de3 = 0.0
         de4 = 0.0
         de5 = 0.0
         de6 = 0.0
c
         dL   = 0.0
         ddL  = 0.0
         p    = pold
c
         s1 = t1
         s2 = t2
         s3 = t3
         s4 = t4
         s5 = t5
         s6 = t6
c
         inum = 0
!-----------------------------------------------------------------------
!     Solve nonlinear equation in plastic parameter
!-----------------------------------------------------------------------
  100 continue
!-----------------------------------------------------------------------
!     Update gradients
!-----------------------------------------------------------------------
         con   = 1.0/PHI
         dfds1 = con*(s1-0.5*(s2+s3))
         dfds2 = con*(s2-0.5*(s3+s1))
         dfds3 = con*(s3-0.5*(s1+s2))
         con   = 3.0*con
         dfds4 = con*s4
         dfds5 = con*s5
         dfds6 = con*s6
!-----------------------------------------------------------------------
!     Update Isotropic hardening
!-----------------------------------------------------------------------
         SIGMAY = SIGMA0+QR1*(1.0-exp(-T1oQ1*p))
     +                  +QR2*(1.0-exp(-T2oQ2*p))
     +                  +QR3*(1.0-exp(-T3oQ3*p))
         dSIGMAYdp = TR1*exp(-T1oQ1*p)
     +              +TR2*exp(-T2oQ2*p)
     +              +TR3*exp(-T3oQ3*p)
!-----------------------------------------------------------------------
!     Update Cutting plane derivatives
!-----------------------------------------------------------------------
         temp1 = Ce(1,1)*dfds1+Ce(1,2)*dfds2+Ce(1,3)*dfds3
         temp2 = Ce(2,1)*dfds1+Ce(2,2)*dfds2+Ce(2,3)*dfds3
         temp3 = Ce(3,1)*dfds1+Ce(3,2)*dfds2+Ce(3,3)*dfds3
         temp4 = Ce(4,4)*dfds4
         temp5 = Ce(5,5)*dfds5
         temp6 = Ce(6,6)*dfds6
c
         dyfds = temp1*dfds1+temp2*dfds2
     +          +temp3*dfds3+temp4*dfds4
     +          +temp5*dfds5+temp6*dfds6
c
         dyfdL = dSIGMAYdp
!-----------------------------------------------------------------------
!     Update yield function and equivalent plastic strain increment
!-----------------------------------------------------------------------
         YF  = PHI-SIGMAY
         ddL = YF/(dyfds+dyfdL)
         dL  = dL+ddL
         p   = p+ddL
!-----------------------------------------------------------------------
!     Update Cutting plane increment
!-----------------------------------------------------------------------
         de1 = de1+ddL*dfds1
         de2 = de2+ddL*dfds2
         de3 = de3+ddL*dfds3
         de4 = de4+ddL*dfds4
         de5 = de5+ddL*dfds5
         de6 = de6+ddL*dfds6
!-----------------------------------------------------------------------
!     Update Stress tensor
!-----------------------------------------------------------------------
         s1 = t1-(Ce(1,1)*de1+Ce(1,2)*de2+Ce(1,3)*de3)
         s2 = t2-(Ce(2,1)*de1+Ce(2,2)*de2+Ce(2,3)*de3)
         s3 = t3-(Ce(3,1)*de1+Ce(3,2)*de2+Ce(3,3)*de3)
         s4 = t4-Ce(4,4)*de4
         s5 = t5-Ce(5,5)*de5
         s6 = t6-Ce(6,6)*de6
!-----------------------------------------------------------------------
!     Update equivalent stress
!-----------------------------------------------------------------------
         PHI = sqrt((s1*s1+s2*s2+s3*s3-s1*s2-s2*s3-s3*s1
     +         +3.0*(s4*s4+s5*s5+s6*s6)))
!-----------------------------------------------------------------------
!     Check for convergence
!-----------------------------------------------------------------------
         inum = inum + 1
c
         IF(inum.gt.1000)THEN
            print*,'NO CONVERGENCE'
            print*,'',PHI,YF
            stop
         ENDIF
c
         IF((YF.gt.TOL*SIGMAY).or.(inum.eq.1))THEN
            GOTO 100
         ELSE
c
         ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Update Damage
!-----------------------------------------------------------------------
      if(dL.gt.0.0)then
         tr = (s1+s2+s3)/3.0
         e1 = s1-tr
         e2 = s2-tr
         e3 = s3-tr
         e4 = s4
         e5 = s5
         e6 = s6
         f1 = 0.5d0*(e1*e1+e2*e2+e3*e3)+e4*e4+e5*e5+e6*e6
         f2 = e1*e5*e5+e2*e6*e6+e3*e4*e4-e1*e2*e3-2.0d0*e4*e5*e6
         f3 =-sqrt(27.d0/f1)*f2*0.5/f1
         f3 = sign(min(abs(f3),1.0d0),f3)
         f4 = acos(f3)/3.0
c
         SIG1P  = tr+2.0*sqrt(f1/3.0)*cos(f4)
         damage = damageold+max(SIG1P,0.0)*(dL/WC)
      else
         damage = damageold
      endif
      STRESS(1) = s1
      STRESS(2) = s2
      STRESS(3) = s3
      STRESS(4) = s4
      IF(NTENS.eq.6)THEN
         STRESS(6) = s5
         STRESS(5) = s6
      ENDIF
!-----------------------------------------------------------------------
!     Update Stresses
!-----------------------------------------------------------------------
      STRESS(1) = s1
      STRESS(2) = s2
      STRESS(3) = s3
      STRESS(4) = s4
      IF(NTENS.eq.6)THEN
         STRESS(6) = s5
         STRESS(5) = s6
      ENDIF
!-----------------------------------------------------------------------
!     Update Equivalent plastic strain
!-----------------------------------------------------------------------
      STATEV(1) = p
      STATEV(2) = damage
!-----------------------------------------------------------------------
!     CHECK FOR FRACTURE
!-----------------------------------------------------------------------
      IF(STATEV(2).ge.DCRIT)THEN
         STATEV(NSTATV) = 0.0
      ELSE
         STATEV(NSTATV) = 1.0
      ENDIF
!-----------------------------------------------------------------------
!     Update Consistent tangent operator
!-----------------------------------------------------------------------
      DO i=1,NTENS
         DO j=1,NTENS
            DDSDDE(i,j) = Ce(i,j)
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
