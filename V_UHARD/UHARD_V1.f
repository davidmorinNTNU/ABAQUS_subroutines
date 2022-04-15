!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine UHARD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,
     + TIME,DTIME,TEMP,DTEMP,
     + NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     + STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      DIMENSION HARD(3),STATEV(NSTATV),TIME(*),
     +          PREDEF(NUMFIELDV),DPRED(*),PROPS(*)
      CHARACTER*80 CMNAME
!-----------------------------------------------------------------------
!-----Declaration material parameters
!-----------------------------------------------------------------------
      real*8 SIGMA0,T1,Q1,T2,Q2,T3,Q3
      real*8 C,PDOT0
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      real*8 T1oQ1,T2oQ2,T3oQ3
      real*8 vp,dvp
!-----------------------------------------------------------------------
!     Beginning of subroutine
!-----------------------------------------------------------------------
c     Isotropic work-hardening
      SIGMA0 = props(1)
      T1     = props(2)
      Q1     = props(3)
      T2     = props(4)
      Q2     = props(5)
      T3     = props(6)
      Q3     = props(7)
c     Strain rate sensitivity
      C      = props(9)
      PDOT0  = props(10)
!-----------------------------------------------------------------------
!     Compute hardening constants to avoid division by zero
!-----------------------------------------------------------------------
      if(abs(Q1).gt.0.0)then
         T1oQ1 = T1/Q1
      else
         T1oQ1 = 0.0
      endif
c
      if(abs(Q2).gt.0.0)then
         T2oQ2 = T2/Q2
      else
         T2oQ2 = 0.0
      endif
c
      if(abs(Q3).gt.0.0)then
         T3oQ3 = T3/Q3
      else
         T3oQ3 = 0.0
      endif
!-----------------------------------------------------------------------
!     Compute yield stress and derivatives
!-----------------------------------------------------------------------
      SYIELD = (SIGMA0+Q1*(1.0-exp(-T1oQ1*EQPLAS))
     +                +Q2*(1.0-exp(-T2oQ2*EQPLAS))
     +                +Q3*(1.0-exp(-T3oQ3*EQPLAS)))
!-----------------------------------------------------------------------
!     Compute yield stress derivatives
!-----------------------------------------------------------------------
c     Derivative with respect to equivalent plastic strain   
      HARD(1) =(T1*exp(-T1oQ1*EQPLAS)
     +         +T2*exp(-T2oQ2*EQPLAS)
     +         +T3*exp(-T3oQ3*EQPLAS))
c     Derivative with respect to equivalent plastic strain rate
      HARD(2) = 0.0
c     Derivative with respect to temperature
      HARD(3) = 0.0
!-----------------------------------------------------------------------
!     Apply rate-sensitivity if needed
!-----------------------------------------------------------------------
      if((C.gt.0.0).and.(PDOT0.gt.0.0))then
c         Compute rate effects and its derivative    
          vp  = (1.0+EQPLASRT/PDOT0)**C
          dvp = (1.0+EQPLASRT/PDOT0)**(C-1.0)*(C/PDOT0)
c         Apply rate effects to the required quantities
          HARD(2) = SYIELD*dvp
          HARD(1) = HARD(1)*vp
          SYIELD  = SYIELD*vp
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END