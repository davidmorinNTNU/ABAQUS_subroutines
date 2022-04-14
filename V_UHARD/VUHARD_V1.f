!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine VUHARD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE VUHARD(NBLOCK,JELEM,KINTPT,KLAYER,KSECPT,
     . LANNEAL,STEPTIME,TOTALTIME,DT,CMNAME,
     . NSTATEV,NFIELDV,NPROPS,
     . PROPS,TEMPOLD,TEMPNEW,FIELDOLD,FIELDNEW,
     . STATEOLD,
     . EQPS,EQPSRATE,
     . YIELD,DYIELDDTEMP,DYIELDDEQPS,
     . STATENEW)
      include 'vaba_param.inc'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      dimension PROPS(NPROPS), TEMPOLD(NBLOCK), TEMPNEW(NBLOCK),
     . FIELDOLD(NBLOCK,NFIELDV), FIELDNEW(NBLOCK,NFIELDV),
     . STATEOLD(NBLOCK,NSTATEV), EQPS(NBLOCK), EQPSRATE(NBLOCK),
     . YIELD(NBLOCK), DYIELDDTEMP(NBLOCK), DYIELDDEQPS(NBLOCK,2),
     . STATENEW(NBLOCK,NSTATEV), JELEM(NBLOCK)
      character*80 CMNAME
!-----------------------------------------------------------------------
!-----Declaration material parameters
!-----------------------------------------------------------------------
      real*8 SIGMA0,T1,Q1,T2,Q2,T3,Q3
      real*8 C,PDOT0
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 T1oQ1,T2oQ2,T3oQ3
      real*8 vp,dvp
!-----------------------------------------------------------------------
!     Read material properties
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
!     Compute yield stress and its derivatives
!-----------------------------------------------------------------------
      do i=1,nblock
         yield(i) =(SIGMA0+Q1*(1.0-exp(-T1oQ1*eqps(i)))
     +                    +Q2*(1.0-exp(-T2oQ2*eqps(i)))
     +                    +Q3*(1.0-exp(-T3oQ3*eqps(i))))
c        Derivative with respect to equivalent plastic strain   
         dyieldDeqps(i,1) =(T1*exp(-T1oQ1*eqps(i))
     +                     +T2*exp(-T2oQ2*eqps(i))
     +                     +T3*exp(-T3oQ3*eqps(i)))
c        Derivative with respect to equivalent plastic strain rate
         dyieldDeqps(i,2) = 0.0
c        Derivative with respect to temperature
         dyieldDtemp(i)   = 0.0
      enddo
!-----------------------------------------------------------------------
!     Apply rate-sensitivity if needed
!-----------------------------------------------------------------------
      if((C.gt.0.0).and.(PDOT0.gt.0.0))then
         do i=1,nblock
c           Compute rate effects and its derivative    
            vp  = (1.0+eqpsRate(i)/PDOT0)**C
            dvp = (1.0+eqpsRate(i)/PDOT0)**(C-1.0)*(C/PDOT0)
c           Apply rate effects to the required quantities
            dyieldDeqps(i,2) = yield(i)*dvp
            dyieldDeqps(i,1) = dyieldDeqps(i,1)*vp
            yield(i) = yield(i)*vp
         enddo
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
