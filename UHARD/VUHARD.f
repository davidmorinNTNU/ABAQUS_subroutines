!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine VUHARD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE VUHARD(
     . nblock,
     . jElem, kIntPt, kLayer, kSecPt,
     . lAnneal, stepTime, totalTime, dt, cmname,
     . nstatev, nfieldv, nprops,
     . props, tempOld, tempNew, fieldOld, fieldNew,
     . stateOld,
     . eqps, eqpsRate,
     . yield, dyieldDtemp, dyieldDeqps,
     . stateNew)
      include 'vaba_param.inc'
!-----Data from ABAQUS
      dimension props(nprops), tempOld(nblock), tempNew(nblock),
     . fieldOld(nblock,nfieldv), fieldNew(nblock,nfieldv),
     . stateOld(nblock,nstatev), eqps(nblock), eqpsRate(nblock),
     . yield(nblock), dyieldDtemp(nblock), dyieldDeqps(nblock,2),
     . stateNew(nblock,nstatev), jElem(nblock)
      character*80 cmname
!-----material parameters
      real*8 sigma0,THETAR1,QR1,THETAR2,QR2,THETAR3,QR3
      real*8 C,epsdot0,m,Tr,Tm,Teff
!-----Internal variables
      integer i
      real*8 T,THETAR1oQR1,THETAR2oQR2,THETAR3oQR3
      integer iflVUHARD
      data iflVUHARD/0/
!-----------------------------------------------------------------------
!     Read material properties
!-----------------------------------------------------------------------
      SIGMA0  = props(1)
      T1      = props(2)
      Q1      = props(3)
      T2      = props(4)
      Q2      = props(5)
      T3      = props(6)
      Q3      = props(7)
      C       = props(9)
      PDOT0   = props(10)
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
c
      do i=1,nblock
!-----------------------------------------------------------------------
!        Compute yield stress
!-----------------------------------------------------------------------
         yield(i) =(SIGMA0+Q1*(1.0-exp(-T1oQ1*eqps(i)))
     +                    +Q2*(1.0-exp(-T2oQ2*eqps(i)))
     +                    +Q3*(1.0-exp(-T3oQ3*eqps(i))))
     +                    *(1.0+eqpsRate(i)/PDOT0)**C
!-----------------------------------------------------------------------
!        Compute yield stress derivatives
!-----------------------------------------------------------------------
c     Derivative with respect to equivalent plastic strain   
         dyieldDeqps(i,1) =(T1*exp(-(T1oQ1)*eqps(i))
     +                     +T2*exp(-(T2oQ2)*eqps(i))
     +                     +T3*exp(-(T3oQ3)*eqps(i)))
     +                     *(1.0+eqpsRate(i)/PDOT0)**C
c     Derivative with respect to equivalent plastic strain rate
         dyieldDeqps(i,2) = (SIGMA0+Q1*(1.0-exp(-T1oQ1*eqps(i)))
     +                            +Q2*(1.0-exp(-T2oQ2*eqps(i)))
     +                            +Q3*(1.0-exp(-T3oQ3*eqps(i))))
     +                     *(1.0+eqpsRate(i)/PDOT0)**(C-1.0)*(C/PDOT0)
c     Derivative with respect to temperature
         dyieldDtemp(i)   = 0.0
      enddo
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
