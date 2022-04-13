!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine VUSDFLD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine VUSDFLD(
     + NBLOCK, NSTATEV, NFIELDV, NPROPS, NDIR, NSHR,
     + JELEM, KINTPT, KLAYER, KSECPT,
     + STEPTIME, TOTALTIME, DT, CMNAME,
     + COORDMP, DIRECT, T, CHARLENGTH, PROPS,
     + STATEOLD,STATENEW,FIELD)
      include 'vaba_param.inc'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      dimension JELEM(nblock),COORDMP(nblock,*),DIRECT(nblock,3,3),
     .          T(nblock,3,3),CHARLENGTH(nblock),PROPS(nprops),
     .          STATEOLD(nblock,nstatev),STATENEW(nblock,nstatev),
     .          FIELD(nblock,nfieldv)
      character*80 cmname
!-----Data from ABAQUS
      dimension stressdata(maxblk*(NDIR+NSHR))
      integer jSData(maxblk*(NDIR+NSHR))
      character*3 cSData(maxblk*(NDIR+NSHR))
      integer jStatus
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 s(NBLOCK,NDIR+NSHR)
      real*8 SMISES(NBLOCK),SIGH(NBLOCK),TRIAX(NBLOCK)
!-----------------------------------------------------------------------
!     Access stress tensor
!-----------------------------------------------------------------------
      call vgetvrm( 'S', stressdata,jSData,cSData,jStatus)
!-----------------------------------------------------------------------
!     Extract data from stressdata
!-----------------------------------------------------------------------
      if(nshr.gt.1)then
         do i=1,nblock
            s(i,1) = stressdata(i)
            s(i,2) = stressdata(i+nblock)
            s(i,3) = stressdata(i+nblock*2)
            s(i,4) = stressdata(i+nblock*3)
            s(i,5) = stressdata(i+nblock*4)
            s(i,6) = stressdata(i+nblock*5)
         enddo
      else
         do i=1,nblock
            s(i,1) = stressdata(i)
            s(i,2) = stressdata(i+nblock)
            s(i,3) = stressdata(i+nblock*2)
            s(i,4) = stressdata(i+nblock*3)
         enddo
      endif
!-----------------------------------------------------------------------
!     Compute von Mises equivalent stress and hydrostatic stress
!-----------------------------------------------------------------------
      if(nshr.gt.1)then
         do i=1,nblock
            SIGH(i)   = (s(i,1)+s(i,2)+s(i,3))/3.0
            SMISES(i) = sqrt(s(i,1)*s(i,1)+s(i,2)*s(i,2)
     +                      +s(i,3)*s(i,3)
     +                      -s(i,1)*s(i,2)-s(i,2)*s(i,3)
     +                      -s(i,3)*s(i,1)
     +                 +3.0*(s(i,4)*s(i,4)+s(i,5)*s(i,5)
     +                      +s(i,6)*s(i,6)))
         enddo
      else
         do i=1,nblock
            SIGH(i)   = (s(i,1)+s(i,2)+s(i,3))/3.0
            SMISES(i) = sqrt(s(i,1)*s(i,1)+s(i,2)*s(i,2)
     +                      -s(i,1)*s(i,2)+3.0*s(i,4)*s(i,4))
         enddo
      endif
!-----------------------------------------------------------------------
!     Compute stress triaxiality
!-----------------------------------------------------------------------
      do i=1,nblock
         if(SMISES(i).gt.0.0)then
            TRIAX(i) = SIGH(i)/SMISES(i)
         else
            TRIAX(i) = 0.0
         endif
      enddo
!-----------------------------------------------------------------------
!     Update field variable
!-----------------------------------------------------------------------
      do i=1,nblock
         if(TRIAX(i).ge.0.0)then
            field(i,1) = 1.05
         else
            field(i,1) =-1.05
         endif
      enddo
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
