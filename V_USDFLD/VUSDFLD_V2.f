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
      dimension JELEM(NBLOCK),COORDMP(NBLOCK,*),DIRECT(NBLOCK,3,3),
     .          T(NBLOCK,3,3),CHARLENGTH(NBLOCK),PROPS(NPROPS),
     .          STATEOLD(NBLOCK,NSTATEV),STATENEW(NBLOCK,NSTATEV),
     .          FIELD(NBLOCK,NFIELDV)
      character*80 CMNAME
!-----Data from ABAQUS
      dimension stressdata(maxblk*(ndir+nshr))
      integer jSData(maxblk*(ndir+nshr))
      character*3 cSData(maxblk*(ndir+nshr))
      integer jStatus
      dimension peeqdata(maxblk)
      integer jPData(maxblk)
      character*3 cPData(maxblk)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 s(NBLOCK,NDIR+NSHR),SIG1(NBLOCK)
      real*8 eigVal(NBLOCK,3)
      real*8 damage(NBLOCK),dp(NBLOCK)
!-----------------------------------------------------------------------
!-----Declaration material parameters
!-----------------------------------------------------------------------
      real*8 Wc
!-----------------------------------------------------------------------
!     Read material properties
!-----------------------------------------------------------------------
      Wc = props(1)
!-----------------------------------------------------------------------
!     Access stress tensor
!-----------------------------------------------------------------------
      call vgetvrm( 'S', stressdata,jSData,cSData,jStatus)
!-----------------------------------------------------------------------
!     Access equivalent plastic strain
!-----------------------------------------------------------------------
      call vgetvrm('PEEQ',   peeqdata,jPData,cPData,jStatus)
!-----------------------------------------------------------------------
!     Extract data
!-----------------------------------------------------------------------
      do i=1,nblock
         dp(i)     = peeqdata(i)-stateOld(i,1)
         damage(i) = stateOld(i,2)
      enddo
!-----------------------------------------------------------------------
!     Extract data from stressdata
!-----------------------------------------------------------------------
      if(NSHR.gt.1)then
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
!     Compute principal stress
!-----------------------------------------------------------------------
      call vsprinc(nblock,s,eigVal,NDIR,NSHR)
      do i=1,nblock
         SIG1(i) = max(eigVal(i,1),eigVal(i,2),eigVal(i,3))
      enddo
!-----------------------------------------------------------------------
!     Update damage variable
!-----------------------------------------------------------------------
      do i=1,nblock
         if(dp(i).gt.0.0)then
            damage(i) = damage(i)+max(0.0,SIG1(i))*(dp(i)/Wc)
         endif
      enddo
!-----------------------------------------------------------------------
!     Update state dependent variables
!-----------------------------------------------------------------------
      do i=1,nblock
         stateNew(i,1) = peeqdata(i)
         stateNew(i,2) = min(1.0,damage(i))
      enddo
!-----------------------------------------------------------------------
!     Check for fracture
!-----------------------------------------------------------------------
      do i=1,nblock
         if(damage(i).ge.1.0)then
            stateNew(i,3) = 0.0
         else
            stateNew(i,3) = 1.0
         endif
      enddo
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end