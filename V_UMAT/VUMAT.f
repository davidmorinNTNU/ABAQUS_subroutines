!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine VUMAT
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      INCLUDE './UMAT_MODEL.f'
!-----------------------------------------------------------------------
      SUBROUTINE VUMAT(
     + NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME,
     + TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY,
     + STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRAdoLD, FIELdoLD,
     + STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW,
     + STRETCHNEW, DEFGRADNEW, FIELDNEW,
     + STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
C
      INCLUDE 'VABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      character*(*) CMNAME
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     + CHARLENGTH(*), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(*),
     + TEMPOLD(*), FIELdoLD(NBLOCK,NFIELDV),FIELDNEW(NBLOCK,NFIELDV),
     + STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV),
     + ENERINTERNOLD(NBLOCK),  ENERINELASOLD(NBLOCK), TEMPNEW(*),
     + STRETCHOLD(NBLOCK,NDIR+NSHR), DEFGRADOLD(NBLOCK,NDIR+2*NSHR),
     + STRETCHNEW(NBLOCK,NDIR+NSHR), DEFGRADNEW(NBLOCK,NDIR+2*NSHR),
     + STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     + ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
!-----------------------------------------------------------------------
!-----Declaration parameters
!-----------------------------------------------------------------------
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     1 third = one/three, half = .5, twoThirds = two/three,
     2 threeHalfs = 1.5 )
!-----------------------------------------------------------------------
!-----Declaration elastic parameters
!-----------------------------------------------------------------------
      real*8 BULK,R2G,E0,NU
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i,k
      integer KINC,NTENS,NSTATV
      real*8 PROPS2(NPROPS),STRESS(NDIR+NSHR),DSTRAN(NDIR+NSHR)
      real*8 STATEV(NSTATEV),DDSDDE(NDIR+NSHR,NDIR+NSHR)
!-----------------------------------------------------------------------
!     Initialization step (elastic)
!-----------------------------------------------------------------------
      if((steptime.eq.totaltime).and.(steptime.eq.zero))then
         E0   = props(1)
         NU   = props(2)
         BULK = E0/(3.0*(1.0-2.0*NU))
         R2G  = E0/(1.0+NU)
         do i=1,NBLOCK
            if(NSHR.eq.3)then
               volstrain = STRAININC(i,1)+STRAININC(i,2)+STRAININC(i,3)
               STRESSNEW(i,1) = R2G*(STRAININC(i,1)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,2) = R2G*(STRAININC(i,2)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,3) = R2G*(STRAININC(i,3)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,4) = R2G*STRAININC(i,4)
               STRESSNEW(i,5) = R2G*STRAININC(i,5)
               STRESSNEW(i,6) = R2G*STRAININC(i,6)
            else
               volstrain = STRAININC(i,1)+STRAININC(i,2)+STRAININC(i,3)
               STRESSNEW(i,1) = R2G*(STRAININC(i,1)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,2) = R2G*(STRAININC(i,2)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,3) = R2G*(STRAININC(i,3)-third*volstrain)
     +                         +BULK*volstrain
               STRESSNEW(i,4) = R2G*STRAININC(i,4)
            endif
         enddo
!-----------------------------------------------------------------------
!     Ordinary increment
!-----------------------------------------------------------------------
      else
         if(totaltime.eq.DT)then
            KINC = 1
         else
            KINC = 2
         endif
         NSTATV = NSTATEV
         NTENS  = NDIR+NSHR
         PROPS2 = PROPS
!-----------------------------------------------------------------------
!        Loop over elements
!-----------------------------------------------------------------------
         do i=1,NBLOCK
!-----------------------------------------------------------------------
c           GRAB STRESSES            
!-----------------------------------------------------------------------
            STRESS(1) = STRESSOLD(i,1)
            STRESS(2) = STRESSOLD(i,2)
            STRESS(3) = STRESSOLD(i,3)
            STRESS(4) = STRESSOLD(i,4)
            if(NSHR.eq.3)then
               STRESS(5) = STRESSOLD(i,6)
               STRESS(6) = STRESSOLD(i,5)
            endif
!-----------------------------------------------------------------------
c           GRAB STRAININC
!-----------------------------------------------------------------------
            DSTRAN(1) = STRAININC(i,1)
            DSTRAN(2) = STRAININC(i,2)
            DSTRAN(3) = STRAININC(i,3)
            DSTRAN(4) = 2.0*STRAININC(i,4)
            if(NSHR.eq.3)then
               DSTRAN(5) = 2.0*STRAININC(i,6)
               DSTRAN(6) = 2.0*STRAININC(i,5)
            endif
!-----------------------------------------------------------------------
c           GRAB HISTORY VARIABLES
!-----------------------------------------------------------------------
            do k=1,NSTATEV
               STATEV(k) = STATEOLD(i,k)
            enddo
!-----------------------------------------------------------------------
c           CALL UMAT_MODEL
!-----------------------------------------------------------------------
            call UMAT_MODEL(DDSDDE,STRESS,STATEV,DSTRAN,PROPS2,
     +                      KINC,NTENS,NSTATV,NPROPS)
!-----------------------------------------------------------------------
c           UNPACK STRESSES
!-----------------------------------------------------------------------
            STRESSNEW(i,1) = STRESS(1)
            STRESSNEW(i,2) = STRESS(2)
            STRESSNEW(i,3) = STRESS(3)
            STRESSNEW(i,4) = STRESS(4)
            if(NSHR.eq.3)then
               STRESSNEW(i,6) = STRESS(5)
               STRESSNEW(i,5) = STRESS(6)
            endif
!-----------------------------------------------------------------------
c           UNPACK HISTORY VARIABLES
!-----------------------------------------------------------------------
            do k=1,NSTATEV
               STATENEW(i,k) = STATEV(k)
            enddo
!-----------------------------------------------------------------------
c           UPDATE INTERNAL ENERGIES
!-----------------------------------------------------------------------
            if(NSHR.eq.3)then
               enerInternNew(i) = enerInternOld(i)
     +     +((stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     +     +(stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     +     +(stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     +     +2.0*(stressOld(i,4)+stressNew(i,4))*strainInc(i,4)
     +     +2.0*(stressOld(i,5)+stressNew(i,5))*strainInc(i,5)
     +     +2.0*(stressOld(i,6)+stressNew(i,6))*strainInc(i,6))
     +     *0.5/DENSITY(i)
            else
               enerInternNew(i) = enerInternOld(i)
     +     +((stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     +     +(stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     +     +(stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     +     +2.0*(stressOld(i,4)+stressNew(i,4))*strainInc(i,4))
     +     *0.5/DENSITY(i)
            endif
         enddo
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end