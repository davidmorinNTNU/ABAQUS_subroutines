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
      real*8 PROPS2(NPROPS),STRESS(6),DSTRAN(6)
      real*8 STATEV(NSTATEV),DDSDDE(6,6)
!-----------------------------------------------------------------------
!-----Declaration plane stress conditions
!-----------------------------------------------------------------------
      integer iterb,maxiter
      parameter(maxiter=100)
      real*8 sig3o,eps3o,deps
      real*8 tol,tolalg,stol
      parameter(tol=1e-6)
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
         PROPS2 = PROPS
!-----------------------------------------------------------------------
!        Modified the number of tensor components
!-----------------------------------------------------------------------
         NTENS = 6
         NU    = props(2)
!-----------------------------------------------------------------------
!        Loop over elements
!-----------------------------------------------------------------------
         do i=1,NBLOCK
!-----------------------------------------------------------------------
c     Start iteration through the thickness
!-----------------------------------------------------------------------
            iterb = 1
            STRESS(3) = 0.0
            DSTRAN(3) = 0.0
  100 continue
!-----------------------------------------------------------------------
c           GRAB STRESSES            
!-----------------------------------------------------------------------
            STRESS(1) = STRESSOLD(i,1)
            STRESS(2) = STRESSOLD(i,2)
c            STRESS(3) = 0.0
            STRESS(4) = STRESSOLD(i,4)
            STRESS(5) = 0.0
            STRESS(6) = 0.0
!-----------------------------------------------------------------------
c           GRAB STRAININC
!-----------------------------------------------------------------------
            DSTRAN(1) = STRAININC(i,1)
            DSTRAN(2) = STRAININC(i,2)
c            DSTRAN(3) = 0.0
            DSTRAN(4) = 2.0*STRAININC(i,4)
            DSTRAN(5) = 0.0
            DSTRAN(6) = 0.0
!-----------------------------------------------------------------------
c           COMPUTE THROUGH-THICKNESS STRAIN INCREMENT
!-----------------------------------------------------------------------
            if(iterb.eq.1)then
               DSTRAN(3) = -(NU/(1.0-NU))*(DSTRAN(1)+DSTRAN(2))
            elseif(iterb.eq.2) then
               sig3o     = STRESS(3)
               eps3o     = DSTRAN(3)
               DSTRAN(3) = -(DSTRAN(1)+DSTRAN(2))
            elseif(abs(stress(3)-sig3o).gt.0.0)then
               deps      = -(DSTRAN(3)-eps3o)/(STRESS(3)-sig3o)
     +                      *STRESS(3)
               sig3o     = STRESS(3)
               eps3o     = DSTRAN(3)
               DSTRAN(3) = DSTRAN(3)+deps
            else
               deps      = 0.5*(DSTRAN(3)+eps3o)
               eps3o     = DSTRAN(3)
               sig3o     = STRESS(3)
               DSTRAN(3) = deps
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
            STRESS(3) = 0.0
            call UMAT_MODEL(DDSDDE,STRESS,STATEV,DSTRAN,PROPS2,
     +                      KINC,NTENS,NSTATV,NPROPS)
!-----------------------------------------------------------------------
c           CHECK PLANE STRESS CONDITIONS
!-----------------------------------------------------------------------
            iterb  = iterb+1
            stol   = (abs(STRESS(1))+abs(STRESS(2))+abs(STRESS(4)))
            tolalg = tol*stol
            if(abs(STRESS(3)).gt.tolalg)then
               if(iterb.lt.maxiter)then
                  GOTO 100
               else
                  print*,'shell iteration did not converge'
                  stop
               endif
            endif
  200 continue
!-----------------------------------------------------------------------
c           UNPACK STRESSES
!-----------------------------------------------------------------------
            STRESSNEW(i,1) = STRESS(1)
            STRESSNEW(i,2) = STRESS(2)
            STRESSNEW(i,3) = 0.0
            STRESSNEW(i,4) = STRESS(4)
            STRAININC(i,3) = DSTRAN(3)
!-----------------------------------------------------------------------
c           UNPACK HISTORY VARIABLES
!-----------------------------------------------------------------------
            do k=1,NSTATEV
               STATENEW(i,k) = STATEV(k)
            enddo
!-----------------------------------------------------------------------
c           UPDATE INTERNAL ENERGIES
!-----------------------------------------------------------------------
            enerInternNew(i) = enerInternOld(i)
     +        +((stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     +         +(stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     +         +(stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     +     +2.0*(stressOld(i,4)+stressNew(i,4))*strainInc(i,4))
     +     *0.5/DENSITY(i)
         enddo
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end