!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine USDFLD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     + TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     + KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      CHARACTER*80 CMNAME,ORNAME
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3)
      DIMENSION T(3,3),TIME(2),COORD(*),JMAC(*),JMATYP(*)
!-----Data from ABAQUS
      DIMENSION ARRAY(15),JARRAY(15)
      CHARACTER*3 FLGRAY(15)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      real*8 SIGH,SMISES,TRIAX
!-----------------------------------------------------------------------
!     Access stress invariants
!-----------------------------------------------------------------------
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      SIGH   = ARRAY(3)
      SMISES = ARRAY(1)
!-----------------------------------------------------------------------
!     Compute the stress triaxiality
!-----------------------------------------------------------------------
      if(SMISES.gt.0.0)then
         TRIAX = -SIGH/SMISES
      else
         TRIAX = 0.0
      endif
!-----------------------------------------------------------------------
!     Update field variable
!-----------------------------------------------------------------------
      if(TRIAX.ge.0.0)then
         FIELD(1) = 1.05
      else
         FIELD(1) =-1.05
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END
