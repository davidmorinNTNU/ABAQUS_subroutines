!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine UVARM
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     + NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     + JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      CHARACTER*80 CMNAME,ORNAME
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION JMAC(*),JMATYP(*),COORD(*)
!-----Data from ABAQUS
      CHARACTER*3 FLGRAY(15)
      DIMENSION ARRAY(15),JARRAY(15)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      real*8 SIGH,SMISES
      real*8 SP1,SP2,SP3
      real*8 TRIAX,LODE
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
      TRIAX = -SIGH/SMISES
!-----------------------------------------------------------------------
!     Access principal stresses
!-----------------------------------------------------------------------
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      SP1 = ARRAY(3)
      SP2 = ARRAY(2)
      SP3 = ARRAY(1)
!-----------------------------------------------------------------------
!     Compute the Lode parameter
!-----------------------------------------------------------------------
      if(abs(SP1-SP3).gt.0.0)then
         LODE = (2.0*SP2-SP1-SP3)/(SP1-SP3)
      else
         LODE = 0.0
      endif
!-----------------------------------------------------------------------
!     Update user-defined variables
!-----------------------------------------------------------------------
      UVAR(1) = TRIAX
      UVAR(2) = LODE
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END
