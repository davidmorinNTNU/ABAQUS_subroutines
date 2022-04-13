      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     + NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     + JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      real*8 SIGH,SMISES
      real*8 SP1,SP2,SP3
!-----------------------------------------------------------------------
!     Compute the stress triaxiality
!-----------------------------------------------------------------------
c     Access stress invariants
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      SIGH    = ARRAY(3)
      SMISES  = ARRAY(1)
c
      UVAR(1) = -ARRAY(3)/ARRAY(1)
!-----------------------------------------------------------------------
!     Compute the Lode parameter
!-----------------------------------------------------------------------
c     Access principal stresses
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      SP1 = ARRAY(3)
      SP2 = ARRAY(2)
      SP3 = ARRAY(1)
c
      if(abs(SP1-SP3).gt.0.0)then
         UVAR(2) = (2.0*SP2-SP1-SP3)/(SP1-SP3)
      else
         UVAR(2) = 0.0
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END
