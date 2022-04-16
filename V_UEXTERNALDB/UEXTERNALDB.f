!########################################################################
!########################################################################
!     Subroutine UEXTERNALDB
!########################################################################
!########################################################################
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
!-----------------------------------------------------------------------
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      DIMENSION TIME(2)
      integer LOP,LRESTART
      integer KSTEP,KINC
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
!     Possible values for the LOP argument
      parameter(j_int_StartAnalysis  = 0,
     *          j_int_StartIncrement = 1,
     *          j_int_EndIncrement   = 2,
     *          j_int_EndAnalysis    = 3, 
     *          j_int_Restart        = 4,
     *          j_int_StartStep      = 5,
     *          j_int_EndStep        = 6)
!-----Declare global variables
      integer myThreadID
      integer i,k
      integer LENJOBNAME,LENOUTDIR
      character*256 JOBNAME,OUTDIR
      character*1000 line
!-----Declare global variables
      REAL*8 PROPS_REAL(24)         ! ARRAY WITH REAL PROPERTIES
      integer ID_REAL               ! ID for pointers
      parameter(ID_REAL = 1)        ! ID for pointers
      pointer(ptr_REAL, PROPS_REAL) ! pointer link
!--------------------------------------------------------------------------
!     Get threadID
!--------------------------------------------------------------------------
      myThreadID = get_thread_id()
!--------------------------------------------------------------------------
!     Start of the analysis
!--------------------------------------------------------------------------
      if(LOP.eq.j_int_StartAnalysis)then
         IF((kInc.eq.0).and.(myThreadID.eq.0))THEN
!--------------------------------------------------------------------------
!          CREATE GLOBAL ARRAYS
!--------------------------------------------------------------------------
           ptr_REAL = SMAFloatArrayCreate(ID_REAL, 24, 0.0)
!          FETCH ABAQUS JOBNAME
           call GETJOBNAME(JOBNAME,LENJOBNAME)
!          FETCH ABAQUS JOBNAME
           call GETOUTDIR(OUTDIR,LENOUTDIR)
!          Try to open jobname.k for extra-parameters
           OPEN(unit=15,file=trim(OUTDIR)//'/'//trim(JOBNAME)//'.k',
     +     STATUS='OLD',IOSTAT=ios)
!          Read MPC parameters from the jobname.k file
           if(ios.gt.0)then
              WRITE(6,*) 'No *.k file submitted'
           else
            k = 0
            do while (ios == 0)
!              Read line
               READ(15,FMT='(A)',end=77) line
!              Read properties from line
               if(line(1:2).ne.'**')then
                   READ(line,*,end=78) (PROPS_REAL(k+i),i=1,8)
                   k = k+8
               endif
   78 CONTINUE
            enddo
   77 CONTINUE
            CLOSE(unit=15)
           endif
         ENDIF
!--------------------------------------------------------------------------
!     Start of the increment
!--------------------------------------------------------------------------
      elseif(LOP.eq.j_int_StartIncrement)then
!--------------------------------------------------------------------------
!     End of the increment
!--------------------------------------------------------------------------
      elseif(LOP.eq.j_int_EndIncrement)then
!--------------------------------------------------------------------------
!     End of the Analysis
!--------------------------------------------------------------------------
      elseif(LOP.eq.j_int_EndAnalysis)then
      endif
!--------------------------------------------------------------------------
!     End of subroutine
!--------------------------------------------------------------------------
      RETURN
      END
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
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
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
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      real*8 T1oQ1,T2oQ2,T3oQ3
!-----Declare global variables
      REAL*8 PROPS_REAL(24)         ! ARRAY WITH REAL PROPERTIES
      integer ID_REAL               ! ID for pointers
      parameter(ID_REAL = 1)        ! ID for pointers
      pointer(ptr_REAL, PROPS_REAL) ! pointer link
!-----------------------------------------------------------------------
!     Beginning of subroutine
!-----------------------------------------------------------------------
c     Access array in the memory
      ptr_REAL = SMAFloatArrayAccess(ID_REAL)
c     Isotropic work-hardening
      SIGMA0 = PROPS_REAL(1)
      T1     = PROPS_REAL(2)
      Q1     = PROPS_REAL(3)
      T2     = PROPS_REAL(4)
      Q2     = PROPS_REAL(5)
      T3     = PROPS_REAL(6)
      Q3     = PROPS_REAL(7)
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
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END
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
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
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
      real*8 SIG1
      real*8 PEEQ
      real*8 damage,dp
!-----------------------------------------------------------------------
!-----Declaration material parameters
!-----------------------------------------------------------------------
      real*8 Wc
!-----Declare global variables
      REAL*8 PROPS_REAL(24)         ! ARRAY WITH REAL PROPERTIES
      integer ID_REAL               ! ID for pointers
      parameter(ID_REAL = 1)        ! ID for pointers
      pointer(ptr_REAL, PROPS_REAL) ! pointer link
!-----------------------------------------------------------------------
!     Read material properties
!-----------------------------------------------------------------------
c     Access array in the memory
      ptr_REAL = SMAFloatArrayAccess(ID_REAL)
c     Fracture parameter
      Wc = PROPS_REAL(9)
!-----------------------------------------------------------------------
!     Access principal stresses
!-----------------------------------------------------------------------
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      SIG1 = ARRAY(3)
!-----------------------------------------------------------------------
!     Access equivalent plastic strain
!-----------------------------------------------------------------------
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,
     +            JMAC,JMATYP,MATLAYO,LACCFLA)
c
      PEEQ = ARRAY(7)
!-----------------------------------------------------------------------
!     Extract data
!-----------------------------------------------------------------------
      dp     = PEEQ-STATEV(1)
      damage = STATEV(2)
!-----------------------------------------------------------------------
!     Update damage variable
!-----------------------------------------------------------------------
      if(dp.gt.0.0)then
         damage = damage+max(0.0,SIG1)*(dp/Wc)
      endif
!-----------------------------------------------------------------------
!     Update state dependent variables
!-----------------------------------------------------------------------
      STATEV(1) = PEEQ
      STATEV(2) = min(1.0,damage)
!-----------------------------------------------------------------------
!     Check for fracture
!-----------------------------------------------------------------------
      if(damage.ge.1.0)then
         STATEV(3) = 0.0
      else
         STATEV(3) = 1.0
      endif
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      RETURN
      END