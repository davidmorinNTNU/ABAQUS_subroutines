!########################################################################
!########################################################################
!     Subroutine VEXTERNALDB
!########################################################################
!########################################################################
      SUBROUTINE VEXTERNALDB(LOP,I_ARRAY,NIARRAY,R_ARRAY,NRARRAY)
      INCLUDE 'vaba_param.inc'
!-----------------------------------------------------------------------
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      dimension i_Array(niArray)
      dimension r_Array(nrArray)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
!     Contents of i_Array
      parameter(i_int_nTotalNodes    = 1,
     *          i_int_nTotalElements = 2,
     *          i_int_kStep          = 3,
     *          i_int_kInc           = 4,
     *          i_int_iStatus        = 5,
     *          i_int_lWriteRestart  = 6)
!     Possible values for i_Array(i_int_iStatus)
      parameter(j_int_Continue          = 0,
     *          j_int_TerminateStep     = 1,
     *          j_int_TerminateAnalysis = 2)
!     Possible values for the lOp argument
      parameter(j_int_StartAnalysis  = 0,
     *          j_int_StartStep      = 1,
     *          j_int_SetupIncrement = 2,
     *          j_int_StartIncrement = 3,
     *          j_int_EndIncrement   = 4,
     *          j_int_EndStep        = 5,
     *          j_int_EndAnalysis    = 6 )
!     Contents of r_Array
      parameter(i_flt_TotalTime = 1,
     *          i_flt_StepTime  = 2,
     *          i_flt_dTime     = 3)
!-----------------------------------------------------------------------
!     Declaration of internal variables
!-----------------------------------------------------------------------
      integer kStep,kInc,kNel
      integer LENJOBNAME,LENOUTDIR
      character*256 JOBNAME,OUTDIR
      character*1000 line
!-----Declare global variables
      REAL*8 PROPS_REAL(24)         ! ARRAY WITH REAL PROPERTIES
      integer ID_REAL               ! ID for pointers
      parameter(ID_REAL = 1)        ! ID for pointers
      pointer(ptr_REAL, PROPS_REAL) ! pointer link
!-----------------------------------------------------------------------
!     Initialization
!-----------------------------------------------------------------------
      kStep = i_Array(i_int_kStep)
      kInc  = i_Array(i_int_kInc)
      kNel  = i_Array(i_int_nTotalElements)
!-----------------------------------------------------------------------
!     Start of the analysis
!-----------------------------------------------------------------------
      if(lOp.eq.j_int_StartAnalysis)then
         IF(kInc.eq.0)THEN
!--------------------------------------------------------------------------
!          CREATE GLOBAL ARRAYS
!--------------------------------------------------------------------------
           ptr_REAL = SMAFloatArrayCreate(ID_REAL, 24, 0.0)
!          FETCH ABAQUS JOBNAME
           call VGETJOBNAME(JOBNAME,LENJOBNAME)
!          FETCH ABAQUS JOBNAME
           call VGETOUTDIR(OUTDIR,LENOUTDIR)
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
!-----------------------------------------------------------------------
!     Setup of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_SetupIncrement)then
!-----------------------------------------------------------------------
!     Start of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_StartIncrement)then
!-----------------------------------------------------------------------
!     End of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_EndIncrement)then
!-----------------------------------------------------------------------
!     End of the analysis
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_EndAnalysis)then
      endif
!-----------------------------------------------------------------------
!     End of the subroutine
!-----------------------------------------------------------------------
      return
      end
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
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
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
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 T1oQ1,T2oQ2,T3oQ3
!-----Declare global variables
      REAL*8 PROPS_REAL(24)         ! ARRAY WITH REAL PROPERTIES
      integer ID_REAL               ! ID for pointers
      parameter(ID_REAL = 1)        ! ID for pointers
      pointer(ptr_REAL, PROPS_REAL) ! pointer link
!-----------------------------------------------------------------------
!     Read material properties
!-----------------------------------------------------------------------
      if((steptime.eq.totaltime).and.(steptime.eq.zero))then
         SIGMA0 = 1e12
         T1     = 0.0
         Q1     = 0.0
         T2     = 0.0
         Q2     = 0.0
         T3     = 0.0
         Q3     = 0.0
      else
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
      endif
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
!     End of subroutine
!-----------------------------------------------------------------------
      return
      end
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
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
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