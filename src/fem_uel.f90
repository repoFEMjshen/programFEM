
 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! THE CHANGES MADE ARE
 ! 1. THE ERROR CALCULATION RELATED PARAMETERS ARE DELETED
 ! 2. MAKEFILE IS MODIFIED AS ERROR.F IS DELETED
 ! 3. AS THE VARIABLES ARE NOW BOLD-FACED, THE INPUT AND HEADER FILES
 !    SHOULD ALSO BE BOLD-FACED
 ! 4. OUTPUT FOR FP IS ADDED
 ! 5. SOME USELESS ARRAYS AND VARIABLES IN ELST1 ARE REMOVED
 ! 6. SIZE OF ARRAYS IN ELST1 ARE CHECKED FOR CONSISTENCY. SIZE OF TWO
 !    ARRAYS, INCLUDING PROPS AND SHP_ST, ARE MODIFIED.
 ! 7. THE TRUE STRAIN FORMULATION IS IMPLEMENTED
 ! 8. IN UMAT, NOW DELTA_GAMMA IS STORED IN STATEV INSTEAD OF TOTAL GAMMA.
 !     YOU NEED TO STORE TOTAL GAMMA ONLY IF YOU HAVE TRANSFORMED BETA PHASE
 ! 9. OUTPUT FOR DELTA_GAMMA IS ADDED
 ! 10. IN UMAT, THE CONVENTION FOR EULER ANGLES IS MODIFIED BASED ON KALIDINDI'S PAPER
 ! 11. IN UMAT, THE LATENT HARDENING PARAMETER IS CHANGED TO 1.4
 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! deniz - changes made since Sept. 2013 --
 !
 ! 1. program options are introduced. read from options.inp
 ! 2. element - grain association is read from elementGrains.inp
 !    (OR) element list for the grains are read from grains.inp
 ! 3. node lists are created for grains, grain boundaries, and grain-grain boundaries:
 !       identifyGrainNodes
 !       identifyBoundaryNodes
 ! 4. anisotropic thermal expansion effects
 ! 5. imports program options
 ! 6. thermally activated flow rule and thermal expansion
 ! 7. improved initial trial displacement, VTAU0. Jiahao and Deniz's approaches.
 ! 8. added line search and BFGS
 ! 9. LOADVC and ELST now takes VT, instead of DVTAU as input
 ! (Dec 2014) 
 ! *  instead of inserting 1 on the diagonal, the original coeff. is retained,
 !    while setting rest of the row to zero. See modified REMOVE_DISPDOF.
 ! *  new force vectors introduced: F_ext & F_rxn (see end of loadvc), FE_appl (set in ELST1) F_appl (assembled in loadvc)
 !    F_ext = F_rxn + F_appl will be used to determine the convergence tolerance
 ! *  switchted to Gfortran 4.6, fixed memory corruptions
 ! *  switched to the new SuperLU_DIST library (4.0), removed the multiple mid-code SuperLU init/destruct.s.
 ! *  removed ICYC. using loadType variable for arbitrary user input for thermo-mechanical loading
 ! *  IROWPTR, IBINDX indexing structure changed to conventional form
 ! (Jan 2015)
 ! *  fixed bug in Loadvc. VAL is corrupted if umat error occurs, may be uncaught.
 ! *  fem_uel solver re-designed. solver flags and cutback triggers added.
 ! *  developed the task_Timer module, for CPU time allocation analysis
 ! *  adaptive time stepping re-arranged
 ! *  added error/info output argument umat_Err to the UMAT subroutine.
 ! *  added macroscopic-stress-controlled displacement BCs
 ! (Mar 2015)
 ! *  Mar 22 - fixed errors in UMAT about rotation matrices
 ! (Apr 2015)
 ! *  now using material modules
 ! (Jun 2015)
 ! * stabilized UMAT's back-stress evolution routine for cases of g --> g_sat
 ! * added support for saving and continuing a solution
 
      PROGRAM DISLOC

      use superlu_mod
      use meshtools
      use task_Timing
      use options
      use grains
      use material_Ti6242
      use nodalRecovery
      use fBarPatch

      IMPLICIT NONE

      INCLUDE 'PARDIS.H'
      INCLUDE 'mpif.h'
      
      logical, parameter :: Ti6242 = .TRUE.
      
      PARAMETER(MAXNRHS = 1)

!*********WITHOUT CONSTRAINT************************************
      REAL(8):: G0XYZ(MAXCRD),VBC1(MBOUND1),VBC3(MBOUND3)
      REAL(8):: VBCS1(MBOUND1),VBCS2(MBOUND2),XUND(MAXCRD)
      REAL(8):: SVARS1(MAXSV),SVARS2(MAXSV),TRAC(MDOFX)
      REAL(8):: stateVarsElem(MAXSTATEV)   ! a temporary local array for state variables of a single element
      REAL(8):: PROMAT2D(MAXPROPS,MAXGRP),TIME1(2)
      REAL(8):: VBC2(MBOUND2),GPXYZ(MAXCRD)
      real(8):: props(MAXPROPS)
      ! fon continuation of solution
      integer :: nElem_Continue, nNodes_Continue, nSvars_Continue
      real(8) :: time_Continue

      integer:: matID,nprops
      ! deniz
      REAL(8), allocatable :: F(:),F0(:),F_ext(:),F_uBC(:),F_BFGS(:)
      ! 	F:	      the whole (reduced) FEM force residual vector. i.e. out of balance FEM forces and RHS of the NR iterations K*ddU=F
      ! 	F_ext:	the whole (reduced) FEM external force vector. includes both support reaction forces and applied external forces.
      !  (note that F (oob forces) = F_free (directly out of UMAT for given u for a free structure) - F_ext
      !         and F_ext = F_rxn + F_appl, where F_rxn is assigned the forces of displ. constrained DOFs from F_free.
      !  (F_ext will be used in determining the convergence tolerance, see Ted Belitscko - Nonline FEM, section on convergence)

      INTEGER:: IMID(3*MAXEL),NBC1(MBOUND1)
      INTEGER:: NBC2(MBOUND2),IDF(MBOUND1)
      INTEGER:: DIR_U(MBOUND1) ! deniz - stores the directions +1/-1 of the displacement imposed by the disp() function. see input flag IFUNC
      INTEGER:: IJK(MNELX),NBC3(MBOUND3)
      INTEGER:: IBCT1(MDOFX),NPROPSG(MAXGRP)
      INTEGER:: IBC1(MDOFX),IBCT(MDOFX)
!      INTEGER:: IBC(MDOFX) -- not used

      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      INTEGER:: IBELEM(MAXEL),ILF(MAXEL),IFACE_TET(12)
      INTEGER:: IFACE_BRICK(24)
      INTEGER:: NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
      INTEGER::NX,NELX,NEQ,NPLANE,NODE
      INTEGER::NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,MCRD1
      INTEGER::NELST,NELEND,NELNOX
      ! added by deniz
      ! each processor stores the  range of elements handled by each processor
      INTEGER::NELST_proc(0:MAXPROC-1),NELEND_proc(0:MAXPROC-1), & 
                                   NELNOX_proc(0:MAXPROC-1)
      integer :: tNELST, tNELEND, tNELNOX, tSIZE
      ! ----------------------------------------
      INTEGER::IERROR, ISTATUS(MPI_STATUS_SIZE)
      INTEGER::IPROC,NPROW,NPCOL,MPIINFSIZE
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      INTEGER::NPROCS,ICYC,LEN,NGAUSS
      INTEGER::NDF,MDIM, NOSTP,NGROUP
      INTEGER::NSTATE,I1,I2,IA,IB
      INTEGER::NSEG,I,J,ICUT
      integer :: iCorner
      INTEGER:: N,NEQ1,ICOMP,JJJ

      INTEGER::II,N1,N_UPDATE,IN_ST,NNZ_LOC,NSTEP
      INTEGER::IFLAG_OVER,IWRITE_FLAG
      INTEGER::MAXREF,NITER_REF
      INTEGER::ITEREF,NUMUPD
      INTEGER:: KINC,IFLAG
      integer(4) :: IAM
      INTEGER:: IGA_ID,M
      INTEGER:: IS,JS,IE
      INTEGER:: NRHS,N_UPDATE_LOC
      INTEGER:: III,NMPD1,INIT,NROW,NCOL
      INTEGER:: NITER,NPOST1,JJ,KK
      INTEGER:: NSHIFT,ITAG,ISTV_START
      INTEGER:: ISTV_START1,ISTV_START2,ISLIP,IST1
      INTEGER:: IKIN_STATE,ISHIFT
      INTEGER:: INDX2,INDX,INFO,KSTEP

      REAL(8):: PRESSB(MAXEL),GTMP(MAXEL/MINPROC,MAXNP,NSLIP)
      REAL(8):: SOLVEINITIALIZED
      REAL(8):: TOTME,XIDT,XMXDT,XMINDT
      REAL(8):: T_NODWELL,T_DWELL
      REAL(8):: SAMP,SMIN,UAMP,UMIN,TSEG,TRAMP
      REAL(8):: TPERIOD,SLOAD,ST_RATE
      REAL(8):: DELTME,TMP,FINTME,STATME,lastConvergedTime
      REAL(8):: PNEWDT
      REAL(8):: T1,T2
      REAL(8):: T_ELAPSED,task_time
      REAL(8):: STEP
      REAL(8):: XMAXV,W2,W1,BERR(MAXNRHS)
      REAL(8):: PNEWDT1
 
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC  ! deniz - modified - NEQ_uBC: # of DOFs with Driclet BCs
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/STATEUPDATE/SVARS2
      COMMON/ABQ/TIME1,NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,SVARS1,MCRD1
      COMMON/PRESS/PRESSB,IBELEM,ILF,IFACE_TET,IFACE_BRICK
      COMMON/EL_PAR/NELST,NELEND,NELNOX
! ******************************************************************
      ! deniz - THERMAL LOADING
      REAL(8):: TEMP(2), TEMP_0
      COMMON/THERM/TEMP,TEMP_0
         
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT

! ******************************************************************
      integer,allocatable:: IBINDX(:)        ! Fortran type / 1-indexed. 
      integer,allocatable:: IROWPTR(:)       ! Fortran type / 1-indexed
      integer,allocatable:: IBINDX_LOC(:)    ! C-type / 0-indexed. Pass to SLU. May be overwritten.
      integer,allocatable:: IROWPTR_LOC(:)   ! C-type / 0-indexed. Pass to SLU. May be overwritten.
      real(8),allocatable:: RESID(:)
      real(8),allocatable:: VAL(:)
      real(8),allocatable:: RESID_LOC(:)
      real(8),allocatable:: DVTAU0(:),DVTAU(:)
      real(8),allocatable:: incVTAU_init(:),VTAU(:),incVTAU(:),VTAU0(:),VT(:)
      ! deniz - settlement stiffness matrix:
      real(8), allocatable :: VAL_rxn(:)
      integer, allocatable :: ICOLPTR_rxn(:),IBINDX_rxn(:)
      integer, allocatable :: EQ_uBC_proc(:)
      real(8), allocatable :: uBC_proc(:)
      integer :: ISIZE_rxn, N_UPDATE_rxn
      
      ! for BFGS update
      real(8),allocatable:: F_prev(:)
      real(8),allocatable:: VAR(:), WAR(:)
      real(8),allocatable:: V(:,:), W(:,:)
      integer            :: N_BFGS        ! # of iteration information stored in V and W

      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV
      COMMON/NODECONNEC/IELEMNO,NINDX
      COMMON/LOAD_DWELL/TPERIOD,TRAMP,T_NODWELL,T_DWELL,SAMP,SMIN,UAMP,UMIN
      COMMON/LOAD_MONO/SLOAD,ST_RATE
      COMMON/ICYC_FLAG/ICYC

      !old FNAME size: 10
      !modified by deniz: new FNAME size: 30
      CHARACTER(len=32) :: FNAME
      CHARACTER(len=50) :: FILNAME
      allocatable::TSEG(:)
      
      real(8),parameter :: PI=3.14159265359D0
      
   ! ******** added by deniz *********************
      ! thermal loading information      
      integer :: temp_loading_type ! 1: temperature time history given
                                   ! 2: cyclic thermal loading around initial temp
                                   ! 3: cyclic thermal loading ABOVE initial temp
      real(8) :: temp_stress_free        ! initial/stress-free temperature
      real(8) :: temp_cyclic_amplitude   ! in Kelvins
      real(8) :: temp_cyclic_period      ! in seconds
      real(8) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                         ! at integration time steps
      integer :: temp_history_N           ! number of temp-time points
 
      integer :: loadType      ! 1: Creep, 2: CSR, 
                               ! 3: stress controlled cyclic
                               ! 4: displacement controlled cyclic
                               ! 5: load-history provided
                               ! 6: displ.-history provided
                               ! 7: stress-history provided
                               ! 8: strain-history provided
      integer :: dummy1,dummy2
      integer*4 :: padding
      real(8) :: P_cyclic_max  ! in MPa
      real(8) :: P_cyclic_min  ! in MPa
      real(8) :: P_cyclic_period ! in seconds
      integer :: P_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: P_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: P_history_repeat
      real(8) :: P_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: P_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: P_history_N    ! number of load-time points
      
      real(8) :: U_cyclic_max  !
      real(8) :: U_cyclic_min  !
      real(8) :: U_cyclic_period ! in seconds
      integer :: U_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: U_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: U_history_repeat    ! 1: repeat provided u-time history, 0: keep u at the last provided value
      real(8) :: U_history_values(MAX_HISTORY_N) ! provide u-time points. 
      real(8) :: U_history_times(MAX_HISTORY_N)  ! program will interpolate at integration time steps
      integer :: U_history_N    ! number of U-time points
      
      integer :: load_axis,iStressLateral,iStressMain
      real(8) :: n_loadaxis(3),phi_load,theta_load,euler_load(3),eulerInv_load(3)
      real(8) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6) 
      logical :: specifiedLoadFrame
      logical :: strainControlledSimulation
      real(8) :: incrVonMises_imposed
      real(8) :: biaxial_StressRatioAngle      
      real(8) :: biaxial_StrainRate
      real(8) :: deformMultiplierTensionCompression
      character(len=100) :: strStrainIterMessage, strExtForces
      real(8) :: macroStrain_t(6),macroStrain_tau(6),macroF_tau(3,3)
      real(8) :: DeltaMacroStrain(6),DeltaMacroStrain_initguess(6),DeltaMacroStrain_prev(6),DeltaDeltaMacroStrain(6)
      real(8) :: stressResidual(6),stressReaction,errorVM,lambda,lambda_temp,DeltaMacroStrain_temp(6)
      
      integer :: strain_history_repeat    ! 1: repeat provided strain-time history, 0: keep u at the last provided value
      real(8) :: strain_history_values(6,MAX_HISTORY_N) ! provide strain-time points. 
      real(8) :: strain_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: strain_history_N    ! number of strain-time points

      real(8) :: xCorner_init(3,8)  ! initial positions of corner nodes
      real(8) :: ksiCorner_init(3,8)  ! initial positions of corner nodes in natural coordinates
      real(8) :: uCorner_history_values(3,8,MAX_HISTORY_N) ! provide displacement-time histories of corners of the cubic domain
      real(8) :: uCorner_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: uCorner_history_N           ! number of displacement-time points
      integer :: uCorner_remove_gradient     ! remove gradient-modes (trapezoids) from the deformation

      integer :: mark_for_export(MAX_HISTORY_N) ! mark time step for data-export. 0:default behavior, 1: force simple export, 2: force export all.
      integer :: export_special                 ! mark the special time points of cyclic loadings for data export

      common/load_cond/loadType,P_history_N,U_history_N,strain_history_N,uCorner_history_N, & 
                       P_history_repeat,U_history_repeat,P_dwell_load,U_dwell_load, &
                       P_dwell_ramptime,U_dwell_ramptime,P_cyclic_max, & 
                       P_cyclic_min,P_cyclic_period,P_history_values, & 
                       P_history_times,U_cyclic_max,U_cyclic_min, & 
                       U_cyclic_period,U_history_values,U_history_times, &
                       macroStrain_t, macroStrain_tau, load_axis,strainControlledSimulation, &
                       biaxial_StressRatioAngle,biaxial_StrainRate,deformMultiplierTensionCompression, &
                       iStressLateral,iStressMain,specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv, &
                       strain_history_values,strain_history_times,strain_history_repeat, &
                       uCorner_history_values,uCorner_history_times,uCorner_remove_gradient, &
                       ksiCorner_init,mark_for_export,export_special
      
      real(8) :: uScale(3),xCorner_range(3)
      real(8) :: dummyTime
      integer :: NFACES, NEDGES
   ! export variables for domain side lengths (for testing the thermal expansion code)
      real(8) :: domainRange0(3,2), domainLen0(3)
      real(8) :: domainRange(3,2), domainLen(3)
      COMMON/DOMAINGEO/domainRange0,domainLen0,domainRange,domainLen
   ! export variables: total plastic slip
      real(8) :: tot_gamma(100)
      real(8) :: max_gamma,norm_tot_gamma
   ! initialize arrays for nodal stress recovery
      integer, parameter :: nDOF_FP = 9
      integer, parameter :: nDOF_Nye = 9
      real(8), allocatable:: gaussValues(:,:)
      real(8), allocatable:: nodalValues(:,:)
      real(8), allocatable:: elemNye(:,:)
      real(8), allocatable:: procNye(:)
      real(8), allocatable:: nodalNye(:,:)
      real(8) :: elemNodalFp_T(3,3,4)
      real(8) :: elemNodalNyeArr(nDOF_FP,4)
      real(8) :: elemNodalNye(3,3,4)
      real(8) :: gradFp_T(3,3,3), gradCurlFp(3,3,3)
      real(8) :: Nye(3,3)
      !real(8), allocatable:: errorIndicators(:)
      !real(8), allocatable:: errorIndicatorsFp(:)
      !integer, allocatable:: elementsToRefine(:)
      real(8) :: maxError
      real(8) :: minError
      real(8) :: frobNorm
      integer :: staGrainIdx,endGrainIdx,nGrainIdx
      integer :: nFailedPatches, error

      ! global state variable array -- reduce(SVARS1)
      real(8), allocatable:: SVARS_all(:)
   !  variables for improving initial trial displacement field, VTAU0
      real(8), allocatable :: dU_BCs(:),U_BCs(:),U_BCs_prev(:)
      real(8) :: DELTME_PREV(2), time_prev(2)  ! these store the previous two: time increments and times
      real(8) :: DELTME_lastBeforeSpecTime     ! stores the last time-increment before applying a special-time cut-back
      real(8), allocatable :: u_prev(:,:)      ! the previous two nodal displacements.. u_prev(:,1) u_prev(:,2) vt(:) vtau(:)
      real(8) :: TT1, TT2, TT3
      
   ! arrays for calculating total and creep strain rates in the loading direction
      real(8), allocatable:: prevStrain(:)
      real(8), allocatable:: prevCreepStrain(:)
   
   ! initialize some matrices&vectors for testing/debugging
      real(8), allocatable:: matA(:,:)
      real(8), allocatable:: vecb(:)
      integer:: matN
      integer:: err
      
   ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)
      integer, allocatable:: elemNonLocalList(:,:)
      integer, allocatable:: elemNonLocalCount(:)
      real(8), allocatable:: elemNonLocalWeights(:,:)
      real(8) :: damageCharacteristicLength
      real(8), allocatable:: procNonLocalDamage(:)    ! temporary variable to store non-local damage variables to be sent to a processor
      real(8) :: weight, localDamage(4)

      ! integer, allocatable:: elemEdgeNeighbors(:)      
      
      character(len=50) :: strMaterialFile
      
   ! line Search
      real(8) :: G0, G
      real(8) :: G_history(2,100)
      integer :: G_n
      ! also see in program options:
      ! integer :: WLF_c2
      ! real(8) :: STEP_max
      
      include 'solver_stages.inc'
 
   ! deniz - umat error flag
      integer :: umat_Err

   ! solver flags:
      logical :: solver_CalledBefore
      logical :: solver_LUisUpToDate
      logical :: solver_reformedK_thisInc
      logical :: solver_Converged
      logical :: solverStrain_Converged
   ! triggers:
      ! error triggers
      logical :: listLocalUmatErrorsOccurred(7)
      logical :: trigger_umatError
      logical :: trigger_BroyDidNotConv
      logical :: trigger_cutBack
      logical :: trigger_reformK
      COMMON/triggers/trigger_umatError,listLocalUmatErrorsOccurred, &
                      trigger_BroyDidNotConv,trigger_cutBack,trigger_reformK
      
   ! some dummy variables
      logical :: calcStiffnessMatrix
      integer :: iElem,jElem,elemID,elemIdx,elemIdx_sta,elemIdx_end,nElem, &
                 iNode,iGloNode,iLocNode,nodeIdx,nElemCount, & 
                 iFace,iDOF,iDir
      real(8) :: facePosMax, facePosMin
      integer :: iEl_proc,tProc
      integer :: iStr,iBlock,iGrain,jGrain,grainIdx,iGrGrBnd
      integer :: Lcube
      logical :: fileExist, importedMaterial, success, successReadNeighbors, successReadNodeConnect
      integer :: printTasksTo
      real(8) :: conv_F, conv_u  ! convergence numbers for F and u. (conv <= 1: convergnece)
      real(8) :: tol_F, tol_u    ! tolerance used for F and u, calculated relative to the external forces
      logical :: tol_BelowNoise  ! external forces negligibly small
      real(8) :: F_ext_max       ! maximum converged external nodal force in the problem history. may be used in convergence criterion.
      real(8) :: lastRefactTime
      integer :: lastNITERafterRefact
      integer :: i_BFGS
      
   ! other
      real(8) :: stressMacroscopic(6),ddsddeMacroscopic(6,6)
      real(8) :: defGradMacroscopic(3,3),dpMacroscopic(6)
      real(8) :: volumeRVE,volume0RVE
      logical :: do_postProc, do_dataExport
      real(8) :: dataExport_lastT, postProc_lastT 
      logical :: EXPORT_this_step, EXPORT_next_step   ! these are flags to mark if the solution time is/will be on a special point on the load waveform
      logical :: prevTimeCutBackForSpecial            ! marks if the previous time step was cut-back to land on a speciel time on the load waveform
      character(10) :: strNSTEP
      character(10) :: strNITER
      character(10) :: strPOSQ
      real(8) :: edgeLength
      real(8) :: xNodes(3,4), xNodes0(3,4)
      real(8) :: elemJacMat(3,3)
      real(8) :: elemJacInv(3,3)
      real(8) :: elemJac, elemVol, elemVol0, totVol
      real(8) :: totStress(6),totStrain(6),totFpZZ,totCrZZ,Ep_Equiv,totWP,avgTemp
      integer, parameter :: NGRAINAVG = 7
      real(8), allocatable :: grainAVG(:,:)   ! grain-averaged quantities are stored here
      real(8), allocatable :: grainRot(:,:)   ! grain-averaged quantities are stored here
      real(8) :: R_cry(9), R_grain(9)
      real(8) :: euler(3), rot(3,3), isliptmp
      integer :: n_HCP_planes
      parameter(n_HCP_planes=8)
      integer :: nslips(n_HCP_planes), ics(n_HCP_planes)
      integer :: isys, icrys ! icrys - a pointless variable
      real(8) :: hard(47)

       ! ************* SUPER LU Structures ************* !
      INTEGER(SUPERLU_PTR) :: grid
      INTEGER(SUPERLU_PTR) :: SLUoptions
      INTEGER(SUPERLU_PTR) :: ScalePermstruct
      INTEGER(SUPERLU_PTR) :: LUstruct
      INTEGER(SUPERLU_PTR) :: SOLVEstruct
      INTEGER(SUPERLU_PTR) :: A
      INTEGER(SUPERLU_PTR) :: STAT
      
      !************* initialize MPI ************** !
      CALL MPI_INIT(IERROR)      ! here the MPI library is initialized, allowing the processes to communicate with each other.
                                 ! note that processes had been created already at the beginning of the program, and they all ran through the file-input procedures.

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)

      !*********** read in # of PROC ************ !
      OPEN(102,FILE='slu_proc_distrib.inp')
      READ(102,*)NPROW,NPCOL
      CLOSE(102)
      
      IF(NPROW*NPCOL.NE.NPROCS)THEN
         WRITE(*,*)'NPROW*NPCOL NOT EQUAL TO NPROCS'
         WRITE(*,*)'NPROW = ',NPROW
         WRITE(*,*)'NPCOL = ',NPCOL
         WRITE(*,*)'NPROCS = ',NPROCS

         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF
      
      IF(NPROCS.LT.MINPROC)THEN
         WRITE(*,*)'NPROCS LESS THAN MINPROC',NPROCS,MINPROC
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF
      
       ! ************* initialize SuperLU handles and structs ************ !
      call f_create_gridinfo_handle(grid)
      call f_create_options_handle(SLUoptions)
      call f_create_ScalePerm_handle(ScalePermstruct)
      call f_create_LUstruct_handle(LUstruct)
      call f_create_SOLVEstruct_handle(SOLVEstruct)
      call f_create_SuperMatrix_handle(A)
      call f_create_SuperLUStat_handle(stat)
      
      ! initialize state flags and triggers
      solver_CalledBefore = .FALSE.
      solver_LUisUpToDate = .FALSE.
      solver_reformedK_thisInc = .FALSE.
      solver_Converged = .FALSE.
      trigger_umatError = .FALSE.
      trigger_BroyDidNotConv = .FALSE.
      trigger_cutBack = .FALSE.
      trigger_reformK = .FALSE.

      
      !******** create superLU computing grid ******!
      !*****************and link to MPI ************!
      CALL F_SUPERLU_GRIDINIT(MPI_COMM_WORLD, NPROW, NPCOL, grid)
      CALL GET_GRIDINFO(grid, IAM=IAM)

      IF ( IAM >= NPROW * NPCOL ) THEN
         IF(IPROC.EQ.0)WRITE(*,*)'IAM>=NPROW*NPCOL'
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF


      !**** arrange PROC connections ******!
      MPIINFSIZE=0

      CALL PROC_COMB(NPROCS,MPIINF,MPIINFSIZE)

      CALL ARRANGE_PROC(IPROCINFO,ISENDORRECV,MPIINF, & 
      		MPIINFSIZE,IPROC)


!---------------------------------------------------------

      ! error handling: floating-point exceptions
      ! CALL SIGNAL(error_FPE, trap_FPE)

!--------------------------------------------------------

      OPEN(201,FILE='loadtime.inp')

      CALL READSTR(201)
      READ(201,*)ICYC   ! ICYC:  specify whether cyclic (1) or creep/constant strain rate(2)

      CALL READSTR(201)
      READ(201,*)FNAME

      CLOSE(201)
      
      PRESSB=0.D0
      IBELEM=0
      ILF=0
      
      ! added by deniz :
      ! ---------------------------------------
      ! IMPORT PROGRAM options from options.inp
      ! ---------------------------------------
      
      ! read options
      call readOptions('options.inp')
      ! print a report on the screen
      if(IPROC.EQ.0) CALL summarizeOptions()
      
! OPEN OUTPUT FILES
!-------------------------------------------

      LEN=LEN_TRIM(FNAME)
      FILNAME=FNAME(1:LEN)//'.inp'

      OPEN(UNIT=LR,FILE=FILNAME)

      IF(IPROC.EQ.0)THEN
         FILNAME=FNAME(1:LEN)//'.out'
         OPEN(UNIT=131,FILE=FILNAME)
         OPEN(UNIT=100,FILE='disp.out')          ! update period: NPOST # of steps
         
         if(DEBUG_linesearch) OPEN(UNIT=995,FILE='debug_linesearch.txt')
         if(DEBUG_timer) OPEN(994,FILE='debug_cputime.txt')
         if(DEBUG_general) OPEN(990,FILE='debug.txt')
         if(DEBUG_RESID) OPEN(970,FILE='debug_resid.txt')

         if(DEBUG_nodalValues) OPEN(997,FILE='debug_nodalValues.txt')
         if(DEBUG_nodalValues) OPEN(987,FILE='nodalValues.txt')
         OPEN(981,FILE='stress_strain_RVEavg.out')
         OPEN(405,FILE='defGrad_RVEavg.out')
         OPEN(410,FILE='dp_RVEavg.out')
         write(981,*) 'time rveVolume stress strain PlasticStrainZZ PlasticWork'
         write(405,*) 'time F11 F12 F13 F21 F22 F23 F31 F32 F33'
         write(410,*) 'time dp_xx dp_yy dp_zz dp_xy dp_xz dp_yz'
         
         if(EXPORT_watmus) then
            call get_subscript(IPROC,fname)
            filname='gfs'//fname
            open(500,file=filname)
            filname='fpfs'//fname
            open(501,file=filname)
            filname='chifs'//fname
            open(502,file=filname)
            filname='dgdtfs'//fname
            open(503,file=filname)
            filname='dfpdtfs'//fname
            open(504,file=filname)
            filname='dchidtfs'//fname
            open(505,file=filname)
            filname='gammafs'//fname
            open(506,file=filname)
            filname='dgammadtfs'//fname
            open(507,file=filname)
         endif
         if(EXPORT_temperature.OR.EXPORT_allAtLast.GT.0.0) OPEN(801,FILE='temperature.out')
         if(EXPORT_cauchy     .OR.EXPORT_allAtLast.GT.0.0) OPEN(982,FILE='cauchy.out')
         if(EXPORT_delta_gamma) OPEN(983,FILE='max_slip_wp.out')
         if(EXPORT_damage     .OR.EXPORT_allAtLast.GT.0.0) OPEN(960,FILE='damage.out')
         if(EXPORT_wp         .OR.EXPORT_allAtLast.GT.0.0) OPEN(961,FILE='wp.out')
         if(EXPORT_backstress) OPEN(962,FILE='backstress.out')
         if(EXPORT_grainAVG   .OR.EXPORT_allAtLast.GT.0.0) OPEN(984,FILE='grain_AVG.out')
         if(EXPORT_grainRot   .OR.EXPORT_allAtLast.GT.0.0) OPEN(980,FILE='grain_Rot.out')
         if(EXPORT_hardening  ) OPEN(986,FILE='g_alpha.out')
         if(EXPORT_nye        .OR.EXPORT_allAtLast.GT.0.0) OPEN(969,FILE='nye.out')
         if(EXPORT_nye        .OR.EXPORT_allAtLast.GT.0.0) OPEN(968,FILE='nye_norm.out')
         if(EXPORT_nodePos    ) OPEN(992,FILE='nodePositions.txt')
         
         if(EXPORT_delta_gamma) OPEN(404,FILE='delta_gamma.out')
         if(EXPORT_fp         .OR.EXPORT_allAtLast.GT.0.0) OPEN(406,FILE='fp.out')
         if(EXPORT_creep      .OR.EXPORT_allAtLast.GT.0.0) OPEN(407,FILE='creep.out')
         if(EXPORT_creepRate  ) OPEN(408,FILE='creepRate.out')
         if(EXPORT_strainRate ) OPEN(409,FILE='strainRate.out')

         !OPEN(UNIT=802,FILE='MOBtoGNDhard.txt') ! ratio of contributions to hardening from SSDs and GNDs 
      ENDIF
      
      if(EXPORT_stiffnessMatrix) then
         write(FILNAME,'(A,I0,A)') 'stiffnessMatrix_compact_',IPROC,'.dat'
         open(971,file=FILNAME)
      endif
      
      if(DEBUG_RESID) then
         write(FILNAME,'(A,I0,A)') 'debug',IPROC,'.txt'
         open(542,file=FILNAME)
      endif
      
      if(DEBUG_UMAT) then
         write(FILNAME,'(A,I0,A)') 'debug-UMAT',IPROC,'.txt'
         open(543,file=FILNAME)
      endif
      
      CALL GET_SUBSCRIPT(IPROC,FNAME)

      FILNAME='time.out'
      OPEN(132,FILE=FILNAME)
!------------------------------------------------------


      CALL READGM(G0XYZ,IMID,IJK,MDIM,NFACES,NEDGES)
      ! GXYZ, G0XYZ, GpXYZ:
      ! array of node positions. stacked: x1 y1 z1 x2 y2 z2 x3 ...
      ! G0XYZ: initial, GpXYZ: current.
      
      !error check - deniz
      if(NODE.GT.MAXELN) then
         if(IPROC.EQ.0) then
            write(*,*) 'Increase MAXELN to ',NODE
            WRITE(131,*)'Increase MAXELN to ',NODE
            CALL MPI_FINALIZE(IERROR)
            stop
         endif
      endif

      IF(IPROC.EQ.0)THEN
      WRITE(*,*)'#DOF #NODES #ELEMENTS', NDF ,NX, NELX
      ENDIF

      CALL READNSET()
      CALL READELSET()
      CALL TMSTEP(TOTME,XIDT,XMXDT,XMINDT,NOSTP)
      
      IF(IPROC.EQ.0)THEN
      WRITE(*,*)'TOTME XIDT XMXDT XMINDT NOSTP', & 
       TOTME, XIDT,XMXDT,XMINDT,NOSTP
      ENDIF

      CALL READMT(NGROUP,PROMAT2D,NPROPSG,NSTATE)
      
      if(Ti6242) NSTATE = mat_nstate
      
      ! determine element partitioning for all processors
      do I=0,NPROCS-1
         CALL partitionElements(NELX,NPROCS,I, & 
                                NELST_proc(I), & 
                                NELEND_proc(I), & 
                                NELNOX_proc(I))
      enddo

      NELST = NELST_proc(IPROC)
      NELEND = NELEND_proc(IPROC)
      NELNOX = NELNOX_proc(IPROC)
      
      ! allocate state variables vector
      if (IPROC==0) then
         allocate(SVARS_all(NELX*NGAUSS*NSTATE))
         SVARS_all(:) = 0.D0
      endif

      CALL INITABQ(NSTATE,MDIM)
      CALL READBC(NBC1,VBC1,IDF,DIR_U,NBC2,VBC2,NBC3,VBC3,NDF,IBC1,IBCT1)
      ! deniz - read thermal loading info
      CALL READTEMP(temp_loading_type, temp_stress_free, & 
                    temp_cyclic_amplitude, temp_cyclic_period, & 
                    temp_history_values, temp_history_times, & 
                    temp_history_N, & 
                    INFO)
      if(IPROC.EQ.0) then
         if(INFO.EQ.1) then
            write(*,*) 'Thermal loading data imported.'
         end if
      endif
            
! ********************************************************************
      ! create NODE -> element connectivity matrices
      ! IELEMNO	
      !  this array stores the node-->element connection info
      !  n1e1 n1e2 n2e1 n2e2 n2e3 n3e1 ......
      !  1         ^NINDX(1)      ^NINDX(2) .. array indices
      if (IPROC==0) write(*,'(A)',advance='no') ' identifying node-to-element connectivity... '
      CALL readNODE_ELEM_CONNEC('nodeElementConnect.dat',successReadNodeConnect)
      if (.not.successReadNodeConnect) then
         if (IPROC==0) write(*,'(A)',advance='no') 'calculating... '
         CALL NODE_ELEM_CONNEC(IJK)
         if (IPROC==0) write(*,*) 'done.'
      else
         if(IPROC==0) write(*,*) 'imported.'
      endif
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)   ! to separate READ above with WRITE below
      if(IPROC==0.and.(.not.successReadNodeConnect)) then
         CALL writeNODE_ELEM_CONNEC('nodeElementConnect.dat',IELEMNO,NINDX, &
                                                         NELX,NX,MAXNODELEM,MAXNODE)
         if(IPROC==0) write(*,*) 'saved connectivity data file: nodeElementConnect.dat'
      endif
! ********************************************************************
      
      ! added by deniz:
      ! ---------------------- IMPORT GRAINS ------------------------!

      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      if(IPROC==0) write(*,*) 'importing grains ...'
      CALL readGrains('grains.inp',NX,NELX)
      if(grainsImported.and.IPROC==0) &
         write(*,'(A,I0,A)') ' ', nGrains,' grains imported from grains.inp'
      
      if(.NOT.grainsImported) then
         CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
         if(grainsImported.and.IPROC==0) &
            write(*,'(A,I0,A)') ' ', nGrains,' grains imported from elementGrains.inp'
      endif
      
      if(.NOT.grainsImported) then
         CALL createSingleGrain(NX,NELX)
         if(IPROC==0) write(*,'(A,I0,A)') 'single grain assumed'
      endif
      
      ! grain-averaged quantities can be stored in this array:
      allocate(grainAVG(NGRAINAVG,nGrains))
      allocate(grainRot(9,nGrains))  ! grain rotations
      
      ! CALCULATE ELEMENT-ELEMENT CONNECTIVITY
      if (IPROC==0) write(*,'(A)',advance='no') ' identifying element-element connectivity... '
      allocate(elemFaceNeighbors(NELX,NFACES))
      CALL calcFaceNeighborsTET(elemFaceNeighbors(:,:),IJK,IELEMNO,NINDX, & 
                                   NELX,NODE,MAXNODELEM,MAXNODE)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)   ! to separate READ above with WRITE below
         if (IPROC==0) write(*,*) 'all processors done.'

      ! allocate(elemEdgeNeighbors(NELX*NEDGES))
      ! CALL calcEdgeNeighborsTET(elemEdgeNeighbors(:),IJK,NELX,NODE)

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      if (IPROC==0) write(*,*) 'identifying grain-node connectivity...'
      CALL identifyGrainNodes(NX,NELX,IJK)

      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      if (IPROC==0) write(*,*) 'identifying boundary nodes...'
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange0)

      ! domainRange0: initial domain range
      ! domainRange: may evolve in time
      domainRange = domainRange0
      domainLen0(:) = domainRange0(:,2) - domainRange0(:,1)
      
      if (IPROC==0) then
         write(*,*) 'Domain range detected - X', domainRange(1,:)
         write(*,*) 'Domain range detected - Y', domainRange(2,:)
         write(*,*) 'Domain range detected - Z', domainRange(3,:)
      endif

      if (IPROC.EQ.0) write(*,*) 'grain connections analysed.'
      
      !read grain sizes: grainSize(:)
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported.and.IPROC.EQ.0) write(*,*) 'grain sizes imported'
      if(.NOT.grainSizesImported.and.IPROC==0) then
         write(*,*) 'WARNING ! no grain size information (grainSizes.inp) found. assuming default grain size of:',grainSize(1)
      endif
            
      !read grain phases: grainPhase(:)
      CALL readGrainPhases('grainPhases.inp')
      if(grainPhasesImported.and.IPROC==0) write(*,*) 'grain phases imported'
      if(.NOT.grainPhasesImported.and.IPROC==0) then
         write(*,*) 'WARNING ! no grain phase information (grainPhases.inp) found. assuming default phase of:',grainPhase(1)
      endif

      !read texture: grainTexture(:)
      CALL readTexture('grainTexture.inp')
      if(textureImported.and.IPROC.EQ.0) write(*,*) 'texture imported'
      if(.NOT.textureImported.and.IPROC==0) then
         write(*,*) 'WARNING ! no texture information (grainTexture.inp) found. assuming single crystal, euler:',grainTexture(:,1)
      endif
      
      ! nodalRecovery - initialize the module,
      ! if post-processing is enabled
      if (.not.(postProcN == 0 .and. postProcT == 0.d0 )) then
         if (recoveryMethod==rec_SPR) then ! for SPR, we need to import the SPR-patches
            CALL initialize_Patches('listPatches.dat', &
                                        NX,NELX,IJK,G0XYZ)
            if (patchesInitialized) then
               if (IPROC==0) write(*,*) 'patches imported:',nSPRPatches
            else
               if (IPROC==0) write(*,*) 'error importing listPatches.dat -- required for SPR recovery'
               CALL MPI_ABORT(IERROR)
            endif
         else
            CALL initialize_GrainNodes(NX,NELX,NODE,nGrains,nTotGrainNodes)
         endif
      endif
      
      ! partition grains on CPUs
      ! this well be needed for nodal recovery
      CALL partitionElements(nGrains,NPROCS,IPROC, & 
                             staGrainIdx, & 
                             endGrainIdx, & 
                             nGrainIdx)
      
      ! calculate or read in element-element pairs for non-local damage/regularisation
      if (IPROC==0.and.damage_nonlocal) then
         
         write(*,'(A)',advance='no') ' reading non-local element pairs for damage regularisation...' 
         allocate(elemNonLocalList(400,NELX))
         allocate(elemNonLocalCount(NELX))
         allocate(elemNonLocalWeights(400,NELX))
         
         CALL readElementNonLocalList('elementNonLocalList.dat', &
                                      elemNonLocalList,elemNonLocalCount, &
                                      elemNonLocalWeights,damageCharacteristicLength, &
                                      400,NELX,success)
                                      
         
         if (.not.success) then
            write(*,*) 'file corrupt/not found. disabling non-local damage'
            ! set default: no non local effects:
            do iElem=1,NELX
               elemNonLocalCount(iElem) = 1
               elemNonLocalList(1,iElem) = iElem
               elemNonLocalWeights(1,iElem) = 1.d0
            enddo
         else
            write(*,*) 'done. Characteristic Length (um):',damageCharacteristicLength
         endif
      endif
      
      if (damage_nonlocal) then
         ! allocate buffer for sending-receiving nonlocal/regularized damage variable
         allocate(procNonLocalDamage((NELNOX+1)*4))   
         procNonLocalDamage = 0.d0
      endif
      
      ! read material properties
      if(IPROC.EQ.0) write(*,*) 'Using Material Propeties: ',materialName
      !write(strMaterialFile,*) 'props_'//trim(materialName)//'.dat'
      strMaterialFile = 'props_'//trim(materialName)//'.dat'
      call readMaterial(props,nprops,strMaterialFile,importedMaterial)
      if(.not.importedMaterial) then
         if(IPROC.EQ.0) write(*,*) 'Can not import material parameters: props_Ti6242.dat corrupt/not found'
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         CALL MPI_FINALIZE(IERROR)
         STOP
      endif
      if(IPROC.EQ.0) write(*,*) '# properties: ',nprops
      if(IPROC.EQ.0) write(*,*) '# state variables: ',NSTATE
      
      ! initialize material parameters
      CALL initialize_Parameters(props,nprops,'Ti6242')
      ! initialize state variables
      do iElem = NELST, NELEND
         iEl_proc = iElem - NELST + 1
         CALL initialize_StateVariables(stateVarsElem,NSTATE,iElem,grainSize,grainTexture,grainPhase,nGrains,grainIDelem,NELX)
         ! copy to the processor-local state variable array:
         SVARS1((iEl_proc-1)*NGAUSS*NSTATE+1:iEl_proc*NGAUSS*NSTATE) = stateVarsElem(1:NSTATE)
      enddo

      ! added by deniz -------------------------------
      ! initialize parameters/arrays for error analysis
      ! ((for master processors))
      if (IPROC.EQ.0) then

         ! arrays below are solely used in error analysis,
         ! and therefore only allocated in the master processor
         
         allocate(gaussValues(9,NELX*NGAUSS))
         allocate(nodalValues(9,nTotGrainNodes))
         allocate(elemNye(nDOF_Nye,NELX))
         allocate(nodalNye(nDOF_Nye,nTotGrainNodes))
         elemNye(:,:)=0.D0
         nodalNye(:,:)=0.D0
!         allocate(errorIndicators(NELX))
!         errorIndicators(:) = 0.D0
         gaussValues(:,:) = 0.D0
         nodalValues(:,:) = 0.D0
         maxError = 0.D0
         minError = 0.D0

      endif
      
      ! allocate arrays for calculating local strain rates
      if (IPROC==0) then
         if(EXPORT_creepRate  ) then
            allocate(prevCreepStrain(NELX))
            prevCreepStrain = 0.d0
         endif
         if(EXPORT_strainRate ) then
            allocate(prevStrain(NELX))
            prevStrain = 0.d0
         endif
      endif
      
      ! ((for all processors))
      allocate(procNye((NELNOX+1)*9))   ! +1, in case some processors have one additional element
      procNye(:)=0.D0
      !allocate(elementsToRefine(NELX))
      !elementsToRefine(:) = 0

!--------------READS TRACTION OR DISPLACEMENT FLAG----------------------


      OPEN(201,FILE='loadtime_creep.inp')

      CALL READSTR(201)
      READ(201,*)TRAMP

      CALL READSTR(201)
      READ(201,*)loadType
      
      strainControlledSimulation = .false. ! unless it is, which will be set below depending on loadType
      specifiedLoadFrame = .false.         ! set to true for loading along arbitrary axis

      IF(loadType==1)THEN       ! CREEP LOADING
      
         CALL READSTR(201)
         ST_RATE=0.D0
         READ(201,*)SLOAD
      ELSEIF(loadType==2)then   ! CONSTANT STRAIN RATE LOADING
      
         CALL READSTR(201)
         SLOAD=0.D0     ! dload() will return 0 pressure for any pressure BCs
         READ(201,*)ST_RATE

         ! loadType > 2 added by Deniz
      ELSEIF(loadType==3)THEN   ! TRACTION CONTROLLED CYCLIC LOADING
      
         CALL READSTR(201)
         READ(201,*) P_cyclic_max, P_cyclic_min, & 
                     P_cyclic_period, P_dwell_load, export_special

         if(P_dwell_load.EQ.1) P_dwell_ramptime = TRAMP
         
         if(EXPORT_specialOnly) export_special = 1 ! if it is requested that data export is done only at special points, 
                                                   ! override this flag to turn on solution+export at those special points

      ELSEIF(loadType==4)THEN   ! DISPLACEMENT CONTROLLED CYCLIC LOADING

         CALL READSTR(201)
         READ(201,*) U_cyclic_max, U_cyclic_min, & 
                     U_cyclic_period, U_dwell_load, export_special
                     
         if(U_dwell_load.EQ.1) U_dwell_ramptime = TRAMP
         
         if(EXPORT_specialOnly) export_special = 1 ! if it is requested that data export is done only at special points, 
                                                   ! override this flag to turn on solution+export at those special points
         
      ELSEIF(loadType==5)THEN   ! EXPLICIT TRACTION-TIME HISTORY GIVEN

         CALL READSTR(201)
         READ(201,*) P_history_N, P_history_repeat
   
          if(P_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N to P_history_N:', & 
                                                     P_history_N
            CALL MPI_FINALIZE(IERROR)
            STOP
         endif
         
         P_history_times = 0.0
         P_history_values = 0.0
         do I=1,P_history_N
            READ(201,*) P_history_times(I), P_history_values(I), mark_for_export(I)
         enddo
         
      ELSEIF(loadType==6)THEN   ! EXPLICIT DISPLACEMENT-TIME HISTORY GIVEN

         CALL READSTR(201)
         READ(201,*) U_history_N, U_history_repeat
         
          if(U_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for U_history_N:', & 
                                                      U_history_N
            CALL MPI_FINALIZE(IERROR)
            STOP
         endif
         
         U_history_times = 0.0
         U_history_values = 0.0
         do I=1,U_history_N
            READ(201,*) U_history_times(I), U_history_values(I), mark_for_export(I)
         enddo
         
      ELSEIF(loadType==7)THEN   ! EXPLICIT (macroscopic)STRESS-TIME HISTORY GIVEN

         write(*,*) 'explicit stress-time history not yet supported'
         CALL MPI_FINALIZE(IERROR)
         stop

      ELSEIF(loadType==8)THEN   ! EXPLICIT (macroscopic) STRAIN-TIME HISTORY GIVEN

         strainControlledSimulation = .true.
         
         CALL READSTR(201)
         READ(201,*) strain_history_N, strain_history_repeat
         
          if(strain_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for strain_history to:', & 
                                                      strain_history_N
            CALL MPI_FINALIZE(IERROR)
            STOP
         endif
         
         strain_history_times(:) = 0.0
         strain_history_values(:,:) = 0.0
         do I=1,strain_history_N
            READ(201,*) strain_history_times(I), strain_history_values(:,I), mark_for_export(I)
         enddo
         
      ELSEIF(loadType==18)THEN   ! EXPLICIT (macroscopic) displacement-TIME HISTORY GIVEN for corners of the cubic domain. 
                                 ! linear interpolation is used to set boundary conditions

         strainControlledSimulation = .true.
         
         CALL READSTR(201)
         READ(201,*) dummyTime,(xCorner_init(:,iCorner),iCorner=1,8)
         
         CALL READSTR(201)
         READ(201,*) uCorner_history_N, uCorner_remove_gradient
         
          if(uCorner_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for strain_history to:', & 
                                                      uCorner_history_N
            CALL MPI_FINALIZE(IERROR)
            STOP
         endif
         
         xCorner_range(1) = maxval(xCorner_init(1,:))-minval(xCorner_init(1,:))
         xCorner_range(2) = maxval(xCorner_init(2,:))-minval(xCorner_init(2,:))
         xCorner_range(3) = maxval(xCorner_init(3,:))-minval(xCorner_init(3,:))
         uScale(1) = domainLen0(1)/xCorner_range(1)
         uScale(2) = domainLen0(2)/xCorner_range(2)
         uScale(3) = domainLen0(3)/xCorner_range(3)
         ksiCorner_init(1,:) = xCorner_init(1,:) - sum(xCorner_init(1,:))/8
         ksiCorner_init(2,:) = xCorner_init(2,:) - sum(xCorner_init(2,:))/8
         ksiCorner_init(3,:) = xCorner_init(3,:) - sum(xCorner_init(3,:))/8
         ksiCorner_init(1,:) = ksiCorner_init(1,:)/xCorner_range(1)*2
         ksiCorner_init(2,:) = ksiCorner_init(2,:)/xCorner_range(2)*2
         ksiCorner_init(3,:) = ksiCorner_init(3,:)/xCorner_range(3)*2
         
         uCorner_history_times(:) = 0.0
         uCorner_history_values(:,:,:) = 0.0
         do I=1,uCorner_history_N
            READ(201,*) uCorner_history_times(I), (uCorner_history_values(:,iCorner,I),iCorner=1,8), mark_for_export(I)
            ! translate displacements to keep the domain centered at origin (a rigid body translation)
            uCorner_history_values(1,:,I) = uCorner_history_values(1,:,I) - sum(uCorner_history_values(1,:,I))/8
            uCorner_history_values(2,:,I) = uCorner_history_values(2,:,I) - sum(uCorner_history_values(2,:,I))/8
            uCorner_history_values(3,:,I) = uCorner_history_values(3,:,I) - sum(uCorner_history_values(3,:,I))/8
         enddo
         ! scale displacements to account for any discrepancies in units (e.g. um/m)
         uCorner_history_values(1,:,:)=uScale(1)*uCorner_history_values(1,:,:)
         uCorner_history_values(2,:,:)=uScale(2)*uCorner_history_values(2,:,:)
         uCorner_history_values(3,:,:)=uScale(3)*uCorner_history_values(3,:,:)
         
      ELSEIF(loadType==9)THEN   ! mixed macroscopic stress-strain controlled 
      
         strainControlledSimulation = .true.
      
         CALL READSTR(201)
         read(201,*) load_axis
         CALL READSTR(201)
         READ(201,*) biaxial_StrainRate, biaxial_StressRatioAngle
                  
         ! convert to radians
         biaxial_StressRatioAngle = biaxial_StressRatioAngle/180.d0*PI
         
         deformMultiplierTensionCompression = +1.0
         
         ! trial strain
         if (load_axis==6) then ! biaxial on Y & Z
            iStressMain = 3
            iStressLateral = 2
         elseif (load_axis==5) then   ! biaxial Z & X
            iStressMain = 1
            iStressLateral = 3
         elseif (load_axis==4) then   ! biaxial X & Y
            iStressMain = 2
            iStressLateral = 1
         else
            write(*,*) 'unsupported load axis for biaxial loading:',load_axis
            stop
         endif
         
      ELSEIF(loadType==10)THEN   ! uniaxial deformation along arbitrary direction with constant VM-strain rate
      
         strainControlledSimulation = .true.
         
         CALL READSTR(201)
         read(201,*) n_loadaxis
         CALL READSTR(201)
         READ(201,*) biaxial_StrainRate, deformMultiplierTensionCompression
                  
         ! define uniaxial loading as a special case of biaxial loading (ZY) with biaxiality angle of 0deg
         biaxial_StressRatioAngle = 0.d0
         load_axis = 6
         iStressMain = 3
         iStressLateral = 2
         
         ! define load frame, calculate associated rotation matrix for the Voigt-form
         specifiedLoadFrame = .true.
         call getPhi1andTheta(n_loadaxis(1),n_loadaxis(2),n_loadaxis(3),phi_load,theta_load)
         euler_load(1) = phi_load
         euler_load(2) = theta_load
         euler_load(3) = 0.d0
         eulerInv_load(1) = -euler_load(3)
         eulerInv_load(2) = -euler_load(2)
         eulerInv_load(3) = -euler_load(1)
         call getRotationMatrixVoigt(R_loadVoigt,euler_load)
         call getRotationMatrixVoigt(R_loadVoigt_Inv,eulerInv_load)
         
      ENDIF      

      CLOSE(201)

      ! overwrite displ. BC input for strain-controlled deformation
      ! disp BCs are assigned to the boundary nodes of the domain
      if (strainControlledSimulation) &
         call setBC_StrainControlled(NBC1,VBC1,VBCS1,IDF,DIR_U,NBC2,VBC2,VBCS2,NBC3,VBC3,NDF,IBC1,IBCT1)
         
      if (IPROC.EQ.0) write(*,*) 'boundary conditions/load cases imported.'

      TPERIOD=0.D0
      T_NODWELL=0.D0
      T_DWELL=0.D0
      SAMP=0.D0
      SMIN=0.D0
      UAMP=0.D0
      UMIN=0.D0

      NSEG=0.D0


!----------------------------------------------------------------------------------

      DELTME=XIDT
      ICUT=0

      DO I=1,MCRD1*NX
         GPXYZ(I)=G0XYZ(I)
         XUND(I)=G0XYZ(I)
      ENDDO

      TRAC(1:NEQ)=0.D0

      NEQ1=NEQ
      
      ! deniz comment
      DO N=1,NBOUD1  !NBOUND1: total number of nodes*DOFs on which displ. BCs were defined
         IF(IDF(N).EQ.0)THEN  ! finite displacements provided in input file. impose displacements gradually.
            DO II=1,NDF
               N1=NDF*(N-1)+II
               VBCS1(N1)=VBC1(N1)/TOTME
               VBC1(N1)=0.D0
            ENDDO
         ENDIF
      ENDDO

      ! NBOUD2: # of defined 'point force' boundary conditions
      DO N=1,NBOUD2
         DO II=1,NDF
            N1=NDF*(N-1)+1
            VBCS2(N1)=VBC2(N1)/TOTME
         ENDDO
      ENDDO
      
! -------- initialize FBar Patches
      CALL initializePatches(NELX,NX,NDF,NODE, &
                             IPROC,NPROCS, &
                             NELST_proc,NELEND_proc,NELNOX_proc, &
                             iError)
                             
      if (iError /= 0) then
         if (IPROC==0) write(*,*) 'error initializing fBar Patches. error #:',iError
         CALL MPI_ABORT(IERROR)
      else
         if (IPROC==0) write(*,*) '# of FBar patches imported:',nPatches
      endif

      ! partitioning of ELEMENTS and EQUATIONS are done separately
      ! ELEMENTS: partitionElements() - for LoadVC parallelization
      ! EQUATIONS: NODE_PARTITION     - for SuperLU parallelization
      !OUT:
      !  IN_ST: first eqn#-1 for this processor
      !  N_UPDATE: # of eqns handled by this processor
      CALL NODE_PARTITION(N_UPDATE,IN_ST,IPROC,NPROCS,NX,NDF)

      NNZ_LOC=0

      CALL INIT_SIZE_SLU(IJK,IPROC,N_UPDATE,IN_ST,NNZ_LOC)
      CALL fBar_INIT_SIZE_SLU(IJK,IPROC,N_UPDATE,IN_ST,NNZ_LOC, & ! find non-zeroes in the stiffness matrix according to FBar couplings
                              IELEMNO,NINDX)

      ALLOCATE(RESID(N_UPDATE))
      ALLOCATE(RESID_LOC(N_UPDATE))
      ALLOCATE(VAL(NNZ_LOC))
      ALLOCATE(IROWPTR(N_UPDATE+1))
      ALLOCATE(IBINDX(NNZ_LOC))
      ALLOCATE(IBINDX_LOC(NNZ_LOC))
      ALLOCATE(IROWPTR_LOC(N_UPDATE+1))

      VAL(:)=0.D0
      RESID(:)=0.D0
      
      ! this subroutine must be re-written when periodic BCs are introduced
      ! constructs sparse matrix structure, using connectivity/coupling information
      ! output: IROWPTR, IBINDX
      CALL INIT_MSR_SLU(IJK,IPROC,N_UPDATE,IN_ST,IROWPTR, & 
            IBINDX,NNZ_LOC)
      CALL fBar_INIT_MSR_SLU(IJK,IPROC,N_UPDATE,IN_ST,IROWPTR,IBINDX,NNZ_LOC, &
                             IELEMNO,NINDX)
                                   
      ! convert IROWPTR, IBINDX (Fortran) --> IROWPTR_LOC, IBINDX_LOC (C)
      CALL convertSPARSE_f2C(IROWPTR,IBINDX,IROWPTR_LOC,IBINDX_LOC, &
                             N_UPDATE,NNZ_LOC)
                             
      !****************** initialize SUPERLU *********************!
      ! **** create SuperLU superMatrix Struct
      N=NEQ
      M=NEQ
      NRHS=1
      
      CALL f_dCreate_CompRowLoc_Mat_dist(A, M, N, NNZ_LOC, & 
      N_UPDATE,IN_ST,VAL,IBINDX_LOC,IROWPTR_LOC,   &
      SLU_NR_LOC,SLU_D,SLU_GE)

      ! make sure superLU has the correct RHS length
      CALL GET_COMPROWLOC_MATRIX(A, NROW_LOC=N_UPDATE_LOC)

      ! **** set options *** see pages 57--59 of the guide.
      CALL f_set_default_options(SLUoptions)    ! set default options. superlu guide page 60:
      ! 1) no ordering
      !CALL set_superlu_options(SLUoptions,COLPERM=NATURAL)    !** these use too much memory
      !CALL set_superlu_options(SLUoptions,ROWPERM=NOROWPERM)  !**
      ! 2) parmetis
      !call set_superlu_options(SLUoptions,ColPerm=PARMETIS)
      !call set_superlu_options(SLUoptions,ParSymbFact=YES) ! <- SymbFact causes a floating point exception ?
      ! 3) Metis A^T + A
      call set_superlu_options(SLUoptions,ColPerm=METIS_AT_PLUS_A)
      CALL set_superlu_options(SLUoptions,ROWPERM=LargeDiag)
      call set_superlu_options(SLUoptions,Equil=NO)      
      ! deniz - SymPattern - support added to SuperLU wrapper: superlu_c2f_dwrap.c and superlu_mod.f90
      !call set_superlu_options(SLUoptions,SymPattern=YES)   
      
      ! initialize ScalePermStruct, LUStruct and the statistics variables
      CALL GET_SUPERMATRIX(A,NROW=M,NCOL=N)
      CALL F_SCALEPERMSTRUCTINIT(M, N, ScalePermstruct)
      !CALL F_LUSTRUCTINIT(N, LUstruct)
      call f_LUstructInit(m, n, LUstruct) ! this was fixed in SuperLU 4.1
      CALL F_PSTATINIT(STAT)             

      CALL assertInt(NNZ_LOC,IROWPTR(N_UPDATE+1)-1,'NNZ_LOC!=IROWPTR')
      CALL assertInt(N_UPDATE,N_UPDATE_LOC,'N_UPDATE_LOC')
      !******************************************************!
                                    
      ! deniz - settlement stiffness matrix in sparse column form
      CALL INIT_SIZE_RXN(N_UPDATE_rxn,ISIZE_rxn, & 
                         IN_ST,N_UPDATE,IROWPTR,IBC1)  ! determines N_UPDATE_rxn and ISIZE_rxn
      
      allocate(VAL_rxn(ISIZE_rxn))           ! reaction stiffness matrix in sparse column format
      allocate(ICOLPTR_rxn(N_UPDATE_rxn+1))  ! data-indices of the first coefficients at each column
      allocate(IBINDX_rxn(ISIZE_rxn))        ! row/EQN IDs of coefficient k, where k is the 1-D data index in column major format
      allocate(EQ_uBC_proc(N_UPDATE_rxn))    ! list of constrained DOF/EQN numbers on this proc.
      allocate(uBC_proc(N_UPDATE_rxn))       ! imposed displacements on this proc.
      
      VAL_rxn(:) = 0.D0
      
      CALL INIT_MSC_rxn(ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn, &    ! constructs ICOLPTR_rxn and IBINDX_rxn
                        ISIZE_rxn,IN_ST,N_UPDATE,NNZ_LOC, &         
                        IROWPTR,IBINDX,IBC1)
      ! output:
      !     ICOLPTR_rxn, IBINDX_rxn       sparse matrix handlers for VAL_rxn
      !     N_UPDATE_rxn                  # of constrained DOFs defined among the EQNs of this proc.
                        
      CALL INIT_uBC_proc(EQ_uBC_proc,N_UPDATE_rxn,IBC1,IN_ST,N_UPDATE)
      !     EQ_uBC_proc(1:N_UPDATE_rxn)   list of constrained DOF/EQN numbers on this proc.
      
!*************************************************************************

      NSTEP= 1
      NITER= 0
      FINTME=0.D0          ! <---- take over FINTME (below). this will be copied into STATME right inside the solver loop

      IFLAG_OVER=0

      IWRITE_FLAG=1

      MAXREF=18
      NITER_REF=10
      
      tol_F = 1.D-6 ! min tolerance to use, in case external forces are zero
      tol_u = 1.D-6 ! min tolerance to use, in case displacements are zero

      IF(DABS(TRAMP-0.D0).LT.1D-10)IWRITE_FLAG=1

      allocate(F(MDOFX),F0(MDOFX),F_ext(MDOFX),F_uBC(MDOFX),F_BFGS(NEQ))
       
      allocate(incVTAU_init(NEQ),DVTAU0(NEQ),DVTAU(NEQ))
      allocate(VTAU(NEQ),incVTAU(NEQ),VTAU0(NEQ),VT(NEQ))
      
      ! for BFGS update
      allocate(VAR(N_UPDATE*MAXREF),WAR(N_UPDATE*MAXREF))
      allocate(V(NEQ,MAXREF+1),W(NEQ,MAXREF+1))
      allocate(F_prev(NEQ))
      
      allocate(dU_BCs(NEQ),U_BCs(NEQ),U_BCs_prev(NEQ))
      allocate(u_prev(NEQ,2))
      U_BCs_prev = 0.0
      U_BCs = 0.0
      dU_BCs = 0.0
      DELTME_PREV(:) = DELTME
      time_prev(:) = 0.D0
      DELTME_lastBeforeSpecTime = 0.D0
      postProc_lastT = 0.D0
      dataExport_lastT = 0.D0
      DVTAU = 0.D0
      DVTAU0 = 0.D0
      incVTAU_init = 0.D0
      VTAU = 0.D0
      incVTAU = 0.d0
      VTAU0 = 0.D0
      VT = 0.D0                     ! <----- init this to the last disp. field
      macroStrain_t = 0.d0
      macroStrain_tau = 0.d0
      
      !---------------------- take over state variables -- START
      if (continueSolution) then    
         if(IPROC==0) write(*,*) 'importing and distributing state variables...'
         open(102,FILE='continue/problem.out')
         read(102,*) nSvars_Continue
         read(102,*) nElem_Continue
         read(102,*) nNodes_Continue
         read(102,*) time_Continue
         close(102)
         if (NSTATE /= nSvars_Continue .or. NELX /= nElem_Continue .or. NX /= nNodes_Continue) then
            if(IPROC==0) write(*,*) 'number of state variables, or elements, or nodes in the saved solution &
                                    &do not match the current mesh.'
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
            CALL MPI_FINALIZE(IERROR)
            STOP
         endif
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         
         FINTME = time_Continue                    ! import the solution time
         
         if (IPROC==0) then                        ! read and send state variables
            open(102,FILE='continue/svars.out')
            do iElem=1,NELX
               read(102,*) SVARS_all((iElem-1)*NGAUSS*NSTATE+1:iElem*NGAUSS*NSTATE)
            enddo
            close(102)
            
            ! copy the part for the master processor
            SVARS1(1:NELNOX*NGAUSS*NSTATE) = SVARS_all((NELST-1)*NGAUSS*NSTATE+1:NELEND*NGAUSS*NSTATE)
            
            ! send to the others:
            do tProc = 1,NPROCS-1
            
               tNELST = NELST_proc(tProc)
               tNELEND = NELEND_proc(tProc)
               tNELNOX = NELNOX_proc(tProc)
               tSIZE = tNELNOX*NGAUSS*NSTATE
               
               CALL MPI_SSEND(SVARS_all((tNELST-1)*NGAUSS*NSTATE+1:tNELEND*NGAUSS*NSTATE), &
                             tSIZE,MPI_DOUBLE_PRECISION, &
                             tProc,969,MPI_COMM_WORLD,IERROR)
                             
               if(IERROR.NE.0) exit
            enddo

         else                                      ! receive state variables
         
               tSIZE = NELNOX*NGAUSS*NSTATE
               CALL MPI_RECV(SVARS1(1:tSIZE), &
                             tSIZE,MPI_DOUBLE_PRECISION, &
                             0,969,MPI_COMM_WORLD,ISTATUS,IERROR)
         endif
   
         if(IERROR.NE.0) then
            write(*,*) 'MPI error at sharing state variables for continuing the solution'
            CALL MPI_FINALIZE(IERROR)
            stop
         endif

      endif
      !---------------------- take over state variables -- END
      
      !---------------------- take over displacements at time t, VT(:)
      if (continueSolution) then    
         if(IPROC==0) write(*,*) 'importing and distributing the displacement field...'
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         
         if (IPROC==0) then                        ! read and send state variables
            open(102,FILE='continue/disp.out')
            CALL READSTR(102)
            do iNode=1,NX
               read(102,*) VT((iNode-1)*3+1:iNode*3)
            enddo
            close(102)
            
            ! send to the others:
            do tProc = 1,NPROCS-1
               
               CALL MPI_SSEND(VT(1:NEQ), &
                             NEQ,MPI_DOUBLE_PRECISION, &
                             tProc,969,MPI_COMM_WORLD,IERROR)
                             
               if(IERROR.NE.0) exit
            enddo

         else                                      ! receive state variables
         
               CALL MPI_RECV(VT(1:NEQ), &
                             NEQ,MPI_DOUBLE_PRECISION, &
                             0,969,MPI_COMM_WORLD,ISTATUS,IERROR)
         endif
   
         if(IERROR.NE.0) then
            write(*,*) 'MPI error at sharing displacement field for continuing the solution'
            CALL MPI_FINALIZE(IERROR)
            stop
         endif

      endif
      !---------------------- take over state variables -- END
      
      u_prev = 0.D0
      EXPORT_this_step = .FALSE.
      EXPORT_next_step = .FALSE.
      prevTimeCutBackForSpecial = .FALSE.
      
      ! initialize timer      
      CALL task_Timing_Initialize()
      CALL task_Start(task_PROGRAM,MPI_WTIME())
      if(IPROC.EQ.0.AND.DEBUG_timer)  then
         printTasksTo = 994    ! master proc. prints timing information to file
         write(994,*) "method for improving displ. field: ",improve_U0
         write(994,*) "0:use raw 1:Jiahaos 2:Deniz's 3:Quad-Jiahao 4:tangentModulus"
         write(994,*) "non-linear solver: ",solverMethod
         write(994,*) "1:newton-Raphson 2:modif Newt 3:modif Newt (Keep K) 4: BFGS 5: BFGS (Keep K) 6: Broyden"
         write(994,*) "linesearch method:",lschMethod
         write(994,*) "0:disabled, 1:deniz 2:matthies, strang, 1979"
         write(994,*) "eps_F, eps_u: ",eps_F,eps_u
         write(994,*) "convergenceNorm: ",convergenceNorm
         write(994,*) "1: L_max, 2: L_2 norm"
      else
         printTasksTo = -1   ! not master proc, do not print
      endif

      ! *******************************************************************************!
		! the main loop over time steps
      ! *******************************************************************************!
      
		! IFLAG_OVER = 1 time integration is complete
      DO WHILE(IFLAG_OVER.EQ.0)
      
         KINC=NSTEP
   
         if(N_UPDATE.EQ.0) then
            write(*,*) 'N_UPDATE 0, beginning of loop'
            CALL MPI_FINALIZE(IERROR)
            stop
         endif

         IF(IPROC.EQ.0)THEN
            WRITE(*,*)'START INCREMENT', NSTEP
         ENDIF
         
         !if there had been a cutback, we just had:
         !ICUT=ICUT+1
         !FINTME=FINTME-DELTME
         !DELTME=PNEWDT*DELTME

         STATME=FINTME
         FINTME=FINTME+DELTME
         TIME1(1)=FINTME
         TIME1(2)=FINTME
         
         ! mark for exporting data on this step, if this step was just marked for exporting (special time in load waveform)
         EXPORT_this_step = EXPORT_next_step
         EXPORT_next_step = .false. ! reset - done its job

         ! reset the cut-back trigger. it has done its job.
         trigger_cutBack = .FALSE.
         ! initialize the flag, reformedK_thisIter, which states whether K_inv has been reformed during this increment
         solver_reformedK_thisInc = .FALSE.
         
        IF(IPROC.EQ.0)THEN
            WRITE(*,*)'FINTME TRAMP DELTME',FINTME,TRAMP,DELTME
         ENDIF

         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'ATTEMPTING NSTEP=',NSTEP
         ENDIF

         KSTEP=1
         PNEWDT=1.D0
         
         VTAU0(1:NEQ) = VT(1:NEQ)
        
         CALL TRAINC(NBC2,VBC2,VBCS2,DELTME,NDF,TRAC)    ! TRAC is not used anywhere..
         
         ! deniz - update temperature state variable, 
         ! according to the prescribed thermal loading
         CALL TEMPINC(temp_loading_type, temp_stress_free, & 
                      temp_cyclic_amplitude, temp_cyclic_period, & 
                      temp_history_values, temp_history_times, & 
                      temp_history_N, TIME1(:), TEMP(:), TEMP_0)
         ! print thermal conditions
         IF(IPROC.EQ.0)THEN
            WRITE(*,*)'AMBIENT TEMPERATURE',TEMP(2)
         ENDIF

         ! DISINC (see below) updates the imposed displacements in VTAU0 for the new time step.
         ! above is the standard method to
         ! introduce the imposed displacements to the initial trial displacement field
         ! the free DOFs remain equal to the previous (converged) values
         ! PROBLEM WITH THIS METHOD:
         ! the jump from previous field, to the updated/prescribed displacements at the boundary nodes
         ! cause excessive initial trial strains at elements near the constrained DOFs, 
         ! causing dt cutbacks, extremely small time increments.
         ! this effect is much more pronounced in fine meshes, as (du / h) > 1 for the boundary elements
         ! and also for high strain rate loadings (large du per increment)
         !==============================================================================
         ! START :             improve initial guess
         !==============================================================================
         
         CALL setFLAGS(IFLAG,solver_LUisUpToDate,stage_initVTAU,NSTEP,&
                       solverMethod,solver_CalledBefore, &
                       solver_reformedK_thisInc, &
                       trigger_cutBack,trigger_reformK)
         
         if(DEBUG_timer.AND.IFLAG.EQ.0) CALL task_Start(task_calcResidual_K,MPI_WTIME())
         if(DEBUG_timer.AND.IFLAG.NE.0) CALL task_Start(task_calcResidual,  MPI_WTIME())
         
         calcStiffnessMatrix = (IFLAG == 0)
         call calcFBarStabilizers(VT,VT, &! displacements over cycle, time increment info
            IJK,G0XYZ,                   &                     ! mesh,parallelism etc.
            calcStiffnessMatrix,success)               
         
         CALL LOADVC(F,F_ext, &
         stressMacroscopic,ddsddeMacroscopic, &
         defGradMacroscopic,dpMacroscopic, &
         volumeRVE,volume0RVE, &
         VT,VT,XUND,IJK,IMID,NPROPSG, & 
         PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE,IN_ST,IROWPTR, & 
         IBINDX,VAL,RESID,NNZ_LOC,VAL_rxn,ICOLPTR_rxn,IBINDX_rxn, & 
         N_UPDATE_rxn,ISIZE_rxn,IFLAG,PNEWDT)
                     
         if(IFLAG.EQ.0) CALL task_Record(task_calcResidual_K,MPI_WTIME(),task_time,printTasksTo)
         if(IFLAG.NE.0) CALL task_Record(task_calcResidual  ,MPI_WTIME(),task_time,printTasksTo)
         
         if (EXPORT_stiffnessMatrix) then
            CALL printMatrix_Sparse(VAL,IBINDX,IROWPTR,  &
                                       N_UPDATE,NEQ,NNZ_LOC,   &
            1,'stiffness matrix, after first loadvc',971)
!               CALL printMatrix_Sparse_Compact(VAL,IBINDX,IROWPTR,  &
!                                       N_UPDATE,NEQ,IN_ST,IPROC,NNZ_LOC, &
!                                       1,'stiffness matrix',971)
            if (IPROC==0) write(*,*) 'EXPORT_stiffnessMatrix = .true. exported stiffness matrix. quitting..'
            flush(971)
            close(971)
            CALL MPI_FINALIZE(IERROR)
            stop
         endif
                                    
         
         F0(1:NEQ) = F(1:NEQ)    ! save initial residual for the convergence test. 
                                 ! what to do if another initial guess is used?

         if(trigger_umatError) then                      ! UMAT error. R is corrupted.
            trigger_cutBack = .TRUE.                     ! cut back time step and recalculate increment
            if(IFLAG.EQ.0) trigger_reformK = .TRUE.      ! K is corrupted too. Mark for reformation next time loadVC is called
            goto 700   ! maybe I dont need this?         ! go to the cut-back handler (right after the q-Newton loop..)
         endif
         
         ! if VM-rate constrained biaxial-stress loading solver, 
         ! calculate the increment of VonMises strain 
         ! and the initial guess for macroscopic strain
         if (loadType==9 .or. &
             loadType==10) then
            
            ! VonMises strain increment:
            incrVonMises_imposed = DELTME*biaxial_StrainRate
            
            ! calculate initial guess for macroscopic strain
            call initialGuessForMacroscopicStrain(DeltaMacroStrain_initguess,ddsddeMacroscopic, &
                     incrVonMises_imposed,biaxial_StressRatioAngle, &
                     iStressMain,iStressLateral,deformMultiplierTensionCompression, &
                     specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv)
            DeltaMacroStrain = DeltaMacroStrain_initguess
            ! update the trial macroscopic strain with the init guess
            macroStrain_tau = macroStrain_t + DeltaMacroStrain
            
            if (IPROC==0 .and. DEBUG_GENERAL) write(*,*) 'initial guess - macroStrain_tau:',macroStrain_tau
            
            ! initialize the Lagrange multiplier associated w/ VM-Constraint
            lambda = 0.d0

         endif
         
         ! get the incremental boundary displacements (displ Boundary conditions)
         CALL DISINC(VTAU0,VBCS1,NBC1,DELTME,IDF,DIR_U,NSTEP,KSTEP, & 
                     G0XYZ,TIME1,NDF,MCRD1)
                     
         !**************************** calculate the 'settlement forces' to 
         !                                impose the incremental displ. BCs
         ! calculate the incremental displacement BCs:
         dU_BCs(1:NEQ) = VTAU0(1:NEQ) - VT(1:NEQ)
         
         if(DEBUG_RESID) &
            CALL printMatrix(RESHAPE(dU_BCs(1:NEQ), & 
            (/NX,3/),(/0.d0/),(/2,1/)),NX,3, & 
            'not improved - dU_BCs',542)
         
         ! extract the imposed values of constrained DOF/EQN for this proc --> uBC_proc.
         CALL extract_uBC_proc(dU_BCs,uBC_proc,EQ_uBC_proc, &
                               N_UPDATE_rxn,NEQ)
         
         ! calculate the 'settlement' forces due to imposed displacements
         CALL MATMUL_ColSparse(VAL_rxn,uBC_proc,F_uBC,   &
              IBINDX_rxn,ICOLPTR_rxn,N_UPDATE_rxn,NEQ,ISIZE_rxn)

         if(DEBUG_RESID) then
            write(542,*) 'uBC_proc',IPROC,uBC_proc(1:N_UPDATE_rxn)
            CALL printMatrix_Sparse(VAL_rxn,IBINDX,IROWPTR,  &
                                    N_UPDATE_rxn,NEQ,NNZ_LOC,   &
                                    1,'VAL_rxn',542)
         endif

         ! mpireduce F_uBC
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,F_uBC,NEQ, & 
                            MPI_DOUBLE_PRECISION, & 
                            MPI_SUM,MPI_COMM_WORLD,IERROR)

         if(DEBUG_RESID) &
         CALL printMatrix(RESHAPE(F(1:NEQ), & 
            (/NX,3/),(/0.d0/),(/2,1/)),NX,3, & 
            ' - F_resid',542)

         if(DEBUG_RESID) &
         CALL printMatrix(RESHAPE(F_uBC(1:NEQ), & 
         (/NX,3/),(/0.d0/),(/2,1/)),NX,3, & 
         ' - F_uBC',542)
         
         ! subtract from RHS (see formulation)
         F(1:NEQ) = F(1:NEQ) - F_uBC(1:NEQ)
          
         ! solve: DVTAU0 <<< SUPER LU <<< RESID (/w uBC)
         !*************************************************** SUPERLU << F
         ! copy RHS
         RESID_LOC(1:N_UPDATE)=F(IN_ST+1:IN_ST+N_UPDATE)
         ! save the [residual + constraint forces] for BFGS update
         F_prev(1:NEQ)=F(1:NEQ)
         ! check local SLU RHS vector length
         !CALL GET_COMPROWLOC_MATRIX(A, NROW_LOC=N_UPDATE_LOC)
         !CALL assertInt(N_UPDATE,N_UPDATE_LOC,'N_UPDATE_LOC')
         ! refresh sparse structures, ONLY IF re-factorization is needed
         if (.NOT.solver_LUisUpToDate) &
                  CALL convertSPARSE_f2C(IROWPTR,IBINDX,IROWPTR_LOC, &
                                         IBINDX_LOC,N_UPDATE,NNZ_LOC)
         ! set options
         CALL setSLUoptions(SLUoptions,solver_CalledBefore,solver_LUisUpToDate)
         ! solve
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         if(DEBUG_timer) CALL task_Start_superLU(solver_LUisUpToDate,MPI_WTIME())
         CALL F_PDGSSVX(SLUoptions, A, ScalePermstruct, RESID_LOC, & 
         N_UPDATE, NRHS, grid, LUstruct, SOLVEstruct, BERR, STAT, INFO)
         if(DEBUG_timer) CALL task_Record_superLU(solver_LUisUpToDate,MPI_WTIME(),lastRefactTime)
         ! update solver flags
         CALL justCalledSLU(solver_CalledBefore,solver_LUisUpToDate)
         CALL printSLUerror(INFO,BERR,NRHS,IPROC)
         if (INFO.NE.0) goto 701    ! terminate program if SLU error
         !****************************************** SUPERLU >> RESID_LOC
         ! initialize and save the solution
         DVTAU0(1:NEQ)=0.D0
         DVTAU0(IN_ST+1:IN_ST+N_UPDATE) = RESID_LOC(1:N_UPDATE)
         ! reduce
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,DVTAU0,NEQ, & 
         MPI_DOUBLE_PRECISION, & 
         MPI_SUM,MPI_COMM_WORLD,IERROR)
         
                                                   ! DVTAU0 = non-homogeneous BCs with improved initial guess
         VTAU0(1:NEQ) = VT(1:NEQ) + DVTAU0(1:NEQ)  ! improved VTAU0 using initial tangent modulus
         
         !==============================================================================
         ! END :             improve initial guess
         !==============================================================================
         ! we are now on a displacement field that satisfies the non-homogeneous u_BC's.

         ! --------------------- INITIALIZATION FOR Q-NEWTON ITERATIONS --------------------- !
         ! Tangent Stiffness Matrix is calculated  (LOADVC)
         ! and first iteration beyond after the initial guess is made (SuperLU)
         ! then the program enters the q-Newton loop
         ! ---------------------------------------------------------------------------------- !

         ! FTOL: calculate this after imposing the disp. BC's (at VTAU0) and calling LOADVC to get F_ext
         ! FTOL=1.0D-5 ! was 1.0D-4, before. 
                       ! CSR tests with 1E-8 confirmed convergence at 384 elem, 40x40x40 um mesh.
                       ! WARNING: finer meshes might require lower absolute force balance tolerance.
         ! now using an adaptive convergence criteria, see checkConvergence()
         ITEREF=0
         N_BFGS=0
         NUMUPD=0
         STEP=1.D0

         VAR(1:N_UPDATE*MAXREF)=0.D0
         WAR(1:N_UPDATE*MAXREF)=0.D0
!         V(:,:)=0.D0 ! you don't really need to initialize. they will be overwritten.
!         W(:,:)=0.D0 ! you don't really need to initialize. they will be overwritten.

!         VTAU(1:NEQ)=VTAU0(1:NEQ)   ! VTAU0 is the initial trial displacement field. 
!                                    ! it satisfies the (incremental) displacement BCs.
!                                    ! from this point on VTAU will be updated through q-Newton iterations.
         incVTAU(1:NEQ) = 0.d0
         VTAU = VT + DVTAU0 + incVTAU
         
         if(DEBUG_RESID) &
         CALL printMatrix(RESHAPE(DVTAU0(1:NEQ), & 
         (/NX,3/),(/0.d0/),(/2,1/)),NX,3, & 
         'improved (?)- DVTAU0',542)
                                    
			! DVTAU: increment in trial displacement field. v_tau(i+1) = v_tau(i) + dvtau
         ! VTAU:	trial displ. field at time tau
         ! DV: trial displ. field increment: v_tau(i) - v_t
         !
         ! ------|----------------|-------------|-----@----o---o---o---o-----------|------
         !  u_prev(:,1)     u_prev(:,2)        VT   VTAU0             VTAU         VTAU* 
         !                                      |\---/             \---/           (soln)
         !                                      |DVTAU0            DVTAU
         !                                       \--------- DV --------/ 

         ! each proc calculates local Stiff Mat and force vec. for their elements
         ! and sends out the coefficients of SKE and F to the processors that are assigned to handle the corresp. Global EQN #.
         ! finally F, the FEM force balance residual, is MPI_REDUCE'd
         ! all such communication is done in loadvc.
         
         CALL setFLAGS(IFLAG,solver_LUisUpToDate,stage_newtLoop,NSTEP,&
                       solverMethod,solver_CalledBefore, &
                       solver_reformedK_thisInc, &
                       trigger_cutBack,trigger_reformK)

         if(IFLAG.EQ.0) CALL task_Start(task_calcResidual_K,MPI_WTIME())
         if(IFLAG.NE.0) CALL task_Start(task_calcResidual,  MPI_WTIME())

         if(DEBUG_RESID) then
            write(542,*) 'after DVTAU0 -- VTAU',VTAU(1:NEQ)
            write(542,*) 'after DVTAU0 -- VT',VT(1:NEQ)
         endif

         calcStiffnessMatrix = (IFLAG == 0)
         call calcFBarStabilizers(VT,VTAU, &    ! displacement fields at time and tau
            IJK,G0XYZ,                     &    ! mesh,parallelism etc.
            calcStiffnessMatrix,success)        
            
         ! deniz modified - added output argument: F_ext
         CALL LOADVC(F,F_ext, &
         stressMacroscopic,ddsddeMacroscopic, & 
         defGradMacroscopic,dpMacroscopic, &
         volumeRVE,volume0RVE, &
         VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
         PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE,IN_ST,IROWPTR, & 
         IBINDX,VAL,RESID,NNZ_LOC,VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,   &
         N_UPDATE_rxn,ISIZE_rxn,IFLAG,PNEWDT)
                  
         if(IFLAG.EQ.0) CALL task_Record(task_calcResidual_K,MPI_WTIME(),task_time,printTasksTo)
         if(IFLAG.NE.0) CALL task_Record(task_calcResidual  ,MPI_WTIME(),task_time,printTasksTo)
         
         if(DEBUG_RESID) then
            write(542,*) 'after DVTAU0 -- F',F(1:NEQ)
            write(542,*) 'after DVTAU0 -- F_ext',F_ext(1:NEQ)
         endif            

         if(trigger_umatError) then                      ! UMAT error. R is corrupted.
            trigger_cutBack = .TRUE.                     ! cut back time step and recalculate increment
            if(IFLAG.EQ.0) trigger_reformK = .TRUE.      ! K is corrupted too. Mark for reformation next time loadVC is called
            goto 700                                     ! go to the cut-back handler (right after the q-Newton loop..)
         endif                                           ! this trigger will be caught there.
         
         G = DOT_PRODUCT(F(1:NEQ),DVTAU0(1:NEQ))

         CALL checkConvergence(solver_Converged,convergenceNorm,  &
                               eps_F,tol_F,F,F0,F_ext,F_ext_max,  &
                               eps_u,tol_u,DVTAU0,DVTAU0,VT,NEQ,  &
                               .TRUE.,   &  ! false: do not check convergence on displacement
                               conv_F, conv_u, tol_BelowNoise)
                               
         if(DEBUG_RESID) then
            write(542,*) 'after DVTAU0 -- conv_F, conv_u:',conv_F,conv_u
         endif    
         
         intData(1) = 0

         IF(IPROC==0) then
            if(tol_BelowNoise) then
               write(*,*) 'RESID_max INIT',MAXVAL(abs(F(1:NEQ))), &
                          ' - EXTERNAL FORCES NEGLIGIBLE'
            else
               write(*,*) 'RESID_max INIT',MAXVAL(abs(F(1:NEQ)))
            endif
         endif

         IF(IPROC.EQ.0) WRITE(*,'(A,I0,A,F8.1,A,F8.1)') &
         'ITER  ',ITEREF,'   conv_F:',conv_F,' conv_u:',conv_u
      
         ! also check if macroscopic strain has converged with initial guess (will happen if elastic)
         if (loadType==9 .or. &
             loadType==10) then
         
            ! calculate the update for macroscopic strain iterator using constrained solver
            DeltaMacroStrain_temp = DeltaMacroStrain
            lambda_temp = lambda
            call updateStrainIterator(DeltaMacroStrain_temp,stressMacroscopic,ddsddeMacroscopic, &
                        stressResidual,stressReaction,errorVM,lambda_temp, &
                        incrVonMises_imposed,biaxial_StressRatioAngle, &
                        iStressMain,iStressLateral, &
                        specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv)
            
            ! check convergence for macroscopic strain
            solverStrain_Converged = &
               (maxval(abs(stressResidual)) < eps_F*abs(stressReaction) .and. &
                errorVM < eps_u*incrVonMises_Imposed)
            
            if (.not.solverStrain_Converged) then
            
               solver_Converged = .false.
               
            endif

         endif
         
         
         !******************************************************************************************!
         !                                         QSINT                                            !
         !******************************************************************************************!
         !
         !   Quasi-Newton iterations start here
         !   in the following iterations, SKE will not be recalculated/inverted
         !   only the residual is calculated, initial SKE_inv will be used/modified by a Broyden method
         !
         !  repeat until convergence, or a signalling of a cut-back
         DO WHILE((.NOT.solver_Converged).AND.(.NOT.trigger_cutBack))

            IF(ITEREF.GE.MAXREF)THEN   ! Broyden updates did not converge
               
               trigger_BroyDidNotConv = .TRUE.  ! broyden did not converge
               exit                             ! this will be handled outside the loop
               
            ENDIF
         
            ITEREF=ITEREF+1
            
           write(strNSTEP,'(I0)') NSTEP
           write(strNITER,'(I0)') ITEREF
           
            ! *************** BFGS(F_prev,RESID,DVTAU) ************
            ! *********************************************************
            ! pre-multiply the residual before the back-substitution
            if (solverMethod.EQ.solmeth_BFGS.OR.   &
                solverMethod.EQ.solmeth_BFGSkeep) then
               
               CALL task_Start(task_BFGS_operation,MPI_WTIME())
               F_BFGS(1:NEQ) = F(1:NEQ)
               CALL BFGS_pre_F(W,V,F_BFGS,NEQ,MAXREF,N_BFGS)
               RESID(1:N_UPDATE) = F_BFGS(IN_ST+1:IN_ST+N_UPDATE)
               CALL task_Record(task_BFGS_operation,MPI_WTIME(),task_time,-1)
            endif
            
            if(DEBUG_RESID) write(542,*) 'BEFORE SLU of ITERATION',ITEREF
            
            ! back-substitution
            ! solve: DVTAU <<< SUPER LU <<< RESID (uBC=0)
            !*************************************************** SUPERLU << RESID
            ! copy RHS
            RESID_LOC(1:N_UPDATE)=RESID(1:N_UPDATE)
            ! copy for BFGS update
            F_prev(1:NEQ)=F(1:NEQ)
            ! check local SLU RHS vector length
            !CALL GET_COMPROWLOC_MATRIX(A, NROW_LOC=N_UPDATE_LOC)
            !CALL assertInt(N_UPDATE,N_UPDATE_LOC,'N_UPDATE_LOC')
            ! refresh sparse structures, ONLY IF re-factorization is needed
            if (.NOT.solver_LUisUpToDate) &
                     CALL convertSPARSE_f2C(IROWPTR,IBINDX,IROWPTR_LOC, &
                                            IBINDX_LOC,N_UPDATE,NNZ_LOC)
            ! set options
            CALL setSLUoptions(SLUoptions,solver_CalledBefore,solver_LUisUpToDate)
            ! solve
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
            if(DEBUG_timer) CALL task_Start_superLU(solver_LUisUpToDate,MPI_WTIME())
            CALL F_PDGSSVX(SLUoptions, A, ScalePermstruct, RESID_LOC, & 
            N_UPDATE, NRHS, grid, LUstruct, SOLVEstruct, BERR, STAT, INFO)
            if(DEBUG_timer) CALL task_Record_superLU(solver_LUisUpToDate,MPI_WTIME(),lastRefactTime)
            ! update solver flags
            CALL justCalledSLU(solver_CalledBefore,solver_LUisUpToDate)
            CALL printSLUerror(INFO,BERR,NRHS,IPROC)
            if (INFO.NE.0) goto 701    ! terminate program if SLU error
            !****************************************** SUPERLU >> RESID_LOC

            ! solution was saved into RESID_LOC. copy to DVTAU (global-size).         
            DVTAU(1:NEQ)=0.D0
            DVTAU(IN_ST+1:IN_ST+N_UPDATE)=RESID_LOC(1:N_UPDATE)

            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,DVTAU,NEQ, & 
            MPI_DOUBLE_PRECISION, & 
            MPI_SUM,MPI_COMM_WORLD,IERROR)
            
            if (solverMethod.EQ.solmeth_BFGS.OR.   &
                solverMethod.EQ.solmeth_BFGSkeep) then

               CALL task_Start(task_BFGS_operation,MPI_WTIME())
               CALL BFGS_post_U(W,V,DVTAU,NEQ,MAXREF,N_BFGS)
               CALL task_Record(task_BFGS_operation,MPI_WTIME(),task_time,-1)

            elseif (solverMethod.EQ.solmeth_Broyden.AND.ITEREF.GT.1) then

               CALL task_Start(task_Broy_update,MPI_WTIME())
               CALL BROY_PARALLEL(DVTAU,VAR,WAR,NUMUPD,NMPD1,MAXREF, &
                  IN_ST,N_UPDATE,NEQ,STEP)
               CALL task_Record(task_Broy_update,MPI_WTIME(),task_time,-1)
            endif
            
            if(IPROC.EQ.0.and.DEBUG_RESID) then
               CALL printMatrix(RESHAPE(DVTAU(1:NEQ), & 
               (/NX,3/),(/0.d0/),(/2,1/)),NX,3, & 
               'Broyden-updated DVTAU',542)
            endif

            STEP=1.D0

            G0 = DOT_PRODUCT(F(1:NEQ),DVTAU(1:NEQ))
            
            incVTAU_init(1:NEQ) = incVTAU(1:NEQ)
            incVTAU(1:NEQ)=incVTAU_init(1:NEQ)+STEP*DVTAU(1:NEQ)	! DVTAU is the displ. increment, i.e. soln of the linearized balance laws
                                                            ! it satisfies homogeneous displ. BC's at all Dirichlet nodes.
            VTAU(1:NEQ) = VT + DVTAU0 + incVTAU
            ! incVTAU_init: trial displacement field at the beginning of line search
            ! DVTAU: direction of iteration step of displacement . v_tau(i+1) = v_tau(i) + step*dvtau
            ! step: iteration step length. start with step=1.0 as this qill give quadratic conv. towards the end
            ! VTAU:	trial displ. field for time tau
            ! DV: trial displ. field increment: v_tau(i) - v_t
            !
            ! ------|----------------|-------------|-----@----o---o----------o-----------|------
            !  u_prev(:,1)     u_prev(:,2)        VT            incVTAU_init    VTAU      VTAU* 
            !                                      |              \----------/         (soln)
            !                                      |                 S*DVTAU
            !                                       \------- DV -------------/ 

            CALL setFLAGS(IFLAG,solver_LUisUpToDate,stage_newtLoop,NSTEP,&
                          solverMethod,solver_CalledBefore, &
                          solver_reformedK_thisInc, &
                          trigger_cutBack,trigger_reformK)
                          
            if(IFLAG.EQ.0) CALL task_Start(task_calcResidual_K,MPI_WTIME())
            if(IFLAG.NE.0) CALL task_Start(task_calcResidual,  MPI_WTIME())

            calcStiffnessMatrix = (IFLAG == 0)
            call calcFBarStabilizers(VT,VTAU, &    ! displacement fields at time and tau
               IJK,G0XYZ,                     &    ! mesh,parallelism etc.
               calcStiffnessMatrix,success)

            ! first, evaluate at STEP=1
            CALL LOADVC(F,F_ext, &
            stressMacroscopic,ddsddeMacroscopic, & 
            defGradMacroscopic,dpMacroscopic, &
            volumeRVE,volume0RVE, &
            VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
            PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE,IN_ST,IROWPTR, & 
            IBINDX,VAL,RESID,NNZ_LOC,VAL_rxn,ICOLPTR_rxn,IBINDX_rxn, & 
            N_UPDATE_rxn,ISIZE_rxn,IFLAG,PNEWDT)
            
            if(IFLAG.EQ.0) CALL task_Record(task_calcResidual_K,MPI_WTIME(),task_time,printTasksTo)
            if(IFLAG.NE.0) CALL task_Record(task_calcResidual  ,MPI_WTIME(),task_time,printTasksTo)

            if(trigger_umatError) then                      ! UMAT error. R is corrupted.
               trigger_cutBack = .TRUE.                     ! cut back time step and recalculate increment
               if(IFLAG.EQ.0) trigger_reformK = .TRUE.      ! K is corrupted too. Mark for reformation next time loadVC is called
               exit                                         ! exit the q-Newton loop. trigger will be caught outside.
            endif
            
            G = DOT_PRODUCT(F(1:NEQ),DVTAU(1:NEQ))
            
            ! ***************************************! >>>>>> line search >>>>>> !
            ! ***************************************! >>>>>>>>>>>>>>>>>>>>>>>>> !
            ! do line search, if Wolfe-2 condition does not hold (i.e., G is too steep)
            if (abs(G).GT.WLF_c2*abs(G0).AND. & 
                lschMethod.GT.0) then

               ! lschMethod = 1 : use Deniz's method
               ! lschMethod = 2 : use the method described in Gilbert Strang's 1979 paper
               if(lschMethod.EQ.2) then
              CALL lineSearchSTRANG(STEP,G0,G,F,F_ext, &
                      stressMacroscopic,ddsddeMacroscopic, & 
                      defGradMacroscopic,dpMacroscopic, &
                      volumeRVE,volume0RVE, &
                      VT,VTAU,DVTAU,incVTAU_init,G0XYZ, & 
                      NEQ,WLF_c2,STEP_max,G_history,G_n,lschMethod,info, & 
                      XUND,IJK,IMID,NPROPSG,PROMAT2D,DELTME, & 
                      NBC1,IBC1,KINC,NDF,N_UPDATE,IN_ST,IROWPTR, & 
                      IBINDX,VAL,RESID,NNZ_LOC, & 
                VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                      NGAUSS,PNEWDT,IPROC)
               else
              CALL lineSearch(STEP,G0,G,F,F_ext, &
                      stressMacroscopic,ddsddeMacroscopic, & 
                      defGradMacroscopic,dpMacroscopic, &
                      volumeRVE,volume0RVE, &
                      VT,VTAU,DVTAU,incVTAU_init,G0XYZ, & 
                      NEQ,WLF_c2,STEP_max,G_history,G_n,lschMethod,info, & 
                      XUND,IJK,IMID,NPROPSG,PROMAT2D,DELTME, & 
                      NBC1,IBC1,KINC,NDF,N_UPDATE,IN_ST,IROWPTR, & 
                      IBINDX,VAL,RESID,NNZ_LOC, & 
                 VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                      NGAUSS,PNEWDT,IPROC)
               endif
            
               if(trigger_umatError) then                      ! UMAT error during line search. R is corrupted.
                  trigger_cutBack = .TRUE.                     ! cut back time step and recalculate increment
                  if(IFLAG.EQ.0) trigger_reformK = .TRUE.      ! K is corrupted too. Mark for reformation next time loadVC is called
                  exit                                         ! exit the q-Newton loop. trigger will be caught outside.
               endif
               
               if(IPROC.EQ.0.AND.DEBUG_linesearch) then
                 write(995,*) 'error status: ',info
                 CALL printMatrix(G_history(:,:),2,G_n, & 
                 'G_history @T_STEP '//strNSTEP//', ITER '//strNITER, & 
                                                                 995)
               endif
            endif
            ! >>>>>> line search >>>>>>!***************************************!
            ! >>>>>>>>>> done >>>>>>>>>!***************************************!
            ! VTAU is the current trial displ. field.
            
            if (solverMethod.EQ.solmeth_Broyden) then ! save the displ. increment for the upcoming Broyden updates
            
               NMPD1=NUMUPD*N_UPDATE

               WAR(NMPD1+1:NMPD1+N_UPDATE) = &
                                 STEP*DVTAU(IN_ST+1:IN_ST+N_UPDATE)
            endif

            if (solverMethod.EQ.solmeth_BFGS.OR.   &
                solverMethod.EQ.solmeth_BFGSkeep) then ! save w_i for the upcoming Broyden updates
                
               CALL task_Start(task_BFGS_operation,MPI_WTIME())
               CALL BFGS_addWV(W,V,DVTAU,F,F_prev,G,G0,STEP,NEQ,MAXREF,N_BFGS)
               CALL task_Record(task_BFGS_operation,MPI_WTIME(),task_time,-1)
               
            endif

            CALL checkConvergence(solver_Converged,convergenceNorm,  &
                                  eps_F,tol_F,F,F0,F_ext,F_ext_max,  &
                                  eps_u,tol_u,DVTAU,DVTAU0,VT,NEQ,         &
                                  .TRUE.,   &  ! check convergence on displacement too
                                  conv_F, conv_u, tol_BelowNoise)
                                  
            ! if strain-controlled biaxial stress Loading:
            ! check if macroscopic strain has converged and 
            ! if not, apply an update and continue QN-iterations
            !==================================== START
            strStrainIterMessage = ''
            if (solver_Converged .and. &
               (loadType==9 .or. loadType==10)) then
               
               ! calculate the update for macroscopic strain iterator using constrained solver
               DeltaMacroStrain_prev = DeltaMacroStrain
               call updateStrainIterator(DeltaMacroStrain,stressMacroscopic,ddsddeMacroscopic, &
                           stressResidual,stressReaction,errorVM,lambda, &
                           incrVonMises_imposed,biaxial_StressRatioAngle, &
                           iStressMain,iStressLateral, &
                           specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv)
               DeltaDeltaMacroStrain = DeltaMacroStrain - DeltaMacroStrain_prev
               
               ! check convergence for macroscopic strain
               solverStrain_Converged = &
                  (maxval(abs(stressResidual)) < eps_F*abs(stressReaction) .and. &
                   errorVM < eps_u*incrVonMises_Imposed)
               
               if (.not.solverStrain_Converged) then
                              
                  ! update the macro strain at time tau
                  macroStrain_tau = macroStrain_t + DeltaMacroStrain
                  
                  if (useTangentStiffnessStrainUpdates) then
                     
                     ! get the incremental boundary displacements (displ Boundary conditions)
                     CALL DISINC(dU_BCs,VBCS1,NBC1,DELTME,IDF,DIR_U,NSTEP,KSTEP, & 
                                 G0XYZ,TIME1,NDF,MCRD1)
                                 
                     !**************************** calculate the 'settlement forces' to 
                     !                                impose the incremental displ. BCs
                     ! calculate the incremental displacement BCs due to strain update
                     dU_BCs(1:NEQ) = dU_BCs - VTAU
                     
                     ! extract the imposed values of constrained DOF/EQN for this proc --> uBC_proc.
                     CALL extract_uBC_proc(dU_BCs,uBC_proc,EQ_uBC_proc, &
                                           N_UPDATE_rxn,NEQ)
                     
                     ! calculate the 'settlement' forces due to imposed displacements
                     CALL MATMUL_ColSparse(VAL_rxn,uBC_proc,F_uBC,   &
                          IBINDX_rxn,ICOLPTR_rxn,N_UPDATE_rxn,NEQ,ISIZE_rxn)

                     ! mpireduce F_uBC
                     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
                     CALL MPI_ALLREDUCE(MPI_IN_PLACE,F_uBC,NEQ, & 
                                        MPI_DOUBLE_PRECISION, & 
                                        MPI_SUM,MPI_COMM_WORLD,IERROR)

                     ! convention: external forces are negative on RHS (see formulation)
                     F(1:NEQ) = - F_uBC(1:NEQ)
                      
                     ! solve: DVTAU0 <<< SUPER LU <<< RESID (/w uBC)
                     !*************************************************** SUPERLU << F
                     ! copy RHS
                     RESID_LOC(1:N_UPDATE)=F(IN_ST+1:IN_ST+N_UPDATE)
                     ! save the [residual + constraint forces] for BFGS update
                     ! F_prev(1:NEQ)=F(1:NEQ)
                     ! check local SLU RHS vector length
                     !CALL GET_COMPROWLOC_MATRIX(A, NROW_LOC=N_UPDATE_LOC)
                     !CALL assertInt(N_UPDATE,N_UPDATE_LOC,'N_UPDATE_LOC')
                     ! refresh sparse structures, ONLY IF re-factorization is needed
                     if (.NOT.solver_LUisUpToDate) &
                              CALL convertSPARSE_f2C(IROWPTR,IBINDX,IROWPTR_LOC, &
                                                     IBINDX_LOC,N_UPDATE,NNZ_LOC)
                     ! set options
                     CALL setSLUoptions(SLUoptions,solver_CalledBefore,solver_LUisUpToDate)
                     ! solve
                     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
                     if(DEBUG_timer) CALL task_Start_superLU(solver_LUisUpToDate,MPI_WTIME())
                     CALL F_PDGSSVX(SLUoptions, A, ScalePermstruct, RESID_LOC, & 
                     N_UPDATE, NRHS, grid, LUstruct, SOLVEstruct, BERR, STAT, INFO)
                     if(DEBUG_timer) CALL task_Record_superLU(solver_LUisUpToDate,MPI_WTIME(),lastRefactTime)
                     ! update solver flags
                     CALL justCalledSLU(solver_CalledBefore,solver_LUisUpToDate)
                     CALL printSLUerror(INFO,BERR,NRHS,IPROC)
                     if (INFO.NE.0) goto 701    ! terminate program if SLU error
                     !****************************************** SUPERLU >> RESID_LOC
                     ! initialize and save the solution
                     dU_BCs(1:NEQ)=0.D0
                     dU_BCs(IN_ST+1:IN_ST+N_UPDATE) = RESID_LOC(1:N_UPDATE)
                     ! reduce
                     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
                     CALL MPI_ALLREDUCE(MPI_IN_PLACE,dU_BCs,NEQ, & 
                     MPI_DOUBLE_PRECISION, & 
                     MPI_SUM,MPI_COMM_WORLD,IERROR)
                     
                     incVTAU(1:NEQ) = incVTAU(1:NEQ) + dU_BCs(1:NEQ)  ! dU_BCs = update of displacement field due to strain update
                                 
                     !==============================================================================
                  else

                     ! update the displacement field due to updated strain
                     ! this cheaper method applies an
                     ! incremental affine transformation (assoc. with strain update) 
                     ! to all the nodes (boundary + bulk)
                     call transformAffineDisplacementField( &
                                 incVTAU, &
                                 macroStrain_tau - DeltaDeltaMacroStrain, &
                                 macroStrain_tau, &
                                 G0XYZ,domainRange0,NX)
                        
                  endif        
                     
                  ! update the displacement field at time tau
                  VTAU(1:NEQ) = VT + DVTAU0 + incVTAU
                        
                  ! update the residual for the new displacement field with the applied strain increment (necessary for healthy BFGS updates)
                  !===============================================
                  CALL setFLAGS(IFLAG,solver_LUisUpToDate,stage_newtLoop,NSTEP,&
                                solverMethod,solver_CalledBefore, &
                                solver_reformedK_thisInc, &
                                trigger_cutBack,trigger_reformK)
                                
                  if(IFLAG.EQ.0) CALL task_Start(task_calcResidual_K,MPI_WTIME())
                  if(IFLAG.NE.0) CALL task_Start(task_calcResidual,  MPI_WTIME())

                  calcStiffnessMatrix = (IFLAG == 0)
                  call calcFBarStabilizers(VT,VTAU, &    ! displacement fields at time and tau
                     IJK,G0XYZ,                     &    ! mesh,parallelism etc.
                     calcStiffnessMatrix,success)

                  ! first, evaluate at STEP=1
                  CALL LOADVC(F,F_ext, &
                  stressMacroscopic,ddsddeMacroscopic, & 
                  defGradMacroscopic,dpMacroscopic, &
                  volumeRVE,volume0RVE, &
                  VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                  PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE,IN_ST,IROWPTR, & 
                  IBINDX,VAL,RESID,NNZ_LOC,VAL_rxn,ICOLPTR_rxn,IBINDX_rxn, & 
                  N_UPDATE_rxn,ISIZE_rxn,IFLAG,PNEWDT)
                  
                  if(IFLAG.EQ.0) CALL task_Record(task_calcResidual_K,MPI_WTIME(),task_time,printTasksTo)
                  if(IFLAG.NE.0) CALL task_Record(task_calcResidual  ,MPI_WTIME(),task_time,printTasksTo)

                  if(trigger_umatError) then                      ! UMAT error. R is corrupted.
                     trigger_cutBack = .TRUE.                     ! cut back time step and recalculate increment
                     if(IFLAG.EQ.0) trigger_reformK = .TRUE.      ! K is corrupted too. Mark for reformation next time loadVC is called
                     exit                                         ! exit the q-Newton loop. trigger will be caught outside.
                  endif
                  !===============================================
               
                  ! mark FEM-solver as not converged
                  solver_Converged = .false.
                  
                  ! reset BFGS records:
                  N_BFGS=0
                  NUMUPD=0
                  VAR(:)=0.D0
                  WAR(:)=0.D0
                     
                  strStrainIterMessage = ' Macro.Strain Updated'
               else
               
                  strStrainIterMessage = ' Macro.Strain Converged'
                  
               endif

            endif

            if(tol_BelowNoise) then
               strExtForces = ' ext.forces~0'
            else
               strExtForces = ''
            endif
            
            if(DEBUG_RESID) then
                  write(542,'(A,I0,A,F8.1,A,F8.1,A,A)') 'ITER  ',ITEREF,'  &
                        conv_F:',conv_F,' conv_u:',conv_u, &
                        trim(strExtForces),trim(strStrainIterMessage)
            endif    
               
            IF(IPROC==0) then
                  write(*,'(A,I0,A,F8.1,A,F8.1,A,A)') 'ITER  ',ITEREF,'  &
                        conv_F:',conv_F,' conv_u:',conv_u, &
                        trim(strExtForces),trim(strStrainIterMessage)
            endif
         ENDDO
         !******************************************************************************************!
         !                                    QSINT END                                             !
         !******************************************************************************************!
         ! END of the QUASI-NEWTON LOOP.
         ! condition satisfied: XMAXV.LE.FTOL .OR. PNEWDT.LT.1.D0
         ! either we have converged on a displacement field: VTAU,
         ! or iterations did not converge, or a umat error occured
         
700      continue ! the cut back handler

         if(trigger_BroyDidNotConv) then  ! the Broyden updates have not converged

            ! if the method is broyden-keep or modified newton-keep, 
            ! mark for the recalculation of the stiffness matrix 
            ! and send back to the beginning of the increment

             ! reset the trigger. it has done its job.
            trigger_BroyDidNotConv = .FALSE.
                           
            ! decide whether to reform K_inv, or to cut back time step
            ! ------------------------------------
            ! 1) if not already reformed in this step,
            if (.NOT.solver_reformedK_thisInc) then
            
               ! report CPU time/convergence performance
               if(IPROC.EQ.0.AND.DEBUG_timer) write(994,*) STATME,'-->',FINTME,'Broy ReformK.'
               IF(IPROC.EQ.0)THEN
                  WRITE(*,*)'THE BROY. ITERATION DID NOT CONVERGE. WILL RECALCULATE STIFFNESS MATRIX.'
               ENDIF
               
               trigger_reformK = .TRUE.   ! trigger reformation for the next loadVC call
               FINTME=FINTME-DELTME       ! and retry the same cycle
               PNEWDT = 1.0d0
               
               CYCLE                      
            else
            ! 2) if already reformed in this step
            
               ! report CPU time/convergence performance
               if(IPROC.EQ.0.AND.DEBUG_timer) write(994,*) STATME,'-->',FINTME,'Broy cutback'
               IF(IPROC.EQ.0)THEN
                  WRITE(*,*)'THE BROY. ITERATION DID NOT CONVERGE EVEN WITH UPDATED STIFFNESS MATRIX. WILL CUTBACK TIMESTEP'
               ENDIF
               
               trigger_cutBack = .TRUE.
               PNEWDT = 0.5d0

            endif
            
         endif 

         if (trigger_cutBack) then
         
            ICUT=ICUT+1
            
            if(IPROC.EQ.0)then
               ! report CPU time/convergence performance
               if(DEBUG_timer.and.trigger_umatError) &
                  write(994,*) STATME,'-->',FINTME,'UMAT error.'
                  
               if(trigger_umatError) then ! list local umat errors coming from all processors
                  write(*,*) '! UMAT ERRORS OCCURRED: ', listLocalUmatErrorsOccurred(:)
               endif
               WRITE(*,*)'CUTBACK',ICUT
            ENDIF
            
            IF(ICUT.GT.MAXCUT)THEN
               IF(IPROC.EQ.0)THEN
                  WRITE(131,*)'TOO MANY CUTBACKS. TERMINATING...'
                  WRITE(*,*)'TOO MANY CUTBACKS. TERMINATING...'
               ENDIF
               goto 701 ! deallocate memory and terminate program
            ENDIF

            FINTME=FINTME-DELTME
            DELTME=PNEWDT*DELTME
            prevTimeCutBackForSpecial = .false. ! reset the adaptive-time step booster for post-special-time
            EXPORT_this_step = .false.

            IF(DELTME.LT.XMINDT)THEN
               IF(IPROC.EQ.0)THEN
                  WRITE(131,*)'TOO SMALL INCREMENT REQUIRED',DELTME
                  WRITE(*,*)'TOO SMALL INCREMENT REQUIRED',DELTME
               ENDIF
               goto 701 ! deallocate memory and terminate program
            ENDIF

            IF(IPROC.EQ.0)THEN
               WRITE(131,*)'USER ROUTINE REQUESTS SMALLER INCREMENT'
            ENDIF
            
            CYCLE  ! retry increment
         else
            ICUT = 0 ! all is good. reset the counter for successive cutbacks
         endif
         
         if(.NOT.solver_Converged) then
            write(*,*) 'unhandled case'
            goto 701 ! terminate program for debugging
         endif

         !*********************************************************************!
         !**************** converged to an incremental solution ***************!
         !*********************************************************************!
         ! convergence failures/umat errors, if any, were caught above 
         ! to re-start the increment. 
         ! at this point, we are guaranteed to have converged.
         !
         ! 1) update variables
         ! 2) do postProcessing
         ! 3) exportData
         !*********************************************************************!
         !             (1)     UPDATE VARIABLES                                !
         !*********************************************************************!
         
         NITER=ITEREF+1

         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'STEP=',NSTEP,'TIME=',FINTME
            WRITE(131,*)'FRAC COMPLETED=',FINTME/TOTME,  &
                        'TIME INCR=',DELTME
            WRITE(131,*)'NITER=',NITER
         ENDIF
         
         
         
         ! report CPU time/convergence performance
         if(IPROC.EQ.0.AND.DEBUG_timer) write(994,*) STATME,'-->',FINTME,'success. #iterations:',NITER+1

         ! 	VTAU:	trial displ. field at time tau
         !	DV: initial trial displ. field increment: v_tau(i) - v_t
         !
         ! ------|----------------|-------------|-----o----o---o---o---o-----------|------
         !  u_prev(:,1)     u_prev(:,2)        VT                     VTAU         VTAU* 
         !                                      |                  \---/           (soln)
         !                                      |                  DVTAU
         !                                       \------- DV ------/ 

         if (DEBUG_RESID.AND.IPROC.EQ.0) then
         write(970,*) 'just converged',NITER
         do iNode=1,NX
            write(970,*) VTAU(3*(iNode-1)+1:3*(iNode-1)+3) - VT(3*(iNode-1)+1:3*(iNode-1)+3)
         enddo
         endif
         
         ! ****** update displacement field solution ****** !
         ! update nodal positions,
         !     v_t <-- v_tau
         GPXYZ(1:NEQ)       = G0XYZ(1:NEQ) + VTAU(1:NEQ)
         u_prev(1:NEQ,1)    = u_prev(1:NEQ,2)
         u_prev(1:NEQ,2)    = VT(1:NEQ)
         VT(1:NEQ)          = VTAU(1:NEQ)
         ! update macroscopic strain
         macroStrain_t = macroStrain_tau
         ! update variables for improving VTAU0
         DELTME_PREV(1)     = DELTME_PREV(2)
         DELTME_PREV(2)     = DELTME
         time_prev(1)       = time_prev(2)
         time_prev(2)       = STATME
         lastConvergedTime  = FINTME

         ! ******* set the new time step ********* !
         ! takes into account:
         !  - pnewdt=2.0 suggestion from umat
         !  - # of steps to convergence (adaptive time-stepping)
         !  - maximum time step, XMXDT
         !  - end of the simulation, TOTME
         !  - load-case specific critical points:
         !     * ramping time in creep loading
         !     * max. amplitude in cyclic displacement/force loading
         !     * all time points in the arbitrary displacement/strain/force-time history         
         ! - also returns if the next time step is marked for data-export
         CALL setDELTME(DELTME,PNEWDT,FINTME,XMXDT,XMINDT,TOTME,  &
                        adaptive_time_stepping,suppress_specialPoints,&
                        NITER,NITER_REF, & 
                        prevTimeCutBackForSpecial,DELTME_lastBeforeSpecTime, &
                        EXPORT_next_step)
      
         ! *********** update state variables *********** !
         ! svars1: state variables from time t (converged) 
         ! svars2: trial state variables at time tau (not converged)
         ! svars1 <- svars2
         DO II=1,NSVARS1*NELNOX*NGAUSS
            SVARS1(II)=SVARS2(II)         ! SVARS2 was updated at the end of ELST1, through common block
         ENDDO
         ! NOTE: SVARS1 and SVARS2 store the state variables of the elements handled by individual processes
         !******* Collect the State Variables of all Elements ********** !
         !        in the master processor for post-proc                  !
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         if (IPROC==0) then   ! receive from individual processors sequentially
         
            SVARS_all(:)=0.D0
            
            ! first, copy your local data into the full array
            SVARS_all((NELST-1)*NGAUSS*NSTATE+1:NELEND*NGAUSS*NSTATE) = SVARS1(1:NELNOX*NGAUSS*NSTATE)

            do tProc = 1,NPROCS-1
            
               tNELST = NELST_proc(tProc)
               tNELEND = NELEND_proc(tProc)
               tNELNOX = NELNOX_proc(tProc)
               tSIZE = tNELNOX*NGAUSS*NSTATE
               
               CALL MPI_RECV(SVARS_all((tNELST-1)*NGAUSS*NSTATE+1:tNELEND*NGAUSS*NSTATE), &
                             tSIZE,MPI_DOUBLE_PRECISION, &
                             tProc,969,MPI_COMM_WORLD,ISTATUS,IERROR)

               if(IERROR.NE.MPI_SUCCESS)  exit
            enddo
         else                 ! send to master processor

               tSIZE = NELNOX*NGAUSS*NSTATE
               CALL MPI_SSEND(SVARS1(1:NELNOX*NGAUSS*NSTATE), &
                             tSIZE,MPI_DOUBLE_PRECISION, &
                             0,969,MPI_COMM_WORLD,IERROR)
         
         endif
         if(IERROR.NE.0) then
            write(*,*) 'MPI Error while collecting state variables'
            write(131,*) 'MPI Error while collecting state variables'
            close(131)               
            CALL MPI_FINALIZE(IERROR)
            stop
         endif
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         
         
         ! SUCCESS _ DEBUG _ REVERT THIS !! ! 
         if(DEBUG_UMAT.and.IPROC==2) then
            write(333,*) FINTME,realData((intData(5000)-1)*6+2001:(intData(5000)-1)*6+2006)
         endif
         
         !********************************************************************!
         !               (2)    POST-PROCESSING                               !
         !********************************************************************!
         ! determine if needed
         CALL doAtPeriod(FINTME,NSTEP,postProcT,postProc_lastT, & 
                         postProcN,do_postProc)
         if (postProcN == 0) do_postProc = .false. ! bypass iF explicitly disabled
                         
         if((dabs(FINTME-TOTME).LT.EXPORT_allAtLast).OR. &
            EXPORT_this_step) then
            ! this time step was marked for data export
            ! do post process, and then export all.
            ! do_postProc = .TRUE. ! -- WHY?
         endif
         
         if (IPROC.eq.0.AND.do_postProc) &
            CALL task_Start(task_PostProcess,MPI_WTIME())

         ! ******* NON-LOCAL DAMAGE / REGULARISATION *********
         ! ***************************************************
         if (damage_nonlocal) then
         
            if(IPROC==0) then
               do iElem=1,NELX
                  ! reset nonlocalDamage
                  !---- not used for this material
                  !SVARS_all((iElem-1)*NGAUSS*NSTATE+327:(iElem-1)*NGAUSS*NSTATE+330) = 0.d0
               enddo

               do iElem=1,NELX
               
                  nElem = elemNonLocalCount(iElem)
                  
                  IS = (iElem-1)*NGAUSS*NSTATE

                  ! read from local damage variable
                  localDamage(1:4) = 0.d0                         !---not used for this material
                  !localDamage(1:4) = SVARS_all(IS+378:IS+381)    !---not used for this material

                  do elemIdx=1,nElem
                     jElem = elemNonLocalList(elemIdx,iElem)
                     JS = (jElem-1)*NGAUSS*NSTATE
                     weight = elemNonLocalWeights(elemIdx,iElem)

                     !---not used for this material
                     !SVARS_all(JS+327:JS+330) = &
                     !   SVARS_all(JS+327:JS+330) + weight*localDamage(1:4)

                  enddo

               enddo
            endif
               
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)   ! non-local damage calculating done.
                                                      ! now, distribute to the other processors
            
            if (IPROC==0) then                        ! send non-local damage
            
               ! save non-local damage to your local state variables 327-330
               ! for the next iteration step
               do iElem=NELST,NELEND
                  IS = (iElem-1)*NGAUSS*NSTATE
                  iEl_proc = iElem - NELST + 1
                  !---not used for this material
                  !SVARS1((iEl_proc-1)*NGAUSS*NSTATE+327: & 
                  !       (iEl_proc-1)*NGAUSS*NSTATE+330) & 
                  !            = SVARS_all(IS+327:IS+330)
               enddo
                 
     
               ! send to the others:
               do tProc = 1,NPROCS-1
               
                  tNELST = NELST_proc(tProc)
                  tNELEND = NELEND_proc(tProc)
                  tNELNOX = NELNOX_proc(tProc)
                  tSIZE = tNELNOX*NGAUSS*4
                  
                  ! prepare package
                  do iElem=tNELST,tNELEND
                     IS = (iElem-1)*NGAUSS*NSTATE
                     iEl_proc = iElem - tNELST + 1
                     !---not used for this material
                     !procNonLocalDamage((iEl_proc-1)*4+1: & 
                     !           (iEl_proc-1)*4+4) & 
                     !            = SVARS_all(IS+327:IS+330)
                  enddo
                  
                  ! send package
                  CALL MPI_SSEND(procNonLocalDamage(1:tSIZE), &
                                tSIZE,MPI_DOUBLE_PRECISION, &
                                tProc,969,MPI_COMM_WORLD,IERROR)
                                
                  if(IERROR.NE.0) exit
               enddo

            else                                 

               ! receive non-local damage from the master processor
               tSIZE = NELNOX*NGAUSS*4
               CALL MPI_RECV(procNonLocalDamage(1:tSIZE), &
                             tSIZE,MPI_DOUBLE_PRECISION, &
                             0,969,MPI_COMM_WORLD,ISTATUS,IERROR)
                             
               ! save non-local damage to your local state variables 327-330
               ! for the next iteration step
               do iElem=NELST,NELEND
                  IS = (iElem-1)*NGAUSS*NSTATE
                  iEl_proc = iElem - NELST + 1
                  !---- not used for this material
                  !SVARS1((iEl_proc-1)*NGAUSS*NSTATE+327: & 
                  !       (iEl_proc-1)*NGAUSS*NSTATE+330) & 
                  !            = procNonLocalDamage((iEl_proc-1)*4+1: &
                  !                          iEl_proc*4)
               enddo

            endif
      
         endif

         ! 1) error analysis
         ! 2) calculate Nye tensor, save it back to state variables
         !
         ! ****(1) FEM ERROR ANALYSIS ************
         ! ***************************************
         
! ---------------- CALCULATE THE NYE TENSOR ----------------- !
! ---------------- ON THE MASTER PROCESSOR ------------------ !
! ----------------- AND SEND TO THE OTHERS ------------------ !

! -------------- RECOVER Fp at NODES --------------- !
         
         if(do_postProc) then
         
            if(IPROC==0) write(*,*) 'Post-processing ... '// & 
                          'calculating Nye tensor field'

            ! extract gauss values of Fp            
            if(IPROC==0) &            
               CALL getGaussValues(gaussValues,SVARS_all, & 
                                   67,nDOF_FP, & ! Fp: 67 to 75 ! WARNING ! Fp - is stored ROW-MAJOR in STATE VARIABLES
                                   NSTATE,NELX,NGAUSS)

            error = 0
            if(recoveryMethod==rec_AVG)then ! nodal averaging

               if(IPROC==0) then ! calculations in master processor
                  CALL nodalAveraging(nodalValues,gaussValues,nDOF_FP, &
                                      G0XYZ,IJK,IJK_grainNodes, &
                                      error)
               endif
                                   
            elseif(recoveryMethod==rec_SPR)then ! SPR
               
               if(IPROC==0) then ! calculations in master processor
                  CALL recoverySPR(nodalValues, &
                                   gaussValues, &
                                   nDOF_FP,p_deg, &
                                   regularizeBoundaries, &
                                   localNonlocal,nFailedPatches, &
                                   error,DEBUG_nodalValues)
               endif

            elseif (recoveryMethod==rec_LSQ) then  ! LSQ
            
               ! LSQ is done in parallel
               if (IPROC/=0) then
                  if(.not.allocated(gaussValues)) allocate(gaussValues(9,NELX*NGAUSS))
                  if(.not.allocated(nodalValues)) allocate(nodalValues(9,nTotGrainNodes))
               endif

               CALL MPI_BCAST(gaussValues,nDOF_FP*NELX*NGAUSS, &
                              MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,IERROR)
                        
               ! this processor will calculate the fields of grains staGrainIdx--endGrainIdx
               ! each grain involves solving a large system (MKL is used)
               CALL recoveryLSQ(nodalValues, &
                                gaussValues, &
                                nDOF_FP, &
                                G0XYZ,IJK, &
                                grainElementIndx,grainElements,IJK_grainNodes, &
                                nGrainNodes,grainNodesIndices, &
                                error,DEBUG_nodalValues,&
                                staGrainIdx,endGrainIdx)

               CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)      
               call MPI_ALLREDUCE(MPI_IN_PLACE,nodalValues,nDOF_FP*nTotGrainNodes,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

            endif
            if (error /= 0) then
               write(*,*) 'nodal Recovery has failed. recoveryMethod',recoveryMethod,'error#',error
               CALL MPI_ABORT(IERROR)
            endif

   !----------- SUCCESSFULLY RECOVERED Fp at NODES ------------- !
   !----------- NOW CALCULATE THE NYE TENSOR on Master PROC ---- !
            if (IPROC==0) then
               ! we have the nodalValues(:,:)
               do iElem=1,NELX
                  !get the grain ID of the element
                  iGrain=grainIDelem(iElem)
                  CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
                  CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol)   
                  CALL MATINV3_UEL(elemJacMat,elemJacInv,elemVol)
                  ! store the nodal values of transpose of Fp
                  CALL getElemNodalValues(iElem,elemNodalNyeArr, &
                                          nodalValues,nDOF_FP, &
                                          IJK_grainNodes,error)
                  if (error /= 0) then
                     write(*,*) 'nodalRecovery: failed reading element nodal values'
                     CALL MPI_ABORT(IERROR)
                  endif
                  do iNode=1,4
                     elemNodalFp_T(1:3,1,iNode) = elemNodalNyeArr(1:3,iNode)
                     elemNodalFp_T(1:3,2,iNode) = elemNodalNyeArr(4:6,iNode) 
                     elemNodalFp_T(1:3,3,iNode) = elemNodalNyeArr(7:9,iNode) 
                  enddo

                  CALL calcGradTensorTET(elemJacInv,elemNodalFp_T, & 
                                                         gradFp_T)
                  ! Nye = - (curl x (Fp)^T)^T
                  CALL calcCurlTensor(-1.0*gradFp_T(:,:,:),Nye)
                  Nye(:,:) = TRANSPOSE(Nye(:,:))
                  
                  ! save to array for further processing (GNDs, error indicators etc)
                  elemNye(1:3,iElem)=Nye(1,1:3) ! ROW-MAJOR
                  elemNye(4:6,iElem)=Nye(2,1:3) ! ROW-MAJOR
                  elemNye(7:9,iElem)=Nye(3,1:3) ! ROW-MAJOR

               enddo
   ! ------------ NYE TENSOR CALCULATED: elemNye ------------ !

   ! ------------- NYE --> state variables ------------------ !
   !          (actually its not too necessary to keep this SVARS_all)
               do iElem=1,NELX
                  SVARS_all((iElem-1)*NGAUSS*NSTATE+357: &  ! ROW-MAJOR
                            (iElem-1)*NGAUSS*NSTATE+365) &  ! ROW-MAJOR
                                  = elemNye(1:9,iElem)      ! ROW-MAJOR
               enddo
   ! ------------- NYE --> OTHER PROCs ---------------------- !         
               do tProc=1,NPROCS-1
                  do iElem=NELST_proc(tProc),NELEND_proc(tProc)
                     iEl_proc = iElem - NELST_proc(tProc) + 1
                     procNye((iEl_proc-1)*9+1: & 
                             (iEl_proc-1)*9+9) = elemNye(1:9,iElem)
                  enddo
                  CALL MPI_SSEND(procNye(:), & 
                              9*NELNOX_proc(tProc),MPI_DOUBLE_PRECISION, & 
                              tProc,969,MPI_COMM_WORLD,IERROR)
                  if(IERROR.NE.MPI_SUCCESS)  & 
                        write(*,*) 'MPI_ERROR_SEND',IPROC, IERROR
               enddo
     
               ! save Nye to state variables 357-365
               ! for the next iteration step
               do iElem=NELST,NELEND
                  iEl_proc = iElem - NELST + 1
                  SVARS1((iEl_proc-1)*NGAUSS*NSTATE+357: & 
                         (iEl_proc-1)*NGAUSS*NSTATE+365) & 
                              = elemNye(1:9,iElem)
               enddo
               
            ELSE   ! IPROC.NE.0
   ! ------------- OTHER PROC RECIEVES < -- NYE  ------------ !         
               CALL MPI_RECV(procNye(:), & 
                              9*NELNOX,MPI_DOUBLE_PRECISION, & 
                              0,969,MPI_COMM_WORLD,ISTATUS,IERROR)
               if(IERROR.NE.MPI_SUCCESS)  & 
                     write(*,*) 'MPI_ERROR_RECV',IPROC, IERROR
     
               ! save Nye (recv'd) to state variables 357-365
               do iEl_proc=1,NELNOX
                  SVARS1((iEl_proc-1)*NGAUSS*NSTATE+357: & 
                         (iEl_proc-1)*NGAUSS*NSTATE+365) & 
                              = procNye((iEl_proc-1)*9+1: & 
                                        (iEl_proc-1)*9+9)
               enddo
               
            ENDIF !(IPROC.EQ.0.)
            
         endif ! do_postProc
         
         if (IPROC.eq.0.AND.do_postProc) &
            CALL task_Record(task_PostProcess,MPI_WTIME(),task_time,-1)
         !*********************************************************************!
         !                     END OF POST-PROCESSING                          !

         !*********************************************************************!
         !             (2)     EXPORT SOLUTION                                 !
         !*********************************************************************!
         ! determine if an export is needed here
         CALL doAtPeriod(FINTME,NSTEP,dataExportT,dataExport_lastT, & 
                         dataExportN,do_dataExport)
            
         if(EXPORT_this_step) do_dataExport =.TRUE.   ! do simple data export for this time step - this is a special point in the load waveform
         if(EXPORT_specialOnly .and. .not.EXPORT_this_step) then 
            do_dataExport = .false.                   ! if the user has requested to export the special load-form points only, 
         endif                                        ! do not export unless we are on one of them
         if((dabs(TOTME-FINTME).LT.EXPORT_allAtLast)) then  ! it is requested that all variables are exported
            do_dataExport = .TRUE.
            !EXPORT_delta_gamma = .TRUE.
            EXPORT_temperature = .TRUE.
            EXPORT_cauchy = .TRUE.
            EXPORT_creep = .TRUE.
            EXPORT_fp = .TRUE.
            EXPORT_nye = .TRUE.
            EXPORT_wp = .TRUE.
            !EXPORT_nodePos = .TRUE.
            EXPORT_grainAVG = .TRUE.
            EXPORT_grainRot = .TRUE.
            !EXPORT_hardening = .TRUE.
         endif
                         
         if (IPROC.eq.0.AND.do_dataExport) &
            CALL task_Start(task_Export,MPI_WTIME())
            
         if (IPROC==0) then
            
! ------------------------- RVE-average stress and strain -----------------------!
! ---------------( this is exported at all integration time steps ---------------!
            totVol = 0.D0
            totStress(:) = 0.D0
            totStrain(:) = 0.D0
            grainAVG(:,:) = 0.D0
            grainRot(:,:) = 0.D0
            totFpZZ = 0.D0
            totWP = 0.d0
            avgTemp = 0.d0
            totCrZZ = 0.D0
            IS = 0
            do iElem=1,NELX
               ! RVE-averaged quantities
               elemVol = SVARS_all(IS+300)
               iGrain = grainIDelem(iElem)
               totVol = totVol + elemVol
               totWP = totWP + elemVol*SVARS_all(IS+299)
               avgTemp = avgTemp + elemVol*SVARS_all(IS+193)
               totStress(:) = totStress(:) + elemVol*SVARS_all(IS+211:IS+216) ! Cauchy Stress: 211-216
               totStrain(:) = totStrain(:) + elemVol*SVARS_all(IS+187:IS+192) ! 187-192 -- Logarithmic Strain E=ln(U)
               totFpZZ = totFpZZ + elemVol*SVARS_all(IS+75) ! 67-75: Fp    ! WARNING ! Fp - is stored ROW-MAJOR in STATE VARIABLES
               totCrZZ = totCrZZ + elemVol*SVARS_all(IS+323) !321-326: Creep Strain
               ! grain-averaged quantities
               grainAVG(1,iGrain) = grainAVG(1,iGrain) + elemVol                                ! 1: grain volume
               grainAVG(2:4,iGrain) = grainAVG(2:4,iGrain) + elemVol*SVARS_all(IS+354:IS+356)   ! 2-4: c-axis direction
               grainAVG(5,iGrain) = grainAVG(5,iGrain) + elemVol*SVARS_all(IS+352)              ! 5: norm stress on basal 
               grainAVG(6,iGrain) = grainAVG(6,iGrain) + elemVol*SVARS_all(IS+353)              ! 6: shear stress on basal
               CALL calcE_Equiv(SVARS_all(IS+67:IS+75),Ep_Equiv)      
               grainAVG(7,iGrain) = grainAVG(7,iGrain) + elemVol*Ep_Equiv                       ! 7: equivalent plastic strain: Sqrt[2/3*||E_pl||]
               grainRot(1:9,iGrain) = grainRot(1:9,iGrain) + elemVol*SVARS_all(IS+369:IS+377)   ! grain crystal rotation matrix (STATEV:369-377, ROW-MAJOR)
               IS = IS + NSTATE
            enddo
            totStress(:) = totStress(:) / totVol
            totStrain(:) = totStrain(:) / totVol
            totFpZZ = totFpZZ / totVol
            totCrZZ = totCrZZ / totVol
            totWP = totWP / totVol
            avgTemp = avgTemp / totVol
            ! also calculate grain-average quantities: c axis, SIG_cc, TAU_c, etc.. while here ..
            do iGrain=1,nGrains  ! calculate grain averages.
               grainAVG(2:NGRAINAVG,iGrain) = grainAVG(2:NGRAINAVG,iGrain) / grainAVG(1,iGrain)
               grainRot(:,iGrain) = grainRot(:,iGrain) / grainAVG(1,iGrain)
               ! renormalize c-axis vector:
               grainAVG(2:4,iGrain) = grainAVG(2:4,iGrain) / DOT_PRODUCT(grainAVG(2:4,iGrain),grainAVG(2:4,iGrain))
            enddo
            write(981,*) FINTME, totVol, totStress(:), totStrain(:),totCrZZ,totWP
            write(405,*) FINTME, defGradMacroscopic(1,1:3), defGradMacroscopic(2,1:3), defGradMacroscopic(3,1:3)
            write(410,*) FINTME, dpMacroscopic(1:6)
            flush(981)
            flush(405)
            flush(410)
            
         endif
         
         CALL MPI_BCAST(totWP,1,MPI_DOUBLE_PRECISION, &
                        0,MPI_COMM_WORLD,IERROR)   
            
         if (IPROC.eq.0.AND.do_dataExport) then ! IMPORTANT ! -- if EXPORT_watmus, integration and data export times 
                                                !             -- must fit on fractions of 2^n.
         
            write(132,*)FINTME   ! time.out
            flush(132)
         
            if(EXPORT_nodePos) then
! ---------------------------- Node Positions ----------------------------------- !
               write(992,*) NSTEP, FINTME
               CALL printVector(GPXYZ(:),NEQ,MDIM, & 
                '',992) 
               flush(992)
               
            endif
            
! ---------------------------- the Displacement Field --------------------------- !
! --------------------------- (Essential for Plotting) -------------------------- !
            if(.not.EXPORT_watmus) then
               write(100,*) NSTEP,FINTME,TEMP(2)   ! watmus import assumes no headers -- perhaps change this in watmus gencsdat?
            endif
            IS = 0
            DO iNode=1,NX
               WRITE(100,'(3ES11.3)') VTAU(IS+1:IS+3) ! DISP.OUT
               IS = IS + 3
            ENDDO
            flush(100)
            
! ---------------------------- Domain side lengths ------------------------------ !
            if(.not.ANY(nFaceNodes(:)==0)) then
               CALL getFacePositionsAVG(domainRange,domainLen,GPXYZ, & 
                    faceNodes,nFaceNodes)
            else
               write(*,*) '**************************************************************************'
               write(*,*) ''
               write(*,*) "WARNING!!! cannot distinguish all the faces of the domain."
               write(*,*) "some surfaces tilted more than 45 degrees. some features might malfunction"
               write(*,*) ''
               write(*,*) '**************************************************************************'
            endif

! ------------------------- EXPORT FOR WATMUS -----------------------------------!
            ! IMPORTANT ! -- for WATMUS integration and data export times 
            !             -- must fit on fractions of 2^n.
            if(EXPORT_watmus) then
               
               ! export hardening variables
               
               ! DECIDE: single/multiple processors?
               ! for now, export single
               
               
            endif


            ! export grain-averaged quantities
            if(EXPORT_grainAVG) then
               write(984,*) NSTEP, FINTME, TEMP(2), 'grainVol c-axis SIGnn TAUn FpEqv'
               do iGrain=1,nGrains
                  write(984,*) grainAVG(1:NGRAINAVG,iGrain)
               enddo
               flush(984)
            endif
            if(EXPORT_grainRot) then
               write(980,*) NSTEP, FINTME, TEMP(2), 'R_cry(R-MAJOR)'
               do iGrain=1,nGrains
                  write(980,*) grainRot(:,iGrain)
               enddo
               flush(980)
            endif
            
            if (EXPORT_temperature) then
               write(801,*) NSTEP, FINTME, TEMP(2), avgTemp
               IS = 0
               do iElem=1,NELX
                  write(801,'(1F8.2)') SVARS_all(IS+193)   ! local temperature (ambient + adiabatic[if enabled])
                  IS = IS + NSTATE
               enddo
               flush(801)
            endif

            
! ------------------------- Cauchy stress ---------------------------------------!
            if(EXPORT_cauchy) then
               write(982,*) NSTEP,FINTME,TEMP(2)
               IS = 0
               do iElem=1,NELX
                  write(982,'(6ES12.3)') SVARS_all(IS+211:IS+216) ! Cauchy Stress: 211-216
                  IS = IS + NSTATE
               enddo
               flush(982)
            endif
! ------------------------- Max Basal and Prism Hardening Parameters ------------!
            if(EXPORT_hardening) then
               write(986,*) NSTEP,FINTME,TEMP(2)
               IS = 0
               do iElem=1,NELX
                  write(986,'(2F6.0)') MAXVAL(SVARS_all(IS+1:IS+3)),MAXVAL(SVARS_all(IS+4:IS+6)) ! Basal-Prismatic hardening parameters: 1-6
                  IS = IS + NSTATE
               enddo
               flush(986)
            endif
! ------------------------- hardening parameters --------------------------------!
!            write(802,*) 'Basal and Prism slip - TIME STEP', NSTEP, & 
!                         'TIME',FINTME,'TEMP',TEMP(2)
!            do iElem=1,NELX
!               write(802,*) SVARS_all((iElem-1)*NGAUSS*NSTATE+1: & 
!                                      (iElem-1)*NGAUSS*NSTATE+6)
!            enddo
! --------------------- max plastic slip, ----------------------------------------!
!                       and effective plastic slip,
!                       plastic work, w_p,
!                       normal shear stress on basal plane
!                       resolved shear stress on basal plane
!                       c-axis orientation
            if(EXPORT_delta_gamma) then
               write(983,*) NSTEP,FINTME,TEMP(2),''
               tot_gamma(:)=0.D0
               IS = 0
               do iElem=1,NELX
                  tot_gamma(1:30)=SVARS_all(IS+31:IS+60) ! HCP total gamma
                  norm_tot_gamma = DSQRT(DOT_PRODUCT(tot_gamma(1:30), & 
                                                     tot_gamma(1:30)))
                  write(983,'(8ES10.2)') SVARS_all(IS+358), & ! norm of delta gamma
                                         norm_tot_gamma, &    ! norm of total gamma
                                         SVARS_all(IS+299), & ! wp
                                         SVARS_all(IS+352), & ! norm stress on basal 
                                         SVARS_all(IS+353), & ! shear stress on basal
                                         SVARS_all(IS+354:IS+356)  ! c-axis (basal normal) direction 
                  IS = IS+NSTATE
               enddo
               flush(983)
            endif
! -----------------------damage related variables--------------------- !
            if(EXPORT_damage) then
               write(960,*) NSTEP,FINTME,TEMP(2),'Ep_equiv,wp,wpBackStr(4),min(delta_g(basals)),min(delta_g(prisms))'
               tot_gamma(:)=0.D0
               IS = 0
               do iElem=1,NELX
                  tot_gamma(1:30)=SVARS_all(IS+31:IS+60) ! HCP total gamma
                  norm_tot_gamma = DSQRT(DOT_PRODUCT(tot_gamma(1:30), & 
                                                tot_gamma(1:30)))
                  CALL calcE_Equiv(SVARS_all(IS+67:IS+75),Ep_Equiv)      

                  write(960,'(8ES10.2)')Ep_Equiv, &    ! equivalent plastic strain
                                         SVARS_all(IS+299), & ! wp
                                         !SVARS_all(IS+221:IS+223), &   ! HCP basal back stresses
                                         minval(SVARS_all(IS+1:IS+3)), &   ! delta_g, accummulated HCP-basal hardening(or softening)
                                         !SVARS_all(IS+327)            ! non-local damage energy on basal plane
                                         minval(SVARS_all(IS+4:IS+6))   ! delta_g, accummulated HCP-prism hardening(or softening)

                  IS = IS+NSTATE
               enddo
               flush(960)
            endif
            if(EXPORT_wp) then
               write(961,*) NSTEP,FINTME,TEMP(2),'wp'
               IS = 0
               do iElem=1,NELX
                  write(961,'(5ES12.4)') SVARS_all(IS+299) ! wp

                  IS = IS+NSTATE
               enddo
               flush(961)
            endif
            if(EXPORT_backstress) then
               write(962,*) NSTEP,FINTME,TEMP(2),'backStress(3xbasal)'
               IS = 0
               do iElem=1,NELX
                  write(962,'(3F8.3)') SVARS_all(IS+221:IS+223)   ! back-stress on basal-<a>

                  IS = IS+NSTATE
               enddo
               flush(962)
            endif

         ! ******* export Fp and delta_gamma ****** !
            ! update period: postProcN or postProcT
            if(EXPORT_delta_gamma) then
               write(404,*) NSTEP,FINTME,TEMP(2)    !GAMMA.out
               IS = 0
               do iElem=1,NELX
                  write(404,'(30ES11.3)')  SVARS_all(IS+31:IS+60) !delta_gamma: 31--60
                  IS=IS+NSTATE
               enddo
               flush(404)
            endif

            if(EXPORT_fp) then
               WRITE(406,*) NSTEP,FINTME,TEMP(2)    !FP.out
               IS = 0
               do iElem=1,NELX
                  write(406,'(9F12.7)')  SVARS_all(IS+67:IS+75) ! FP: 67--75   ! WARNING ! Fp - is stored ROW-MAJOR in STATE VARIABLES
                  IS=IS+NSTATE
               enddo
               flush(406)
            endif
            
            if(EXPORT_creep) then
               WRITE(407,*) NSTEP,FINTME,TEMP(2)    !creep.out
               IS = 0
               do iElem=1,NELX
                  write(407,'(6F9.5)')  SVARS_all(IS+321:IS+326) ! G-L Creep Strain: 321-326
                  IS=IS+NSTATE
               enddo
               flush(407)
            endif
            if(EXPORT_creepRate) then
               WRITE(408,*) NSTEP,FINTME,TEMP(2)    !creepRate.out
               IS = 0
               do iElem=1,NELX
                  write(408,*)  (SVARS_all(IS+323)-prevCreepStrain(iElem))/DELTME_PREV(2)  ! G-L Creep Strain: 321-326
                  prevCreepStrain(iElem) = SVARS_all(IS+323)
                  IS=IS+NSTATE
               enddo
               flush(408)
            endif
            if(EXPORT_strainRate) then
               WRITE(409,*) NSTEP,FINTME,TEMP(2)    !strainRate.out
               IS = 0
               do iElem=1,NELX
                  write(409,*)  (SVARS_all(IS+189)-prevStrain(iElem))/DELTME_PREV(2) ! Total Strain: 187-192
                  prevStrain(iElem) = SVARS_all(IS+189)
                  IS=IS+NSTATE
               enddo
               flush(409)
            endif

            ! print Nye to file
            if(EXPORT_nye) then
               WRITE(969,*) NSTEP,FINTME,TEMP(2)    !nye.out
               do iElem=1,NELX
                  write(969,'(9ES11.3)')  elemNye(1:9,iElem)
               enddo
               flush(969)
               WRITE(968,*) NSTEP,FINTME,TEMP(2)    !Nye_norm.out
               do iElem=1,NELX
                  write(968,'(1ES11.3)')  DSQRT(DOT_PRODUCT(elemNye(:,iElem),elemNye(:,iElem)))
               enddo
               flush(968)
            endif

         endif !IPROC.EQ.0.AND.do_dataExport
         
         if (IPROC.eq.0.AND.do_dataExport) &
            CALL task_Record(task_Export,MPI_WTIME(),task_time,-1)
            
      !*************************************************************************!
      !                           END OF EXPORT STAGE
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         if(IERROR.NE.0)then
            write(*,*) 'MPI BARRIER error at 1462'
            write(131,*) 'MPI BARRIER error at 1462'
            close(131)               
            CALL MPI_FINALIZE(IERROR)
            stop
         endif

!***********************************************************************************!
!        THE END OF CALCULATIONS FOR THIS TIME STEP. PROCEED TO THE NEXT.           !
!***********************************************************************************!

         NSTEP=NSTEP+1

         IF(NSTEP.GT.NOSTP)IFLAG_OVER=1
         
         IF(DABS(FINTME-TOTME).LT.1D-15.OR.FINTME.GT.TOTME)IFLAG_OVER=2
         
         if (endSimulationAtWp > 0.d0 .and. &
             totWP >= endSimulationAtWp) then
         
            IFLAG_OVER = 2
             
         endif

         if(N_UPDATE.EQ.0) then
            write(*,*) 'N_UPDATE 0, end of loop'
            CALL MPI_FINALIZE(IERROR)
            stop
         endif
      ENDDO
      
701   CONTINUE ! program termination point
      
      ! End of calculation
      
      ! EXPORT SOLUTION for Continuation
      if (IPROC==0 .and. exportSolution) then
         write(*,'(A)',advance='no') ' exporting solution for continuation... '
         CALL execute_command_line('mkdir continue')
                  
         open(102,FILE='continue/problem.out')
         write(102,*) NSTATE
         write(102,*) NELX
         write(102,*) NX
         write(102,*) lastConvergedTime
         close(102)
         
         open(102,FILE='continue/disp.out')
         write(102,*) NSTEP,lastConvergedTime,TEMP(2)
         IS = 0
         DO iNode=1,NX
            write(102,*) VT(IS+1:IS+3) ! DISP.OUT
            IS = IS + 3
         ENDDO
         close(102)

         open(102,FILE='continue/svars.out')
         do iElem=1,NELX
            write(102,*) SVARS_all((iElem-1)*NGAUSS*NSTATE+1:iElem*NGAUSS*NSTATE)
         enddo
         close(102)
         write(*,*) 'done.'
      endif

      ! record total time
      if(DEBUG_timer) CALL task_Record(task_PROGRAM,MPI_WTIME(),T_ELAPSED,printTasksTo)
      
      ! extract CPU-time allocation report
      if(DEBUG_timer) CALL task_Report(printTasksTo)

      DEALLOCATE(VAL,RESID,IBINDX,IROWPTR)
      DEALLOCATE(RESID_LOC,IBINDX_LOC,IROWPTR_LOC)
      
      deallocate(VAL_rxn,ICOLPTR_rxn,IBINDX_rxn)
      deallocate(EQ_uBC_proc,uBC_proc)

      IF(ICYC.EQ.1)THEN
         DEALLOCATE(TSEG)
      ENDIF

      deallocate(F,F0,F_ext,F_uBC,F_BFGS)
      deallocate(incVTAU_init,DVTAU0,DVTAU)
      deallocate(VTAU,VTAU0,VT)
      deallocate(u_prev)
      
      ! arrays for BFGS update
      deallocate(VAR,WAR)
      deallocate(V,W)
      deallocate(F_prev)
      
      ! deallocate variables for dvtau0 improvement
      deallocate(U_BCs_prev,U_BCs,dU_BCs)

      ! *********** added by deniz **********
      ! deallocate gradient recovery arrays
      
      
      if (allocated(gaussValues)) deallocate(gaussValues)
      if (allocated(nodalValues)) deallocate(nodalValues)
      if (allocated(nodalNye)) deallocate(nodalNye)
      if(IPROC.EQ.0) then
         !deallocate(errorIndicators)
         !deallocate(elementsToRefine)
         deallocate(elemNye)
      endif
      
      deallocate(elemFaceNeighbors)
      !deallocate(elemEdgeNeighbors)

      ! element pairs for non-local damage/regularisation
      if (IPROC==0.and.damage_nonlocal) &
         deallocate(elemNonLocalList,elemNonLocalCount,elemNonLocalWeights)
      if (damage_nonlocal) deallocate(procNonLocalDamage)

      CALL grains_Destruct()
      deallocate(procNye)
      if(allocated(SVARS_all)) deallocate(SVARS_all)   ! state variables array (only allocated in master proc)
      if(allocated(prevStrain)) deallocate(prevStrain)
      if(allocated(prevCreepStrain)) deallocate(prevCreepStrain)
      deallocate(grainAVG)    ! array for grain-averaged quantities
      deallocate(grainRot)    ! array for grain rotations
      ! **************************************
      
      IF(IPROC.EQ.0)THEN
         WRITE(131,*)'TIME OF OPERATION=',T_ELAPSED,'SECS'
         WRITE(131,*)'PROGRAM ENDS'
         CLOSE(131)
         CLOSE(100)
         
         if(DEBUG_nodalValues) CLOSE(987) ! nodal value recovery data
         if(DEBUG_nodalValues) CLOSE(997) ! nodal value recovery debug
         if(DEBUG_timer) CLOSE(994) ! CPU-time information
         if(DEBUG_linesearch) CLOSE(995) ! debug line search
         if(DEBUG_general) close(990)  ! debug/error log file
         if(DEBUG_RESID) close(970)    ! debug/export residual
         CLOSE(981) ! stress, strain (RVE-averaged)
         close(405) ! macroscopic deformation gradient
         close(410) ! macroscopic rate of plastic deformation (current config)
         if(EXPORT_temperature.OR.EXPORT_allAtLast.GT.0.0) CLOSE(801) ! local temperature
         if(EXPORT_cauchy     .OR.EXPORT_allAtLast.GT.0.0) CLOSE(982)! Cauchy/2nd PK ??
         if(EXPORT_damage     .OR.EXPORT_allAtLast.GT.0.0) CLOSE(960)! damage related variables
         if(EXPORT_wp         .OR.EXPORT_allAtLast.GT.0.0) CLOSE(961)! dissipated / plastic energy
         if(EXPORT_backstress) CLOSE(962)! basal back stress
         if(EXPORT_delta_gamma) CLOSE(983)! plastic slip, plastic work, c-axis orientation
         if(EXPORT_delta_gamma) CLOSE(404)! delta gamma
         if(EXPORT_grainAVG   .OR.EXPORT_allAtLast.GT.0.0) CLOSE(984)! grain-averaged information
         if(EXPORT_grainRot   .OR.EXPORT_allAtLast.GT.0.0) CLOSE(980)! grain-averaged information
         if(EXPORT_hardening) CLOSE(986)! basal/prismatic hardening parameters
         if (EXPORT_nye    .OR.EXPORT_allAtLast.GT.0.0) CLOSE(969) ! the Nye tensor
         if (EXPORT_nye    .OR.EXPORT_allAtLast.GT.0.0) CLOSE(968) ! the Nye tensor norm
         if (EXPORT_nodePos) CLOSE(992) ! node positions
         
         if (EXPORT_fp     .OR.EXPORT_allAtLast.GT.0.0) CLOSE(406) ! Fp
         if (EXPORT_creep  .OR.EXPORT_allAtLast.GT.0.0) CLOSE(407) ! Creep strain
         if (EXPORT_creepRate) CLOSE(408)  ! local rate of creep strain
         if (EXPORT_strainRate) CLOSE(409) ! local rate of total strain
         !CLOSE(802) ! ratio of contributions to hardening from SSDs and GNDs 
      ENDIF
      
      if(DEBUG_RESID) close(542)    ! debug/export processor specific residual/search direction/rxn forces
      if(DEBUG_UMAT) close(543)  ! debug/export processor specific UMAT debug info

      IF(IPROC.EQ.0)THEN
         CLOSE(132)
      ENDIF

      ! from superLU user guide page 71.
      ! DEALLOCATE THE STORAGE ALLOCATED BY SUPERLU_DIST
      CALL F_PSTATFREE(STAT)
      CALL F_DESTROY_SUPERMAT_STORE_DIST(A)
      CALL F_SCALEPERMSTRUCTFREE(ScalePermstruct)
      CALL F_DESTROY_LU(N, grid, LUstruct)
      CALL F_LUSTRUCTFREE(LUstruct)
      CALL GET_SUPERLU_OPTIONS(SLUoptions, SOLVEINITIALIZED=INIT)
      IF (INIT == YES) THEN
         CALL F_DSOLVEFINALIZE(SLUoptions, SOLVEstruct)
      ENDIF

      ! release superLU process grid
      CALL F_SUPERLU_GRIDEXIT(grid)
      
      !     DEALLOCATE THE C STRUCTURES POINTED TO BY THE FORTRAN HANDLES
      CALL F_DESTROY_GRIDINFO_HANDLE(grid)
      CALL F_DESTROY_OPTIONS_HANDLE(SLUoptions)
      CALL F_DESTROY_SCALEPERM_HANDLE(ScalePermstruct)
      CALL F_DESTROY_LUSTRUCT_HANDLE(LUstruct)
      CALL F_DESTROY_SOLVESTRUCT_HANDLE(SOLVEstruct)
      CALL F_DESTROY_SUPERMATRIX_HANDLE(A)
      CALL F_DESTROY_SUPERLUSTAT_HANDLE(STAT)
      
      ! terminate MPI execution environment      
      CALL MPI_FINALIZE(IERROR)

      END PROGRAM

!*********************************************************

      SUBROUTINE READGM(GXYZ,IMID,IJK,MDIM,NFACES,NEDGES)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION GXYZ(MAXCRD),IMID(3*MAXEL)
      DIMENSION IJK(MNELX)
      DIMENSION IJKNEL(MAXELN),XC(3)
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV
      integer :: IERROR

      CALL READSTR(LR)
      READ(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD
      !deniz
      IF(NDF.GT.NDOFELX)THEN
         IF(IPROC.EQ.0) write(*,*) 'INSUFFICIENT MEM-NDOFELX'
         IF(IPROC.EQ.0) write(131,*)'INSUFFICIENT MEM-NDOFELX'
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF

      CALL READSTR(LR)
      READ(LR,*) NX
      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         IF(IPROC.EQ.0) write(*,*) 'Increase MAXNODE to ', NX
         IF(IPROC.EQ.0) write(131,*)'Increase MAXNODE to ', NX
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF
      NEQ=NDF*NX
      IF(NEQ.GT.MDOFX)THEN
         IF(IPROC.EQ.0) write(131,*)'INSUFFICIENT MEMORY-MDOFX'
         IF(IPROC.EQ.0) write(*,*)'INSUFFICIENT MEMORY-MDOFX'
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF
      IF(NX*MDIM.GT.MAXCRD)THEN
         IF(IPROC.EQ.0) write(131,*) ' INSUFFICIENT MEM -MAXCRD'
         IF(IPROC.EQ.0) write(*,*) 'INSUFFICIENT MEM -MAXCRD'
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF
      IF(MDIM.GT.MAXDIM)THEN
         IF(IPROC.EQ.0) WRITE(131,*)'INSUFFICIENT MEM-MAXDIM'
         IF(IPROC.EQ.0) WRITE(*,*)'INSUFFICIENT MEM-MAXDIM'
         CALL MPI_FINALIZE(IERROR)
         STOP
      ENDIF

      DO I=1,NX
         READ(LR,*) NNUM,(XC(II),II=1,MDIM)
         DO II=1,MDIM
            GXYZ(MDIM*(NNUM-1)+II)=XC(II)
         ENDDO
      ENDDO
!     ( READ ELEMENT CONNECTIVITY )

      CALL READSTR(LR)

      !-----------IF 1 THEN BRICK, IF 2 THEN TET-------------
      IF(IBRICKORTET.EQ.1)THEN

         NGAUSS=8
         NODE=8
         !modified: deniz
         !used by SPR
         NFACES=6
         NEDGES=12
         
      ELSEIF(IBRICKORTET.EQ.2)THEN

         NGAUSS=1
         NODE=4
         !modified: deniz
         !used by SPR
         NFACES=4
         NEDGES=6
         
      ENDIF

      MPE=1
      IEP=1
      READ(LR,*) NELX
      IF(NELX.GT.MAXEL)THEN
         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'INSUFFICIENT MEMORY -MNEL'
         ENDIF
         STOP
      ENDIF
      JJ=0
      DO N=1,NELX
         READ(LR,*) NEL,(IJKNEL(J),J=1,NODE)
         IMID(NEL)=NODE
         IMID(NELX+NEL)=MPE
         IMID(2*NELX+NEL)=IEP
         DO J=1,NODE
            JJ=JJ+1
            IJK(JJ)=IJKNEL(J)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE

!************************************************************
      ! NEW:
      SUBROUTINE readMaterial(props,nprops,strPropsFile,importedMaterial)
      implicit none
      include 'PARDIS.H'
      real(8), intent(out) :: props(MAXPROPS)
      integer, intent(out) :: nprops
      character(len=*), intent(in) :: strPropsFile
      logical, intent(out) :: importedMaterial
      
      logical :: fileExists, endOfFile
      integer :: nRows, posProp, error
      
      importedMaterial =.false.
      
      
      inquire(file=trim(strPropsFile),exist=fileExists)
      if (.NOT.fileExists) return
      
      open(105,file=trim(strPropsFile))
      
      read(105,*) nprops
      if (nprops > MAXPROPS) then
         write(*,*) 'increase MAXPROPS to ',nprops
         stop
      endif
      nRows = 0
      posProp = 0
      endOfFile = .false.
      do while (posProp < nprops .and.(.not.endOfFile))
         read(105,*,iostat=error) props(posProp+1:posProp+8)
         if(error.ne.0) then 
            endOfFile = .true.
         else
            posProp = posProp + 8
         endif
      enddo
      
      if(posProp.ne.0.and.posProp > nprops-8) importedMaterial = .true.
      
      END SUBROUTINE
      
      ! OLD:
      SUBROUTINE READMT(NGROUP,PROMAT2D,NPROPSG,NSTATE)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION PROMAT2D(MAXPROPS,MAXGRP)
      DIMENSION NPROPSG(MAXGRP)
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV

      CALL READSTR(LR)
      READ(LR,*)NGROUP
      IF(NGROUP.GT.MAXGRP)THEN
         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'INSUFFICIENT MEM -MAXGRP'
         ENDIF
         STOP
      ENDIF

      READ(LR,*)(NPROPSG(N),N=1,NGROUP)
      DO I=1,NGROUP
         IF(NPROPSG(I).GT.MAXPROPS)THEN
            IF(IPROC.EQ.0)THEN
               WRITE(131,*)'INSUFFICIENT MEM - MAXPROPS'
            ENDIF
            STOP
         ENDIF
      ENDDO

      DO N=1,NGROUP
         MATGRP=NPROPSG(N)
         NCOUNT=MATGRP/8
         NREMAIN=MOD(MATGRP,8)
         NSHIFT=0
         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'NUM OF PROPS: ',MATGRP
         ENDIF
         DO J=1,NCOUNT
            READ(LR,*)(PROMAT2D(NSHIFT+I,N),I=1,8)
            NSHIFT=NSHIFT+8
         ENDDO
         READ(LR,*)(PROMAT2D(NSHIFT+I,N),I=1,NREMAIN)
      ENDDO
      READ(LR,*)NSTATE
      IF(NSTATE.GT.MSV)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY -MSV'
         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'INSUFFICIENT MEMORY -MSV'
         ENDIF
         STOP
      ENDIF
      RETURN
      END SUBROUTINE

                    
      ! evaluate user-defined functions for displacement BCs if applicable
      ! (const strain rate deformation, strain-history, ...)
      subroutine calcFunctionDisp(time,loadType_out,functionDisp,strain,macroF_tau, &
                                  uCorner_tau,ksiCorner_init_out)
      use options, only : useLogStrainRate
      implicit none
      include 'PARDIS.H'   ! need this just because MAX_HISTORY_N
      
      real(8), intent(in) :: time(2)
      integer, intent(out):: loadType_out
      real(8), intent(out):: functionDisp,strain,macroF_tau(3,3)
      real(8), intent(out):: uCorner_tau(3,8),ksiCorner_init_out(3,8)
      
      ! locals/common blocks
      real(8) :: du(3)
      real(8) :: tperiod,tramp,t_nodwell,t_dwell,samp,smin, & 
                      uamp,umin
      real(8) :: sload,st_rate
      integer :: icyc
      real(8), parameter :: Identity(3,3) = &
         reshape((/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /), (/3,3/))
      real(8), parameter :: PI=3.14159265359D0
      
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin, & 
                      uamp,umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc

      ! deniz
      integer :: loadType      ! 1: creep, 2: CSR, 3: stress cyclic, 4: displ cyclic. 5/6: load/displ. history
      integer*4 :: padding
      real(8) :: P_cyclic_max  ! in MPa
      real(8) :: P_cyclic_min  ! in MPa
      real(8) :: P_cyclic_period ! in seconds
      integer :: P_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: P_dwell_ramptime
      integer :: P_history_repeat
      real(8) :: P_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: P_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: P_history_N    ! number of load-time points
      
      real(8) :: U_cyclic_max  !
      real(8) :: U_cyclic_min  !
      real(8) :: U_cyclic_period ! in seconds
      integer :: U_dwell_load    
      real(8) :: U_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: U_history_repeat    ! 1: repeat provided u-time history, 0: keep u at the last provided value
      real(8) :: U_history_values(MAX_HISTORY_N) ! provide u-time points. 
      real(8) :: U_history_times(MAX_HISTORY_N)  ! program will interpolate at integration time steps
      integer :: U_history_N    ! number of U-time points
      
      integer :: load_axis,iStressLateral,iStressMain
      real(8) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6) 
      logical :: specifiedLoadFrame
      logical :: strainControlledSimulation
      real(8) :: biaxial_StressRatioAngle      
      real(8) :: biaxial_StrainRate
      real(8) :: deformMultiplierTensionCompression
      real(8) :: macroStrain_t(6),macroStrain_tau(6)
      
      integer :: strain_history_repeat    ! 1: repeat provided strain-time history, 0: keep u at the last provided value
      real(8) :: strain_history_values(6,MAX_HISTORY_N) ! provide strain-time points. 
      real(8) :: strain_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: strain_history_N    ! number of strain-time points
      
      real(8) :: ksiCorner_init(3,8)  ! initial positions of corner nodes
      real(8) :: uCorner_history_values(3,8,MAX_HISTORY_N) ! provide displacement-time histories of corners of the cubic domain
      real(8) :: uCorner_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: uCorner_history_N           ! number of displacement-time points
      integer :: uCorner_remove_gradient     ! remove gradient-modes (trapezoids) from the deformation

      integer :: mark_for_export(MAX_HISTORY_N) ! mark time step for data-export. 0:default behavior, 1: force simple export, 2: force export all.
      integer :: export_special                 ! mark the special time points of cyclic loadings for data export

      common/load_cond/loadType,P_history_N,U_history_N,strain_history_N,uCorner_history_N, & 
                       P_history_repeat,U_history_repeat,P_dwell_load,U_dwell_load, &
                       P_dwell_ramptime,U_dwell_ramptime,P_cyclic_max, & 
                       P_cyclic_min,P_cyclic_period,P_history_values, & 
                       P_history_times,U_cyclic_max,U_cyclic_min, & 
                       U_cyclic_period,U_history_values,U_history_times, &
                       macroStrain_t, macroStrain_tau, load_axis,strainControlledSimulation, &
                       biaxial_StressRatioAngle,biaxial_StrainRate,deformMultiplierTensionCompression, &
                       iStressLateral,iStressMain,specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv, &
                       strain_history_values,strain_history_times,strain_history_repeat, &
                       uCorner_history_values,uCorner_history_times,uCorner_remove_gradient, &
                       ksiCorner_init,mark_for_export,export_special

      real(8) :: domainRange0(3,2), domainLen0(3)
      real(8) :: domainRange(3,2), domainLen(3)
      COMMON/DOMAINGEO/domainRange0,domainLen0,domainRange,domainLen
      
      real(8) :: umax
      real(8) :: U_cyclic_amp, U_cyclic_mean, t_cyc, t2_cyc, tramp_cyc
      real(8) :: w1,w2
      logical :: found
      integer :: iTime, i_cyc
      
      loadType_out = loadType
      
      if(loadType==2) then    ! --------- CSR loading

         if (useLogStrainRate) then
            strain = dexp(st_rate*time(1))-1.d0
         else
            strain = st_rate*time(1)
         endif
         
      elseif(loadType==4) then    ! --------- displacement-controlled cyclic loading
         U_cyclic_amp = 0.5D0*(U_cyclic_max-U_cyclic_min)
         U_cyclic_mean = 0.5D0*(U_cyclic_max+U_cyclic_min)
         tramp_cyc = U_cyclic_mean /  & 
                     (U_cyclic_amp*2*PI/U_cyclic_period)
  
         if (time(1).LT.tramp_cyc) then
            functionDisp = time(1)/tramp_cyc*U_cyclic_mean
         else
            t_cyc = time(1) - tramp_cyc
            functionDisp = U_cyclic_mean +  & 
                     U_cyclic_amp*SIN(2*PI*t_cyc/U_cyclic_period)
         endif
      
      elseif(loadType==6) then    ! --------- disp-time provided
         t_cyc = time(1)
         if(U_history_repeat.EQ.1) then
            t_cyc = DMOD(time(1),U_history_times(U_history_N))
         endif
         
         found=.FALSE.
         do iTime=1,U_history_N
            if(t_cyc.LT.U_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               functionDisp = U_history_values(1)
            else
               w1 = (U_history_times(iTime)-t_cyc) / & 
            (U_history_times(iTime)-U_history_times(iTime-1))
               w2 = (t_cyc-U_history_times(iTime-1)) / & 
            (U_history_times(iTime)-U_history_times(iTime-1))
               functionDisp = w1*U_history_values(iTime-1) & 
                    + w2*U_history_values(iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            functionDisp = U_history_values(U_history_N)
         endif
         
      elseif(loadType==8) then    ! --------- strain-time provided

         t_cyc = time(1)
         if(strain_history_repeat.EQ.1) then
            t_cyc = DMOD(time(1),strain_history_times(strain_history_N))
         endif
         
         found=.FALSE.
         do iTime=1,strain_history_N
            if(t_cyc.LT.strain_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               macroStrain_tau(:)=strain_history_values(:,1)
            else
               w1 = (strain_history_times(iTime)-t_cyc) / & 
            (strain_history_times(iTime)-strain_history_times(iTime-1))
               w2 = (t_cyc-strain_history_times(iTime-1)) / & 
            (strain_history_times(iTime)-strain_history_times(iTime-1))
               macroStrain_tau(:) = w1*strain_history_values(:,iTime-1) & 
                         + w2*strain_history_values(:,iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            macroStrain_tau(:)=strain_history_values(:,strain_history_N)
         endif

         if (useLogStrainRate) then
            call getFfromE_log(macroF_tau,macroStrain_tau)
         else
            call getFfromE_GL(macroF_tau,macroStrain_tau)
         endif
         

      elseif(loadType==18) then    ! --------- displacement-time provided for the corner nodes

         t_cyc = time(1)
         
         found=.FALSE.
         do iTime=1,uCorner_history_N
            if(t_cyc.LT.uCorner_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               uCorner_tau(:,:)=uCorner_history_values(:,:,1)
            else
               w1 = (uCorner_history_times(iTime)-t_cyc) / & 
            (uCorner_history_times(iTime)-uCorner_history_times(iTime-1))
               w2 = (t_cyc-uCorner_history_times(iTime-1)) / & 
            (uCorner_history_times(iTime)-uCorner_history_times(iTime-1))
               uCorner_tau(:,:) = w1*uCorner_history_values(:,:,iTime-1) & 
                         + w2*uCorner_history_values(:,:,iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            uCorner_tau(:,:)=uCorner_history_values(:,:,uCorner_history_N)
         endif
         
         ksiCorner_init_out = ksiCorner_init

      elseif(loadType==9 .or. &  ! --------- VM-strain controlled, biaxial stress loading
             loadType==10) then  ! --------- VM-strain controlled, uniaxial loading along arbitrary direction

         if (useLogStrainRate) then
            call getFfromE_log(macroF_tau,macroStrain_tau)
         else
            call getFfromE_GL(macroF_tau,macroStrain_tau)
         endif
         
      endif
      
      end subroutine

!*******************************************************************
      ! output: V(iDOF)<--  imposed displacements, 
      !                     for all displacement imposed DOFs.
      !                     other DOFs are untouched.
      SUBROUTINE DISINC(V,VBCS1,NBC1,DELT,IDF,DIR_U,KINC,KSTEP, & 
                   G0XYZ,TIME,NDF,MCRD)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION V(NEQ),VBCS1(MBOUND1),NBC1(MBOUND1),U(3),TIME(2)
      DIMENSION IDF(MBOUND1),G0XYZ(MAXCRD)
      INTEGER:: DIR_U(MBOUND1)   ! stores the directions +1/-1 of the displacement imposed by the disp() function. see input flag IFUNC
      real(8) :: nodePosRel(3)
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      real(8) :: bendingMultiplier
      real(8) :: functionDisp,strain,macroF_tau(3,3)
      real(8) :: uCorner_tau(3,8),ksiCorner_init(3,8)
      integer :: loadType

      integer::NX,NELX,NEQ,NPLANE,NODE,NDFF,NGAUSS
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDFF,NGAUSS
      !deniz - line added:
      real(8) :: domainCenter(3)
      real(8) :: domainRange0(3,2), domainLen0(3)
      real(8) :: domainRange(3,2), domainLen(3)
      COMMON/DOMAINGEO/domainRange0,domainLen0,domainRange,domainLen
      ! domain range < -----------
      
      domainCenter(:) = 0.5D0*(domainRange0(:,1) + domainRange0(:,2))
      
      ! evaluate user-defined functions for displacement BCs if applicable
      ! (const strain rate deformation, strain-history, ...)
      call calcFunctionDisp(time,loadType,functionDisp,strain,macroF_tau, &
                            uCorner_tau,ksiCorner_init)
      
      NDF1=NDF+1
      DO N=1,NBOUD1
         NODE1=NBC1(NDF1*(N-1)+1)   ! ID of the constrained node
         DO  I=1,NDF
            INDEX=NBC1(NDF1*(N-1)+1+I)
            IF(INDEX.EQ.1) THEN
               NN=NDF*(NODE1-1)+I   ! EQN number of the constrained displacement
!               NN=IEQNO(NN)
               IF(IDF(N).EQ.0)THEN  ! displacements provided. if some displ. BCs are imposed, 
                  N1=NDF*(N-1)+I    ! they are GRADUALLY increased to prescribed values 
                  V(NN)=V(NN)+VBCS1(N1)*DELT ! throughout the simulation. (see this line. DELT=time increment)
               ELSEIF(IDF(N).EQ.1)THEN  ! displacement imposed by an internal function (DISP)
                  ICN=MCRD*(NODE1-1)
                  DO II=1,MCRD
                     nodePosRel(II)=G0XYZ(ICN+II) - domainCenter(II)
                  ENDDO
                  ! DISP modified by deniz: added 3rd argument: disp DOF
                  call disp(u,I,nodePosRel,loadType,functionDisp,strain, &
                            macroF_tau,uCorner_tau,ksiCorner_init,domainLen0)
                  V(NN)=U(1)*SIGN(1,DIR_U(N))   ! deniz - modify the sign of the disp()-imposed displacement by the sign of IFUNC provided in the input file
               ELSEIF(IDF(N).EQ.4)THEN  ! added by deniz: impose bending
                  ICN=MCRD*(NODE1-1)
                  DO II=1,MCRD
                     nodePosRel(II)=G0XYZ(ICN+II) - domainCenter(II)
                  ENDDO
                  ! predetermined bending direction according to
                  ! the direction of imposed displacement
                  J=0
                  if(I.EQ.1) J=3
                  if(I.EQ.2) J=3
                  if(I.EQ.3) J=2
                  bendingMultiplier =  & 
                     nodePosRel(J)/domainLen0(J) - 0.5D0
                  call disp(u,I,nodePosRel,loadType,functionDisp,strain, &
                            macroF_tau,uCorner_tau,ksiCorner_init,domainLen0)
                  V(NN)=U(1)*bendingMultiplier
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
      
      ! calculates the initial guess of macroscopic strain increment for biaxial (stress) loading
      
      ! output: DeltaMacroStrain_initguess
      ! definition and usage:
      ! macroStrain_initguess = macroStrain_t + DeltaMacroStrain_initguess
      subroutine initialGuessForMacroscopicStrain(DeltaMacroStrain_initguess, &
                 DDSDDE_RVE,incrVonMises_imposed,biaxial_StressRatioAngle, &
                 iStressMain,iStressLateral,deformMultiplierTensionCompression, &
                 specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv)
               
      implicit none   
      
      real(8), intent(out):: DeltaMacroStrain_initguess(6)
      real(8), intent(in) :: DDSDDE_RVE(6,6)
      real(8), intent(in) :: incrVonMises_imposed,biaxial_StressRatioAngle
      integer, intent(in) :: iStressMain,iStressLateral
      real(8), intent(in) :: deformMultiplierTensionCompression
      logical, intent(in) :: specifiedLoadFrame
      real(8), intent(in) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6)
      
      !locals
      real(8) :: vmBiaxialStressLoading
      real(8) :: strainBiaxialStressLoading(6), unitBiaxialStressLoading(6)
      integer :: iError
      
      ! unit loading expressed in the load frame
      unitBiaxialStressLoading = 0.d0
      unitBiaxialStressLoading(iStressMain) = +1.d0*cos(biaxial_StressRatioAngle)  ! main
      unitBiaxialStressLoading(iStressLateral) = +1.d0*sin(biaxial_StressRatioAngle)  ! lateral
      
      ! apply Tension/Compression factor
      unitBiaxialStressLoading = unitBiaxialStressLoading * deformMultiplierTensionCompression
   
      ! (if loading is in arbitrary direction)
      ! rotate the loading condition to the desired load frame
      !----------------------------------------------------------
      if (specifiedLoadFrame) then
         unitBiaxialStressLoading = MATMUL(R_loadVoigt,unitBiaxialStressLoading)
      endif

      CALL solveAXb_LU(DDSDDE_RVE,6,unitBiaxialStressLoading,iError)
      if (ierror /= 0) then
         write(*,*) 'lapack error while solving for biaxial initial guess'
         CALL printMatrix(DDSDDE_RVE,6,6,'DDSDDE_RVE',0)
         stop
      endif
      
      strainBiaxialStressLoading = unitBiaxialStressLoading
!         strainBiaxialStressLoading = solveMatrixEquation(DDSDDE_RVE,unitBiaxialStressLoading)

      call getVonMisesStrainVoigt_MathStrain(vmBiaxialStressLoading,strainBiaxialStressLoading)
      
      DeltaMacroStrain_initguess = strainBiaxialStressLoading * incrVonMises_imposed / vmBiaxialStressLoading
      
      end subroutine

!***********************************************

      SUBROUTINE READSTR(LR)
      IMPLICIT REAL*8(A-H,O-Z)
      !old version:
      !DIMENSION FLAG(20)
      !READ(LR,505) (FLAG(I),I=1,20)
!505   FORMAT(20A4)

! **** modified by deniz
      ! this code skips 1 line of text in the input file.
      character(150)::dummyString
      read(LR,'(A)') dummyString

      RETURN
      END SUBROUTINE
!*************************************************

      SUBROUTINE INITABQ(NSTATE,MDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION TIME1(2),SVARS1(MAXSV),SVARS2(MAXSV)
      COMMON/ABQ/TIME1,NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,SVARS1,MCRD1
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/DOF/NDOF
      COMMON/EL_PAR/NELST,NELEND,NELNOX
      COMMON/STATEUPDATE/SVARS2

      MCRD1=MDIM
      NDOF=NEQ/NX
      NUMN=4
      IF(MCRD1.EQ.2)NUMN=4
      NDOFEL1=NDOF*NUMN
      MDLOAD1=1
      NPREDF1=1
      TIME1(1)=0.D0
      TIME1(2)=0.D0
      NSVARS1=NSTATE
      DO I=1,MAXSV
         SVARS1(I)=0.D0
         SVARS2(I)=0.D0
      ENDDO
      RETURN
      END SUBROUTINE

!******************************************************

      SUBROUTINE TMSTEP(TOTME,XIDT,XMXDT,XMINDT,NOSTP)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION FLAG(20)
      READ(LR,505) (FLAG(I),I=1,20)
      READ(LR,*) TOTME,XIDT,XMXDT,XMINDT,NOSTP
505   FORMAT(20A4)
      RETURN
      END SUBROUTINE

!*****************************************************************


! deniz -- read thermal loading information
      SUBROUTINE READTEMP(temp_loading_type,temp_stress_free, & 
                          temp_cyclic_amplitude, temp_cyclic_period, & 
                          temp_history_values, temp_history_times, & 
                          temp_history_N, & 
                          info)
     
      IMPLICIT NONE

      INCLUDE 'PARDIS.H'
      
      integer, intent(out) :: temp_loading_type ! 1: temperature time history given
                                                ! 2: cyclic thermal loading around initial temp
                                                ! 3: cyclic thermal loading ABOVE initial temp
      real(8), intent(out) :: temp_stress_free  ! initial/stress-free temperature
      real(8), intent(out) :: temp_cyclic_amplitude   ! in Kelvins
      real(8), intent(out) :: temp_cyclic_period      ! in seconds
      real(8), intent(out) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8), intent(out) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
      integer, intent(out) :: temp_history_N                          ! for integration time steps
      integer, intent(out) :: info 
      ! info=1: temp loading imported, info=0: error, or no temp. loading
      
      ! locals
      integer :: iTemp
      logical :: fileExists
      CHARACTER(len=32) :: fileName
      
      fileName='thermal_loading.inp'
      
      info = 0
            
      INQUIRE(file=fileName,exist=fileExists)
      if (.NOT.fileExists) then
         ! no thermal loading info. default conditions:
         ! no temperature change throughout simulation
         temp_loading_type=1
         temp_stress_free=273
         temp_history_values(1)=273
         temp_history_times(1)=0
         temp_history_N=1
         return
      end if
      
      OPEN(201,file=fileName)

      CALL READSTR(201)
      READ(201,*) temp_loading_type
      
      CALL READSTR(201)
      READ(201,*) temp_stress_free

      CALL READSTR(201)
      CALL READSTR(201)
      if (temp_loading_type.EQ.2.OR. & 
          temp_loading_type.EQ.3) then       ! CYCLIC THERMAL LOADING
      
         READ(201,*) temp_cyclic_amplitude, temp_cyclic_period
         
      elseif (temp_loading_type.EQ.1) then   ! TEMP. HISTORY PROVIDED
      
         READ(201,*) temp_history_N
         
         if(temp_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N to ',temp_history_N
            STOP
         endif
         
         do iTemp = 1,temp_history_N
            READ(201,*) temp_history_times(iTemp), & 
                        temp_history_values(iTemp)
         enddo
         if(temp_history_N.EQ.0)then
            ! default: T = T0 (const)
            temp_history_N = 1
            temp_history_times(1) = 0.0
            temp_history_values(1) = temp_stress_free
         end if

      else
      
         ! unexpected input. adopt default conditions:
         ! no temperature change throughout simulation
         temp_loading_type=1
         temp_stress_free=273.0
         temp_history_values(1)=273.0
         temp_history_times(1)=0.0
         temp_history_N=1
         
      endif
      
      ! success
      info = 1
      CLOSE(201)
      
      END SUBROUTINE

      ! deniz - comment
      ! NBC1 : in blocks of 4 integers, indicates whether DOFs are fixed.
      !        stores: nodeID, DOF1_fixed?, DOF2_fixed?, DOF3_fixed?, nodeID, DOF1_fixed?, DOF2_fixed?, ...
      ! VBC1 : in blocks of 3 reals,
      !        stores: u_DOF1, u_DOF2, u_DOF3, u_DOF1, u_DOF2, u_DOF3, ...
      ! IBC(EQN) : indicates whether DOF corresp to the EQN number is fixed (i.e. a displ. BC is defined).
      ! DIR_U    : stores the directions +1/-1 of the displacement imposed by the disp() function. see input flag IFUNC
      SUBROUTINE READBC(NBC1,VBC1,IDF,DIR_U,NBC2,VBC2,NBC3,VBC3,NDF,IBC,IBCT)
      IMPLICIT NONE

      INCLUDE 'PARDIS.H'

      REAL(8):: VBC1(MBOUND1),VBC2(MBOUND2)
      REAL(8):: VBC3(MBOUND3)
      REAL(8):: DK(MAXDIM),PRESSB(MAXEL)
      INTEGER:: NBC3(MBOUND3),IK(MAXDIM)
      INTEGER:: NBC1(MBOUND1),NBC2(MBOUND2)
      INTEGER:: NNSET(MAXNSET),NLIST(MAXNS,MAXNSET)
      INTEGER:: IDF(MBOUND1), IBC(MDOFX),IBCT(MDOFX)
      INTEGER:: DIR_U(MBOUND1)   ! stores the directions +1/-1 of the displacement imposed by the disp() function. see input flag IFUNC
      INTEGER:: IPRESS_SET(MAXPRESS),IPRESS_SF(MAXPRESS)
      INTEGER:: ILF(MAXEL),IFACE_TET(12),IFACE_BRICK(24)
      INTEGER:: IBELEM(MAXEL)
      INTEGER:: NELSET(MAXELSET),NELLIST(MAXELS,MAXELSET)
      INTEGER:: MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER:: ISENDORRECV(MAXPROC)
      INTEGER:: NDF,I,IFUNC,IEND,II
      INTEGER:: IPROC,INDEX,IPF,J
      INTEGER:: ISTART,ISET,MPIINFSIZE,NCONSTR
      INTEGER:: N,NBOU2,NBN,NBOU1
      INTEGER:: NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
      INTEGER:: NCONST, NPLANE, NODE1,NN,NESET
      INTEGER:: NDF1,NEL,NEQ,NELX,NODE,NDFF,NGAUSS
      INTEGER:: NNOD, NODES,NX,NSET
      INTEGER:: NPROCS,NPRESS
      REAL(8):: PRESS_VAL

      COMMON/CBELSET/NESET,NELSET,NELLIST
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDFF,NGAUSS
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC ! NEQ_uBC: total # of DOFs with displ. BCs
      COMMON/CBNSET/NSET,NNSET,NLIST
      COMMON/PRESS/PRESSB,IBELEM,ILF,IFACE_TET,IFACE_BRICK
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV

      CALL READSTR(LR)

      !    imposed displacements
      !----------------------------
      READ(LR,*) NBOU1  ! # of entries to read

      NDF1=NDF+1
      NBOUD1=0
      DO N=1,NBOU1
         READ(LR,*) NBN,ISET,IFUNC,ISTART,IEND
         DO I=1,NDF
            IK(I)=0
            DK(I)=0.D0
         ENDDO

         IF(IFUNC.EQ.0)THEN            ! displacements given as inputs
            READ(LR,*)(DK(I),I=ISTART,IEND)
         ENDIF
         DO I=ISTART,IEND
            IK(I)=1
         ENDDO
         NNOD=1                        ! NBN is a node
         IF(ISET.EQ.1)NNOD=NNSET(NBN)  ! NBN is a node set

         DO I=1,NNOD
            NBOUD1=NBOUD1+1
            IF(NBOUD1.GT.MBOUND1)THEN  ! for IDF()
               IF(IPROC.EQ.0)THEN
                  WRITE(*,*)'INSUFFICIENT MEMORY-MBOUND1'
               ENDIF
               STOP
            ENDIF
            IF(NDF1*NBOUD1.GT.MBOUND1)THEN  ! for NBC1()
               IF(IPROC.EQ.0)THEN
                  WRITE(*,*)'INSUFFICIENT MEMORY-MBOUND1'
               ENDIF
               STOP
            ENDIF
            NODE1=NBN                        ! if a node
            IF(ISET.EQ.1)NODE1=NLIST(I,NBN)  ! if a node set
            IDF(NBOUD1) = ABS(IFUNC)         ! deniz - detect direction of the code-controlled displacement
            DIR_U(NBOUD1) = SIGN(1,IFUNC)    ! deniz - sign of IFUNC, determines the direction. -( ) to invert.
            DO II=1,NDF1
               ! NBC1 : in blocks of 4 integers, indicates whether DOFs are fixed.
               !        stores: nodeID, DOF1_fixed?, DOF2_fixed?, DOF3_fixed?, nodeID, DOF1_fixed?, DOF2_fixed?, ...
               ! VBC1 : in blocks of 3 reals,
               !        stores: u_DOF1, u_DOF2, u_DOF3, u_DOF1, u_DOF2, u_DOF3, ...
               IF(II.EQ.1)THEN
                  NBC1(NDF1*(NBOUD1-1)+II)=NODE1
               ELSE
                  NBC1(NDF1*(NBOUD1-1)+II)=IK(II-1)
                  IF(IFUNC.EQ.0)VBC1(NDF*(NBOUD1-1)+II-1)=DK(II-1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DO N=1,NEQ
         IBC(N)=0
      ENDDO
      NEQ_uBC = 0
      DO N=1,NBOUD1
         NODES=NBC1(NDF1*(N-1)+1)
         DO I=1,NDF
            INDEX=NBC1(NDF1*(N-1)+1+I)
            IF(INDEX.EQ.1)THEN   ! if this DOF is fixed
               NN=NDF*(NODES-1)+I
               ! counter for number of constrained displacements:
               if(IBC(NN).EQ.0) NEQ_uBC = NEQ_uBC + 1
               IBC(NN)=1         ! mark the EQN number as fixed
               ! IBC(EQN) : indicates whether DOF corresp to the EQN number is fixed (i.e. a displ. BC is defined).
            ENDIF
         ENDDO
      ENDDO
!
!     READ THE 2-ND TYPE BOUNDARY CONDITION : POINT FORCES
!
      CALL READSTR(LR)
      READ(LR,*) NBOU2
      NBOUD2=0
      DO N=1,NBOU2
         READ(LR,*) NBN,ISET,IFUNC,ISTART,IEND
         DO I=1,NDF
            IK(I)=0
            DK(I)=0.D0
         ENDDO

         IF(IFUNC.EQ.0)THEN
            READ(LR,*)(DK(I),I=ISTART,IEND)
         ENDIF
         DO I=ISTART,IEND
            IK(I)=1
         ENDDO
         NNOD=1
         IF(ISET.EQ.1)NNOD=NNSET(NBN)

         DO I=1,NNOD
            NBOUD2=NBOUD2+1
            IF(NBOUD2.GT.MBOUND2)THEN
               IF(IPROC.EQ.0)THEN
                  WRITE(131,*)'INSUFFICIENT MEMORY-MBOUND2'
               ENDIF
               STOP
            ENDIF
            NODE1=NBN
            IF(ISET.EQ.1) THEN
            NODE1=NLIST(I,NBN)
            ENDIF
            IDF(NBOUD1) = ABS(IFUNC)        ! deniz - detect direction of the code-controlled displacement
            DIR_U(NBOUD1) = SIGN(1,IFUNC)   ! deniz - sign of IFUNC, determines the direction. -( ) to invert.
            DO II=1,NDF1
               IF(II.EQ.1)THEN
!                  NBC2(NDF1*(NBOUD1-1)+II)=NODE1                        <-----------------------------------
                   NBC2(NDF1*(NBOUD1-1)+II)=0
               ELSE
                  IF(NDF1*(NBOUD2-1)+II.GT.MBOUND2)THEN
                     IF(IPROC.EQ.0)THEN
                        WRITE(131,*)'INSUFFICIENT MEMEORY-MBOUND1'
                     ENDIF
                     STOP
                  ENDIF

! deniz - fix this - creep BC's are not set!

!                  NBC2(NDF1*(NBOUD2-1)+II)=IK(II-1)
                  NBC2(NDF1*(NBOUD2-1)+II)=0
!                  IF(IFUNC.EQ.0)VBC2(NDF*(NBOUD2-1)+II-1)=DK(II-1)
                  IF(IFUNC.EQ.0)VBC2(NDF*(NBOUD2-1)+II-1)=0.0D0
                  
! deniz - fix this - creep BC's are not set!

               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DO N=1,NEQ
         IBCT(N)=0
      ENDDO
      DO N=1,NBOUD2
         NODES=NBC2(NDF1*(N-1)+1)
         DO I=1,NDF
            INDEX=NBC2(NDF1*(N-1)+1+I)
            IF(INDEX.EQ.1)THEN
               NN=NDF*(NODES-1)+I
               IBCT(NN)=1
            ENDIF
         ENDDO
      ENDDO
!
!      READ THE THIRD TYPE BOUNDARY CONDITION (PRESSURE)
!
      CALL READSTR(LR)
      READ(LR,*) NBOUD3

!-----------------------------------
      ! face 1 (CODE) ==== face 4 (deniz)
      IFACE_TET(1)=1
      IFACE_TET(2)=2
      IFACE_TET(3)=3

      ! face 2 (CODE) ==== face 3 (deniz)
      IFACE_TET(4)=1
      IFACE_TET(5)=2
      IFACE_TET(6)=4

      ! face 3 (CODE) ==== face 1 (deniz)
      IFACE_TET(7)=2
      IFACE_TET(8)=3
      IFACE_TET(9)=4

      ! face 4 (CODE) ==== face 2 (deniz)
      IFACE_TET(10)=1
      IFACE_TET(11)=3
      IFACE_TET(12)=4
!------------------------------------

      IFACE_BRICK(1)=1
      IFACE_BRICK(2)=2
      IFACE_BRICK(3)=3
      IFACE_BRICK(4)=4

      IFACE_BRICK(5)=5
      IFACE_BRICK(6)=6
      IFACE_BRICK(7)=7
      IFACE_BRICK(8)=8

      IFACE_BRICK(9)=2
      IFACE_BRICK(10)=3
      IFACE_BRICK(11)=7
      IFACE_BRICK(12)=6

      IFACE_BRICK(13)=1
      IFACE_BRICK(14)=4
      IFACE_BRICK(15)=8
      IFACE_BRICK(16)=5

      IFACE_BRICK(17)=4
      IFACE_BRICK(18)=3
      IFACE_BRICK(19)=7
      IFACE_BRICK(20)=8

      IFACE_BRICK(21)=1
      IFACE_BRICK(22)=2
      IFACE_BRICK(23)=6
      IFACE_BRICK(24)=5


!------------------------------------
      DO N=1,NBOUD3
         READ(LR,*)NPRESS,IPF,(IPRESS_SET(I),I=1,NPRESS), & 
              (IPRESS_SF(I),I=1,NPRESS)
         if (NPRESS.GT.MAXPRESS.AND.IPROC.EQ.0)THEN
            WRITE(131,*)'Increase MAXNSET to ', NPRESS
            WRITE(*,*)'Increase MAXNSET to ', NPRESS
            STOP
         ENDIF


         IF(IPF.EQ.0)THEN
            READ(LR,*)PRESS_VAL
         ENDIF

         DO I=1,NPRESS
            ISET=IPRESS_SET(I)
            DO J=1,NELSET(ISET)
               NEL=NELLIST(J,ISET)
               IBELEM(NEL)=IPRESS_SF(I)
               IF(IPF.EQ.0)THEN
                  ILF(NEL)=0
                  PRESSB(NEL)=PRESS_VAL   ! WARNING! this option is ignored in FEM UEL assembly
                                          ! WARNING! this stores one single pressure value per element-
                                          ! will only support uni-axial creep !
                  WRITE(*,*) 'EXPLICIT PRESSURE B.C. NOT SUPPORTED.'
                  WRITE(*,*) 'SET flag=1, use dload()'
                  STOP
               ELSE
                  ILF(NEL)=1  ! flag=1: loading applied by dload()
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      CALL READSTR(LR)
!     READ CONSTRAINTS
      READ(LR,*)NCONST
      NCONSTR=0

      IF(NCONST.GT.0)THEN
         WRITE(*,*)'CONSTRAINTS NOT IN CODE'
         STOP
      ENDIF


      DO I=1,MBOUND2
      NBC2(I)=0
      VBC2(I)=0.0D0
      ENDDO

      DO I=1,MBOUND3
      NBC3(I)=0
      VBC3(I)=0.0D0
      ENDDO

!-------------------------------------------------------------

      RETURN
      END SUBROUTINE

!************************************************************
      SUBROUTINE TEMPINC(temp_loading_type, temp_stress_free, & 
                      temp_cyclic_amplitude, temp_cyclic_period, & 
                      temp_history_values, temp_history_times, & 
                      temp_history_N, TIME, TEMP, TEMP_0)
      implicit none
      INCLUDE 'PARDIS.H'
      real(8) :: PI
      parameter(PI=3.14159265359D0)
      
      ! INPUTS
      integer, intent(in) :: temp_loading_type ! 1: temperature time history given
                                               ! 2: cyclic thermal loading around initial temp
                                               ! 3: cyclic thermal loading ABOVE initial temp
      real(8), intent(in) :: temp_stress_free  ! initial/stress-free temperature
      real(8), intent(in) :: temp_cyclic_amplitude   ! in Kelvins
      real(8), intent(in) :: temp_cyclic_period      ! in seconds
      real(8), intent(in) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8), intent(in) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
      integer, intent(in) :: temp_history_N                          ! for integration time steps
      real(8), intent(in) :: TIME(2)
      
      ! OUTPUTS
      real(8), intent(out) :: TEMP(2), TEMP_0
      
      !locals
      integer :: iTemp, t_or_tau
      real(8) :: w1,w2
      logical :: found
      
      !statements
      
      TEMP_0 = temp_stress_free
      TEMP(:)=TEMP_0 ! init. to default

      if(temp_loading_type.EQ.2) then  ! thermal cycling, +/- around initial temp
         TEMP(1) = TEMP_0+ & 
             temp_cyclic_amplitude*SIN(2*PI*TIME(1)/temp_cyclic_period)  ! temp. at time t
         TEMP(2) = TEMP_0+ & 
             temp_cyclic_amplitude*SIN(2*PI*TIME(2)/temp_cyclic_period)  ! temp. at time tau
      elseif (temp_loading_type.EQ.3) then  ! thermal cycling, above initial temp
         TEMP(1) = TEMP_0+ & 
           temp_cyclic_amplitude*SIN(PI*TIME(1)/temp_cyclic_period)**2.0  ! temp. at time t
         TEMP(2) = TEMP_0+ & 
           temp_cyclic_amplitude*SIN(PI*TIME(2)/temp_cyclic_period)**2.0  ! temp. at time tau
      elseif(temp_loading_type.EQ.1) then
         do t_or_tau=1,2
            ! interpolate temperature at integration time steps
            ! from the provided temp-time points
            
            found=.FALSE.
            do iTemp=1,temp_history_N
               if(TIME(t_or_tau).LT.temp_history_times(iTemp))then
                  found=.TRUE.
                  exit
               endif
            enddo
            if(found)then
               if(iTemp.EQ.1)then
                  TEMP(t_or_tau)=temp_history_values(1)
               else
                  w1 = (temp_history_times(iTemp)-TIME(t_or_tau)) / & 
               (temp_history_times(iTemp)-temp_history_times(iTemp-1))
                  w2 = (TIME(t_or_tau)-temp_history_times(iTemp-1)) / & 
               (temp_history_times(iTemp)-temp_history_times(iTemp-1))
                  TEMP(t_or_tau)=w1*temp_history_values(iTemp-1) & 
                               + w2*temp_history_values(iTemp)
               endif 
            else
               ! out of provided data range. assign the last provided value
               TEMP(t_or_tau)=temp_history_values(temp_history_N)
            endif
         enddo
      endif

      
!      REAL(8) :: SVARS1(MAXSV),TIME1(2)
!      INTEGER::NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,MCRD1
!      COMMON/ABQ/TIME1,NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,SVARS1,MCRD1

      END SUBROUTINE

      SUBROUTINE TRAINC(NBC2,VBC2,VBCS2,DELT,NDF,TRAC)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV
      DIMENSION VBC2(MBOUND2),VBCS2(MBOUND2),NBC2(MBOUND2), & 
           TRAC(MDOFX)

      NDF1=NDF+1
      DO N=1,NBOUD2
         NODE=NBC2(NDF1*(N-1)+1)
         DO I=1,NDF
            INDEX=NBC2(NDF1*(N-1)+I)
            IF(INDEX.EQ.1)THEN
               NN=NDF*(NODE-1)+I
!               NN=IEQNO(NN)
               N1=NDF*(N-1)+I
               TRAC(NN)=TRAC(NN)+VBCS2(N1)*DELT
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE

!*******************************************************************

      ! deniz - modifications:
      ! NDF, NGAUSS removed from arguments
      ! VT, instead of VTAU is being passed to LOADVC
      ! F and RESID are now properly initialized to zero
      ! added output vector, F_ext(:) -- stores the constraint/reaction forces.
      !  (note that F (oob forces) = F_free (directly out of UMAT for given u for a free structure) - F_ext
      !         and F_ext = F_rxn + F_appl, where F_rxn is assigned the forces of displ. constrained DOFs from F_free.
      SUBROUTINE LOADVC(F,F_ext, &
      stressAverageElements,ddsddeAverageElements, &
      FtauAverageElements,dpAverageElements, &
      volumeTotalElements,volume0TotalElements, &
      VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
      PROMAT2D,DELTME,NBC1,IBC,KINC,N_UPDATE,IN_ST,IROWPTR, & 
      IBINDX,VAL,RESID,NNZ_LOC,VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,   &
      N_UPDATE_rxn,ISIZE_rxn,IFLAG,PNEWDT)
               ! IN
         ! 	VTAU:	   trial displ. field for time tau,  v_tau
         !	VT:	   converged displ. field for time t, v_t
         !	XUND: 	nodal positions in underformed config.
         !  IN_ST:   first EQN number - 1 within the range assigned to this proc.
         !  IFLAG:   0 -> re-calculate stiffness matrix.
         ! OUT
         !	VAL:	   stiffness matrix corresp to this processor
         ! 	RESID:	portion of the FEM force residual (RHS), corresp to this processor
         ! 	F:	      the whole (reduced) FEM force residual vector. i.e. out of balance FEM forces and RHS of the NR iterations K*ddU=F
         ! 	F_ext:	the whole (reduced) FEM external force vector. includes both support reaction forces and applied external forces.
         !  (note that F (oob forces) = F_free (directly out of UMAT for given u for a free structure) - F_ext
         !         and F_ext = F_rxn + F_appl, where F_rxn is assigned the forces of displ. constrained DOFs from F_free.
         !	PNEWDT:	modifier for the displacement increment factor
         !			0.5 : ignore calculation, cut, repeat
         !			>1	: use larger steps in following calculations
         ! now every PROC has their portion of SK and F
         ! move on to SuperLU.

      use task_Timing   ! process timer
      use options, only: DEBUG_timer,DEBUG_RESID, &
                  StiffnessMatrixFormulation,SM_delPdelF,SM_delSdelLnU
      use fBarPatch
      
      IMPLICIT NONE
      INCLUDE 'PARDIS.H'
      INCLUDE 'mpif.h'

      INTEGER:: N_UPDATE
      ! deniz modify
      
      real(8), intent(inout)::VAL_rxn(ISIZE_rxn)
      integer, intent(in)   ::IBINDX_rxn(ISIZE_rxn), & 
                              ICOLPTR_rxn(N_UPDATE_rxn)
      integer, intent(in)   ::ISIZE_rxn,N_UPDATE_rxn
      
      ! macropscopic stress calculation
      real(8), intent(out) :: stressAverageElements(6),ddsddeAverageElements(6,6)
      real(8), intent(out) :: FtauAverageElements(3,3),dpAverageElements(6)
      real(8), intent(out) :: volumeTotalElements,volume0TotalElements
      real(8) :: stressElement(6),ddsddeElement(6,6)
      real(8) :: FtauElement(3,3),dpElement(6)
      real(8) :: volumeElement,volumeElement0
      
      ! deniz - FBar      
      integer :: patchElemIdx,patchElemID,patchID,elemIdx
      integer, allocatable :: mapDofLocalToGlobal_ePatch(:,:)
      real(8), allocatable :: EXTRNVAL_Patch(:,:,:)
      real(8), allocatable :: EXTRNMPI_Patch(:)
      real(8), allocatable :: RECVVAL_Patch(:)
      
      real(8)::SKE(MDOFE,MDOFE),FE(MDOFE),FE_appl(MDOFE)
      real(8)::SKE_Patch(MDOFE,MDOFE,MaxElemInPatch)
      REAL(8)::XEUN(MAXDIM,MAXELN),XUND(MAXCRD)
      REAL(8)::VTAU(NEQ),VT(NEQ)
      REAL(8)::PROMAT2D(MAXPROPS,MAXGRP),VTAU_E(MAXELN*NDOFELX)
      REAL(8)::PROPS1(MAXPROPS),TIME1(2),SVARS1(MAXSV),TIME(2)
      REAL(8)::VT_E(MAXELN*NDOFELX)
      REAL(8)::EXTRNVAL(MAXDIM*8,MDOFPROC),EXTRNRHS(MDOFPROC,MAXPROC)
      REAL(8)::RHSARR(MDOFPROC),EXTRNMPI(MDOFPROC*MAXDIM*MAXELN)
      REAL(8)::RECVVAL(MDOFPROC*8*MAXDIM),RECVRHS(MDOFPROC)
      REAL(8)::DELTME,PNEWDT,PNEWDT_LOC,PNEWDT_RECV

      INTEGER::IMID(MAXDIM*MAXEL),IJK(MNELX),NBC1(MBOUND1)
      INTEGER::NPROPSG(MAXGRP),IBC(MDOFX)
      INTEGER::NUMROWSPROC(MAXPROC),IROWNO(MDOFPROC,MAXPROC)
      INTEGER::IELNO(MDOFPROC,MAXPROC),INDXCOLM(MDOFPROC,MAXPROC)
      INTEGER::IROWINDX(MDOFPROC,MAXPROC),mapDofLocalToGlobal_e(MDOFX*MAXELN)
      INTEGER:: IROWARR(MDOFPROC),IELARR(MDOFPROC),IRINDX(MDOFPROC)
      INTEGER::NDF,IN_ST,NNZ_LOC
      INTEGER::NGAUSS,IROW,IIB,idEqnGlobal
      INTEGER::IE,ICOLM,IEL,IERROR,IEND
      INTEGER::IEXTRN,IIA,II,III

      INTEGER:: IRECVROW(MDOFPROC),IRECVELNO(MDOFPROC)
      INTEGER:: IRECVRINDX(MDOFPROC),ISTATUS(MPI_STATUS_SIZE)
      INTEGER:: NNODE
      INTEGER:: NELX,NEQ,NPLANE,NODE
      INTEGER:: NX,NSVARS1,NDOFEL1,MDLOAD1
      INTEGER:: NPREDF1,MCRD1,MPIINF(MPIINFSIZEMAX)
      INTEGER:: IPROCINFO(MAXPROC), ISENDORRECV(MAXPROC)
      INTEGER:: NELST,NELEND,NELNOX
      INTEGER:: IBINDX(NNZ_LOC)
      INTEGER:: IROWPTR(N_UPDATE+1)
      INTEGER:: INDX1,I,J,IRNO,ISEND
      INTEGER:: IS,JELEM,ITEMP,JKIA,JJ
      INTEGER:: JINDX,JKIA1,K,LL,KK
      INTEGER:: MCRD,MDLOAD,MLVARX,MPIPROCNO
      INTEGER:: NDOFEL,NNO,NN
      ! deniz modify
      INTEGER:: matID
      INTEGER:: NODEINPROC,NPREDF,IPROC
      INTEGER:: MPIINFSIZE, NEL, NROWS
      INTEGER:: NPROCS,NRECVROWS,NPROPS,NSVARS
      
      ! deniz - modification
      real(8) :: F_appl(MDOFX)    ! applied forces for all processors. will be reduced
      real(8) :: F_rxn(MDOFX)     ! constraint/reaction forces for all processors. will be reduced
      real(8), intent(out):: F(MDOFX),F_ext(MDOFX)
      real(8), intent(out):: RESID(N_UPDATE),VAL(NNZ_LOC)
      integer, intent(in) :: IFLAG,KINC
      
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/ABQ/TIME1,NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,SVARS1,MCRD1
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV
      COMMON/EL_PAR/NELST,NELEND,NELNOX
      
      ! error triggers:
      logical :: listLocalUmatErrorsOccurred(7)
      logical :: trigger_umatError
      logical :: trigger_BroyDidNotConv
      logical :: trigger_cutBack
      logical :: trigger_reformK
      COMMON/triggers/trigger_umatError,listLocalUmatErrorsOccurred, &
                      trigger_BroyDidNotConv,trigger_cutBack,trigger_reformK
      
      ! deniz - error flag
      integer :: umat_Err
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT  
      
      real(8) :: task_time
      
      character(len=32) :: FILNAME
      
      ! initialize the list of local umat errors
      listLocalUmatErrorsOccurred(:) = .false.
      
      
      ! initialize macropscopic stress calculation variables
      stressAverageElements = 0.d0
      dpAverageElements = 0.d0
      FtauAverageElements = 0.d0
      volumeTotalElements = 0.d0
      volume0TotalElements = 0.d0
      if (IFLAG == 0) then ! reconstruct stiffness matrix
         ddsddeAverageElements = 0.d0
      endif

      RESID(1:N_UPDATE)=0.D0
      F_appl(1:NEQ)=0.D0
      F_rxn(1:NEQ)=0.D0
      IF(IFLAG.EQ.0)THEN   ! will reconstruct stiffness matrix. initialize it.
         VAL = 0.D0
      endif
      
      ! initialize error trigger. UMAT will set this to TRUE is an error occurs.
      trigger_umatError = .FALSE.

      TIME(1)=TIME1(1)
      TIME(2)=TIME1(2)
      NSVARS=NSVARS1
      NDOFEL=NDOFEL1
      MDLOAD=MDLOAD1
      NPREDF=NPREDF1
      MCRD=MCRD1
      MLVARX=NDOFEL

      allocate(mapDofLocalToGlobal_ePatch(MAXDIM*MAXELN,MaxElemInPatch))
      allocate(EXTRNVAL_Patch(MAXDIM*8,MaxElemInPatch,MDOFPROC))
      allocate(EXTRNMPI_Patch(MDOFPROC*MAXDIM*MAXELN*MaxElemInPatch))
      allocate(RECVVAL_Patch(MDOFPROC*8*MAXDIM*MaxElemInPatch))
      
!****************************LINEAR PARTITIONING OF ELEMENTS***************************

      DO I=1,NPROCS-1
         NUMROWSPROC(I)=0
      ENDDO

      IEXTRN=1
!***************************************************************************************
      PNEWDT_LOC=10.D0
      
      ! start ELST timer
      if(DEBUG_timer) CALL task_Start(task_loadVC_elst,MPI_WTIME())

      DO NEL=NELST,NELEND

         NNODE=IMID(NEL)			! # of nodes in this element
         matID=IMID(NELX+NEL)	! mat. ID
         NPROPS=NPROPSG(matID)

		! get material parameters for the element
         DO II=1,NPROPS
            PROPS1(II)=PROMAT2D(II,matID)
         ENDDO

         IIB=0
         DO IIA=1,NODE
            JKIA=IJK(NODE*(NEL-1)+IIA)
            JKIA1=MCRD*(JKIA-1)
            DO III=1,MCRD
               XEUN(III,IIA)=XUND(JKIA1+III)	! xeun(i,n) = coordinate i of the initial position of the element node n.
            ENDDO
            DO II=1,NDF
               NN=NDF*(JKIA-1)+II
               IIB=IIB+1
               VTAU_E(IIB)=VTAU(NN)
               VT_E(IIB)=VT(NN)
            ENDDO
         ENDDO

         JELEM=NEL
         NDOFEL=NODE*MCRD
         
         curElemUMAT = JELEM
         
         
         !VT_umat = t_umat*VT_E + (1.d0-t_umat)*VTAU_E
         !VTAU_umat = tau_umat*VT_E + (1.d0-tau_umat)*VTAU_E
         ! GOTTO CHANGE THE ARGUMENTS OF ELST1 -- in SVARS_t, out: svars_tau, instead of common blocks..
         ! follow your notes

         CALL ELST1(JELEM,NODE,XEUN,SKE,SKE_Patch,FE,FE_appl,TIME,DELTME, &
         stressElement,ddsddeElement,FtauElement,dpElement, &
         volumeElement,volumeElement0, &
         PROPS1,NPROPS,NSVARS,NDOFEL,MDLOAD,MLVARX,NPREDF,MCRD,VTAU_E,VT_E, & 
         KINC,IFLAG,PNEWDT,umat_Err)
         !OUT
         !  SKE1: matrix, corresp to a single elem.
         !  FE:   residual correspt to a single element
         !  PNEWDT
         !        step modifier
         
         ! if a UMAT error has occurred, PNEWDT is set to 0.5.
         ! in that case, proceed until the end of the communication 
         ! to retain processor synchronization, then set trigger_umatError
         ! and leave the subroutine
         
         ! record any local umat errors have occurred
         if(umat_Err /= 0) then
            listLocalUmatErrorsOccurred(umat_Err) = .true.
         endif

         IF(PNEWDT.LT.PNEWDT_LOC)PNEWDT_LOC=PNEWDT
         
         ! calculation of macroscopic (RVE-average) quantities.
         ! they will be mpi-reduced and normalized at the end of the subroutine 
         volumeTotalElements = volumeTotalElements + volumeElement
         volume0TotalElements = volume0TotalElements + volumeElement0
         stressAverageElements = stressAverageElements + stressElement*volumeElement
         dpAverageElements = dpAverageElements + dpElement*volumeElement
         FtauAverageElements = FtauAverageElements + FtauElement*volumeElement0
         if (IFLAG == 0) then ! stiffness matrix recalc
            ddsddeAverageElements = ddsddeAverageElements + ddsddeElement*volumeElement
         endif
         
         ! CREATING TRANSPOSE MATRIX
         ! only in the old (dSdlnU) formulation
         IF(IFLAG.EQ.0 .and. StiffnessMatrixFormulation==SM_delSdelLnU) THEN
            SKE(1:NDF*NODE,1:NDF*NODE) = TRANSPOSE(SKE(1:NDF*NODE,1:NDF*NODE))
         ENDIF

         ! assemble the local element applied forces to the global appl. force array (FULL array, reduced at the end of ELST1)
         DO J=1,NODE
            NNO=IJK((NEL-1)*NODE + J)
            F_appl((NNO-1)*NDF+1:(NNO-1)*NDF+NDF) =  & 
                  F_appl((NNO-1)*NDF+1:(NNO-1)*NDF+NDF) & 
                + FE_appl((J-1)*NDF+1:(J-1)*NDF+NDF)
         enddo

!***************************************************************************************
         DO II=1,NODE
            NNO=IJK((NEL-1)*NODE+II)
            DO JJ=1,NDF
               mapDofLocalToGlobal_e((II-1)*NDF+JJ)=(NNO-1)*NDF+JJ   ! store the EQN numbers for this element
            ENDDO
         ENDDO

         ! FBar - save the local->global EQN maps for the elements in the patch
         patchID = patchIDfromElemID(NEL)
         do patchElemIdx = 1,NElementsInPatch(patchID)
            patchElemID = ElementsInPatch(patchID,patchElemIdx)
            DO II=1,NODE
               NNO=IJK((patchElemID-1)*NODE+II)
               DO JJ=1,NDF
                  mapDofLocalToGlobal_ePatch((II-1)*NDF+JJ,patchElemIdx)=(NNO-1)*NDF+JJ   ! store the EQN numbers for this element
               ENDDO
            ENDDO
         enddo
         
         DO J=1,NODE
            NNO=IJK((NEL-1)*NODE + J)
            CALL FINDPROC(NNO,NPROCS,NODEINPROC)
            ! nodes are also assigned to processors (again, linear fashion, according to node IDs)

            DO K=1,NDF

               IROW=(NNO-1)*NDF+K   ! EQN number of this DOF of this NODE of this ELEM

               IF(NODEINPROC.EQ.IPROC)THEN   ! no need to send

                  INDX1=IROW-IN_ST
                  if(INDX1.GT.N_UPDATE+1)then
                     write(*,*) INDX1, N_UPDATE, IPROC, IROW, IN_ST
                     continue
                  endif
                  IS=IROWPTR(INDX1)
                  IE=IROWPTR(INDX1+1)-1
                  ! VAL(:), IROWPTR(:), IBINDX(:) -- Stiffness Matrix stored in Sparse Row Format. 
                  ! See SuperLU Users Guide.

                  IF(IFLAG.EQ.0)THEN

                     DO KK=1,NODE*NDF  ! looping over the row (node J DOF K) of the element stiffness matrix
                        idEqnGlobal=mapDofLocalToGlobal_e(KK)
                        DO LL=IS,IE ! nonzero coefficients at this row
                           IF(IBINDX(LL).EQ.idEqnGlobal)THEN
                              VAL(LL)=VAL(LL)+SKE((J-1)*NDF+K,KK)
                              GOTO 110
                           ENDIF
                        ENDDO
 110                 ENDDO
 
                     ! now loop over the same row, but the parts coming from each patch element
                     do patchElemIdx = 1,NElementsInPatch(patchID)
                        patchElemID = ElementsInPatch(patchID,patchElemIdx)
                        if (patchElemID == NEL) cycle ! this is 'me' - skip.
                        DO KK=1,NODE*NDF  ! looping over the row (node J DOF K) of the K of the patch element
                           idEqnGlobal=mapDofLocalToGlobal_ePatch(KK,patchElemIdx)
                           DO LL=IS,IE ! nonzero coefficients at this row
                              IF(IBINDX(LL).EQ.idEqnGlobal)THEN
                                 VAL(LL)=VAL(LL)+SKE_Patch((J-1)*NDF+K,KK,patchElemIdx)
                                 GOTO 210
                              ENDIF
                           ENDDO
 210                    ENDDO
                     enddo


                  ENDIF

               ! assemble the local element residual to the global residual array for this PROC
               RESID(INDX1)= RESID(INDX1) + FE((J-1)*NDF + K)              

               ELSE  ! DOF in different processor, need to send K coefficients

                  IF(NODEINPROC.LT.IPROC)THEN
                     INDX1 = NODEINPROC+1
                  ELSE
                     INDX1 = NODEINPROC
                  ENDIF

                  ! NODEINPROC: receiver processor ID
                  ! INDX1: receiver proc index for the data-buffer
                  ! JINDX: counter for number of rows/EQNs being sent to INDX1
                  JINDX=NUMROWSPROC(INDX1)+1 ! NUMROWSPROC how many rows to send to the proc
                  NUMROWSPROC(INDX1)=JINDX   ! increase it by one.
                  IROWNO(JINDX,INDX1)=IROW   ! global EQN/ROW ID
                  IELNO(JINDX,INDX1)=NEL     ! check subroutine update local array. it receives the package and updates necessary local vars.
                  INDXCOLM(JINDX,INDX1)=IEXTRN  ! the start index of the data package (in number of total rows sent)
                  EXTRNRHS(JINDX,INDX1)=FE((J-1)*NDF+K)  ! RHS that I'm sending to the proc. INDX1: receiver procID
                  IROWINDX(JINDX,INDX1)=(J-1)*NDF+K   ! local EQN/DOF ID
                  
                  IF(IFLAG.EQ.0)THEN ! re-calc stiffness matrix

                     DO KK=1,NODE*NDF
                        EXTRNVAL(KK,IEXTRN)=SKE((J-1)*NDF+K,KK)   !save the row of the element stiffness matrix to send-buffer
                     ENDDO
                     
                     ! loop over the FBar patch elements to save their contributions to the row
                     ! for sending out
                     do patchElemIdx = 1,NElementsInPatch(patchID)
                        patchElemID = ElementsInPatch(patchID,patchElemIdx)
                        
                        EXTRNVAL_Patch(:,patchElemIdx,IEXTRN) = 0.d0
                        if (patchElemID == NEL) cycle ! this is 'me' - skip.
                        
                        DO KK=1,NODE*NDF
                           EXTRNVAL_Patch(KK,patchElemIdx,IEXTRN) = &
                              SKE_Patch((J-1)*NDF+K,KK,patchElemIdx)   !save the row of the element stiffness matrix to send-buffer
                        ENDDO

                     enddo
                     
                     IEXTRN=IEXTRN+1   ! index counter for the rows sent from my processor

                  ENDIF ! re-calc stiffness matrix

                  IF(JINDX.GT.MDOFPROC.OR.IEXTRN.GT.MDOFPROC)THEN
                     WRITE(*,*)'EXTERNAL COMMUNICATION VARIABLE & 
                     &MDOFPROC SHOULD BE INCREASED'
                     WRITE(*,*)JINDX,IEXTRN,MDOFPROC
                     STOP
                  ENDIF

               ENDIF

            ENDDO ! DOFs of the node
         ENDDO    ! nodes of this element

      ENDDO ! elements handled by this processor
      
      ! record ELST timer
      if (DEBUG_timer) CALL task_Record(task_loadVC_elst,MPI_WTIME(),task_time,-1)
      
!********************************************************************

      ! all processes have looped over their elements
      ! Jacobian and residual is complete
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(PNEWDT_LOC,PNEWDT_RECV,1, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_MIN,MPI_COMM_WORLD,IERROR)

      PNEWDT=PNEWDT_RECV
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,listLocalUmatErrorsOccurred,7, &  ! OR-Reduce the boolean flags of processor-local umat errors
         MPI_LOGICAL,MPI_LOR, & 
         MPI_COMM_WORLD,IERROR)

      IF(PNEWDT.LT.1.D0) then       ! an error had occured, VAL and RESID are corrupt. 
         trigger_umatError = .TRUE. ! this will trigger a recalculation of K
         RETURN
      endif

      ! start Communication timer
      if(DEBUG_timer) CALL task_Start(task_loadVC_comm,MPI_WTIME())
      
      DO I=1,NPROCS-1



         MPIPROCNO=IPROCINFO(I)


         ISEND=ISENDORRECV(I)

         IF(IPROC.GT.MPIPROCNO) THEN
            ITEMP=MPIPROCNO+1
         ELSE
            ITEMP=MPIPROCNO
         ENDIF

         IF(ISEND.EQ.1)THEN      ! See PROC_COMB . generates the send/recv pairs. this determines ISEND.

            NROWS=NUMROWSPROC(ITEMP)   ! ITEMP: receiver PROC ID
                                       ! NROWS = # of rows to send to that processor


            CALL MPI_SSEND(NROWS,1,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

            EXTRNMPI = 0.d0
            EXTRNMPI_Patch = 0.d0

            DO J=1,NROWS   ! nRows being sent to this processor
               IRNO=IROWNO(J,ITEMP)
               IROWARR(J)=IRNO
               ICOLM=INDXCOLM(J,ITEMP)
               IEL=IELNO(J,ITEMP)
               IELARR(J)=IEL
               IRINDX(J)=IROWINDX(J,ITEMP)
               RHSARR(J)=EXTRNRHS(J,ITEMP)

               IF(IFLAG.EQ.0)THEN
                  ! 'pack' the row to be sent, into 1D array
                  EXTRNMPI((J-1)*NODE*NDF+1:J*NODE*NDF)= &
                        EXTRNVAL(1:NODE*NDF,ICOLM)
                  ! 'pack' the rows fro patch elements into 1D array
                  patchID = patchIDfromElemID(IEL)
                  do patchElemIdx = 1,NElementsInPatch(patchID)
                     patchElemID = ElementsInPatch(patchID,patchElemIdx)
                     if (patchElemID == IEL) cycle ! this is 'me' - skip.
                     
                     EXTRNMPI_Patch((J-1)*NODE*NDF*MaxElemInPatch+ &
                                    (patchElemIdx-1)*NODE*NDF+1: &
                                    (J-1)*NODE*NDF*MaxElemInPatch+ &
                                    (patchElemIdx)*NODE*NDF) = &
                           EXTRNVAL_Patch(1:NODE*NDF,patchElemIdx,ICOLM)
                  enddo
               ENDIF
            ENDDO

            IF(NROWS.GT.0)THEN

               ! SEND OUT THE PACKAGES

               CALL MPI_SSEND(IROWARR,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

               CALL MPI_SSEND(IELARR,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

               CALL MPI_SSEND(IRINDX,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)


               CALL MPI_SSEND(RHSARR,NROWS,MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,IERROR)

               IF(IFLAG.EQ.0)THEN
            CALL MPI_SSEND(EXTRNMPI,NROWS*NODE*NDF,MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,IERROR)
                  ! send FBar coefficients
            CALL MPI_SSEND(EXTRNMPI_Patch,NROWS*NODE*NDF*MaxElemInPatch, &
                         MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,IERROR)
               ENDIF
            ENDIF

            CALL MPI_RECV(NRECVROWS,1,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,ISTATUS,IERROR)

            IF(NRECVROWS.GT.0)THEN

!               WRITE(*,*)'RECV1',IPROC,MPIPROCNO

               CALL MPI_RECV(IRECVROW,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                         100,MPI_COMM_WORLD,ISTATUS,IERROR)

               CALL MPI_RECV(IRECVELNO,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)

               CALL MPI_RECV(IRECVRINDX,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)

         CALL MPI_RECV(RECVRHS,NRECVROWS,MPI_DOUBLE_PRECISION,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)

               IF(IFLAG.EQ.0)THEN
          CALL MPI_RECV(RECVVAL,NRECVROWS*NDF*NODE,MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,ISTATUS,IERROR)
                  ! receive FBar patch coupled element coeffs
         CALL MPI_RECV(RECVVAL_Patch,NRECVROWS*NODE*NDF*MaxElemInPatch, &
                  MPI_DOUBLE_PRECISION, & 
                  MPIPROCNO,100,MPI_COMM_WORLD,ISTATUS,IERROR)
               ENDIF

            CALL UPDATE_LOCAL_ARRAY(IJK,NRECVROWS,RECVVAL,RECVRHS,IRECVROW, & 
         IRECVELNO,IRECVRINDX,IPROC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL, & 
         RESID,NNZ_LOC,IFLAG)
            CALL UPDATE_couplings_FBar(IJK,NRECVROWS,RECVVAL_Patch,RECVRHS, &
                  IRECVROW,IRECVELNO,IRECVRINDX,IPROC,N_UPDATE,IN_ST,IROWPTR, &
                  IBINDX,VAL,NNZ_LOC,IFLAG)





            ENDIF

         ELSE  ! ISEND /= 0

            CALL MPI_RECV(NRECVROWS,1,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,ISTATUS,IERROR)


            IF(NRECVROWS.GT.0)THEN

!               WRITE(*,*)'RECV2',IPROC,MPIPROCNO

               CALL MPI_RECV(IRECVROW,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                         100,MPI_COMM_WORLD,ISTATUS,IERROR)

               CALL MPI_RECV(IRECVELNO,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)

               CALL MPI_RECV(IRECVRINDX,NRECVROWS,MPI_INTEGER,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)


         CALL MPI_RECV(RECVRHS,NRECVROWS,MPI_DOUBLE_PRECISION,MPIPROCNO, & 
                        100,MPI_COMM_WORLD,ISTATUS,IERROR)

               IF(IFLAG.EQ.0)THEN
         CALL MPI_RECV(RECVVAL,NRECVROWS*NODE*NDF,MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,ISTATUS,IERROR)
                  ! receive FBar patch coupled element coeffs
         CALL MPI_RECV(RECVVAL_Patch,NRECVROWS*NODE*NDF*MaxElemInPatch, &
                  MPI_DOUBLE_PRECISION, & 
                  MPIPROCNO,100,MPI_COMM_WORLD,ISTATUS,IERROR)
               ENDIF

               CALL UPDATE_LOCAL_ARRAY(IJK,NRECVROWS,RECVVAL,RECVRHS,IRECVROW, & 
       IRECVELNO,IRECVRINDX,IPROC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL, & 
       RESID,NNZ_LOC,IFLAG)
               CALL UPDATE_couplings_FBar(IJK,NRECVROWS,RECVVAL_Patch, &
                     RECVRHS,IRECVROW,IRECVELNO,IRECVRINDX,IPROC, &
                     N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL,NNZ_LOC,IFLAG)

            ENDIF

            NROWS=NUMROWSPROC(ITEMP)

            CALL MPI_SSEND(NROWS,1,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

            EXTRNMPI = huge(EXTRNMPI)  ! for debugging - none of these huge numbers ..
            EXTRNMPI_Patch = huge(EXTRNMPI_Patch) ! should make up to the final stifffness matrix

            DO J=1,NROWS   ! nRows being sent to this processor
               IRNO=IROWNO(J,ITEMP)    ! save coeffs to this global EQN/ROW ID
               IROWARR(J)=IRNO         ! pack eqn IDs as array
               ICOLM=INDXCOLM(J,ITEMP) ! the index of the rows sent from my processor
               IEL=IELNO(J,ITEMP)
               IELARR(J)=IEL
               RHSARR(J)=EXTRNRHS(J,ITEMP)
               IRINDX(J)=IROWINDX(J,ITEMP)
               IF(IFLAG.EQ.0)THEN
                  ! 'pack' the row to be sent, into 1D array
                  EXTRNMPI((J-1)*NODE*NDF+1:J*NODE*NDF)= &
                        EXTRNVAL(1:NODE*NDF,ICOLM)
                  ! 'pack' the rows fro patch elements into 1D array
                  patchID = patchIDfromElemID(IEL)
                  do patchElemIdx = 1,NElementsInPatch(patchID)
                     patchElemID = ElementsInPatch(patchID,patchElemIdx)
                     if (patchElemID == IEL) cycle ! this is 'me' - skip.
                     
                     EXTRNMPI_Patch((J-1)*NODE*NDF*MaxElemInPatch+ &
                                    (patchElemIdx-1)*NODE*NDF+1: &
                                    (J-1)*NODE*NDF*MaxElemInPatch+ &
                                    (patchElemIdx)*NODE*NDF) = &
                           EXTRNVAL_Patch(1:NODE*NDF,patchElemIdx,ICOLM)
                  enddo
               ENDIF
            ENDDO


            IF(NROWS.GT.0)THEN

!               WRITE(*,*)'SEND2',IPROC,MPIPROCNO

               CALL MPI_SSEND(IROWARR,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

               CALL MPI_SSEND(IELARR,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

               CALL MPI_SSEND(IRINDX,NROWS,MPI_INTEGER,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

         CALL MPI_SSEND(RHSARR,NROWS,MPI_DOUBLE_PRECISION,MPIPROCNO,100, & 
                         MPI_COMM_WORLD,IERROR)

               IF(IFLAG.EQ.0)THEN
            CALL MPI_SSEND(EXTRNMPI,NROWS*NODE*NDF,MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,IERROR)
                  ! send FBar coefficients
            CALL MPI_SSEND(EXTRNMPI_Patch,NROWS*NODE*NDF*MaxElemInPatch, &
                         MPI_DOUBLE_PRECISION, & 
                         MPIPROCNO,100,MPI_COMM_WORLD,IERROR)
               ENDIF

            ENDIF

         ENDIF


      ENDDO
      
      ! record Communication timer
      if (DEBUG_timer) CALL task_Record(task_loadVC_comm,MPI_WTIME(),task_time,-1)
      
      ! extract rxn matrix, only if a re-assembly was done
      if(IFLAG.EQ.0)THEN
         
         CALL extractRxnMatrix(IBC,N_UPDATE,IN_ST, & 
                               IROWPTR,IBINDX,NNZ_LOC,VAL, & 
                               ICOLPTR_rxn,IBINDX_rxn,VAL_rxn, & 
                               N_UPDATE_rxn,ISIZE_rxn)
                               
      endif
      
      ! extract reaction forces from the residual. purpose: F_rxn + F_appl = F_ext
      CALL extractRxnForce(RESID,F_rxn,IBC,N_UPDATE,IN_ST,NEQ)

      if(DEBUG_RESID.AND.IFLAG.EQ.0)THEN
         
         CALL printMatrix_Sparse(VAL,IBINDX,IROWPTR,  &
                                    N_UPDATE,NEQ,NNZ_LOC,   &
         1,'VAL at STEP 1 before REMOVEDISPDOF',542)

      endif

      ! IMPOSE boundary conditions by directly modifying K and f.
      ! also, set residual (reaction) forces at the displacement-constrained DOFs to zero, at each q-Newton iter
      ! deniz - modified. added argument F_ext.
      CALL REMOVE_DISPDOF(IBC,N_UPDATE,IN_ST, & 
                          IROWPTR,IBINDX,VAL, & 
                          RESID,NNZ_LOC,IFLAG)

      ! note: traction/pressure BCs are added to element force vectors in ELST1:
      ! see ELST1 line: IF(ILF(JELEM).EQ.1)THEN CALL DLOAD(JELEM,PR,TIME)
      ! these will balance with internal element forces, as R -> 0.

      ! added by Deniz
!      CALL IMPOSE_PERIODIC_BCs(coupledDOFs,couplingWeights,nCoupledDOFs,
!     &                         IBC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL,
!     &                         RESID,NNZ_LOC,IFLAG)

      ! ******* REDUCE F_RXN ************
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,F_rxn,NEQ, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
     
     ! ******* REDUCE F_appl ***********
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,F_appl,NEQ, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
       
      ! construct F_ext = F_appl + F_rxn
      F_ext(1:NEQ) = F_appl(1:NEQ) + F_rxn(1:NEQ)
     
     ! ******* REDUCE F ***********
      ! make sure properly initialized before ALL_REDUCE
      F(1:NEQ) = 0.D0
      DO I=1,N_UPDATE
         F(I+IN_ST)=RESID(I)
      ENDDO

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,F,NEQ, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
      
      ! calculate RVE-averaged macroscopic stress
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,stressAverageElements,6, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
      ! calculate RVE-averaged plastic deformation gradient (current config)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,dpAverageElements,6, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
      ! calculate RVE-averaged macroscopic Deformation Gradient
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,FtauAverageElements,9, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
      ! calculate RVE-averaged DDSDDE
      if (IFLAG == 0) then ! stiffness matrix recalc
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,ddsddeAverageElements,6*6, & 
          MPI_DOUBLE_PRECISION, & 
          MPI_SUM,MPI_COMM_WORLD,IERROR)
      endif
      ! calculate undeformed and current volumes of the domain
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,volumeTotalElements,1, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,volume0TotalElements,1, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)
       
      stressAverageElements = stressAverageElements / volumeTotalElements  ! normalize the average stress
      dpAverageElements = dpAverageElements / volumeTotalElements          ! normalize the average dp
      FtauAverageElements = FtauAverageElements / volume0TotalElements     ! normalize the average Deformation Gradient
      if (IFLAG == 0) then ! stiffness matrix recalc
         ddsddeAverageElements = ddsddeAverageElements / volumeTotalElements  ! normalize the average ddsdde
      endif

      deallocate(mapDofLocalToGlobal_ePatch)
      deallocate(EXTRNVAL_Patch)
      deallocate(EXTRNMPI_Patch)
      deallocate(RECVVAL_Patch)

      RETURN
      END SUBROUTINE

!************************************************
      SUBROUTINE calc_TOL(R,NEQ,TOL,tolNorm)
      implicit none
      real(8), intent(in) :: R(NEQ)
      integer, intent(in) :: NEQ
      real(8), intent(out):: TOL
      integer, intent(in) :: tolNorm
      
      if (tolNorm.EQ.1) then   ! L_inf (max) norm
         TOL = MAXVAL(abs(R(1:NEQ)))
      elseif(tolNorm.EQ.2) then! L_2 norm
         TOL = DSQRT(DOT_PRODUCT(R(1:NEQ),R(1:NEQ)))
      endif
      END SUBROUTINE

!**************************************************

      SUBROUTINE READNSET()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      INTEGER:: NNSET(MAXNSET),NLIST(MAXNS,MAXNSET)
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)
      INTEGER:: N
      COMMON/CBNSET/NSET,NNSET,NLIST
      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                           IPROCINFO,ISENDORRECV

      CALL READSTR(LR)
      READ(LR,*)NSET

      IF(NSET.GT.MAXNSET)THEN
         IF(IPROC.EQ.0)THEN
            WRITE(131,*)'INSUFFICIENT MEMORY-MAXNSET'
            WRITE(*,*)'INSUFFICIENT MEMORY-MAXNSET'
         ENDIF
         STOP
      ENDIF

      DO I=1,NSET
         READ(LR,*)NSNO,NGEN
         IF(NGEN.EQ.0)THEN
            READ(LR,*)NUMN,(NLIST(J,NSNO),J=1,NUMN)
         ENDIF
         ! error check added - deniz
         IF(NUMN.GT.MAXNS) THEN
            WRITE(*,*)'Increase MAXNS to ', NUMN
            IF(IPROC.EQ.0)THEN
               WRITE(131,*)'Increase MAXNSET to ', NSET
               STOP
            ENDIF
         END IF
         IF(NSET.GT.MAXNSET) THEN
            WRITE(*,*)'Increase MAXNSET to ', NSET
            IF(IPROC.EQ.0)THEN
               WRITE(131,*)'Increase MAXNSET to ', NSET
               STOP
            ENDIF
         END IF
         IF(NGEN.EQ.1)THEN
            READ(LR,*)ISTART,IEND,INC
            NUMN=0
            DO N=ISTART,IEND,INC
               NUMN=NUMN+1
               NLIST(NUMN,NSNO)=N
            ENDDO
         ENDIF
         NNSET(NSNO)=NUMN
      ENDDO
      RETURN
      END SUBROUTINE

!*******************************************************
      SUBROUTINE READELSET()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION NELSET(MAXELSET),NELLIST(MAXELS,MAXELSET)
      COMMON/CBELSET/NESET,NELSET,NELLIST
      CALL READSTR(LR)
      READ(LR,*)NESET

      IF(NESET.GT.MAXELSET)THEN
         WRITE(131,*)'INSUFFICIENT MEMORY-MAXELSET'
         WRITE(*,*)'INSUFFICIENT MEMORY-MAXELSET'
         STOP
      ENDIF

      DO I=1,NESET
         READ(LR,*)NSNO,NGEN
         IF(NGEN.EQ.0)THEN
            READ(LR,*)NUMN,(NELLIST(J,NSNO),J=1,NUMN)
         ENDIF
         if(NUMN.GT.MAXELS) then
            WRITE(131,*)'INSUFFICIENT MEMORY-MAXELS'
            WRITE(*,*)'INSUFFICIENT MEMORY-MAXELS'
            STOP
         endif
         IF(NGEN.EQ.1)THEN
            READ(LR,*)ISTART,IEND,INC
            NUMN=0
            DO N=ISTART,IEND,INC
               NUMN=NUMN+1
               NELLIST(NUMN,NSNO)=N
            ENDDO
         ENDIF
         NELSET(NSNO)=NUMN
      ENDDO
      RETURN
      END SUBROUTINE

      ! added by deniz
      ! --------------- readTexture --------------------------
      ! read in:
      !     Euler angles,
      !     phase ID
      !     initial hardness params of elements
      ! and save into state vars
      ! assumes these are done:
      !     element partitioning,
      !     INITABQ()
      !     READMT()
      !--------------------------------------------------------
      SUBROUTINE readTexture_OLD(SVARS1,NELX,NELST,NELNOX,NSVARS,NGAUSS)
      implicit none
      
      real(8), intent(inout) :: SVARS1(NSVARS*NELNOX*NGAUSS)
      integer, intent(in) :: NELX,NELST,NELNOX,NSVARS,NGAUSS
      
      integer :: iElem,iEL_proc,iPhase,NELEND
      real(8) :: euler(3)
      
      NELEND = NELST + NELNOX - 1
      
      OPEN(UNIT=201,FILE='texture.dat')

      DO iElem=1,NELX
         READ(201,*) euler(1),euler(2),euler(3),iPhase
         if (iElem.GE.NELST.AND.iElem.LE.NELEND) then
            iEL_proc = iElem - NELST + 1
            ! save EULER
            SVARS1((iEL_proc-1)*NSVARS*NGAUSS + 301)=euler(1)
            SVARS1((iEL_proc-1)*NSVARS*NGAUSS + 302)=euler(2)
            SVARS1((iEL_proc-1)*NSVARS*NGAUSS + 303)=euler(3)
            ! save PHASE
            if(iPhase.EQ.0) iPhase = 1 ! primary alpha, by default
            SVARS1((iEL_proc-1)*NSVARS*NGAUSS + 304)=iPhase
         end if
      ENDDO
      
      CLOSE(201)
      END SUBROUTINE
      
      SUBROUTINE readHardening_OLD(SVARS1,NELX,NELST,NELNOX,NSVARS,NGAUSS)
      
      implicit none


      integer, intent(in) :: NELX,NELST,NELNOX,NSVARS,NGAUSS
      real(8), intent(inout) :: SVARS1(NSVARS*NELNOX*NGAUSS)
      
      integer :: I, iElem, iEL_proc, iPhase, NELEND
      real(8) :: euler(3), hard(47)
      
      NELEND = NELST + NELNOX - 1
      
      OPEN(800,FILE='statevariable.dat')
      DO I=1,NELX
         READ(800,*) iElem,euler(1),euler(2),euler(3),iPhase, & 
                     hard(1:47)
                     
         if(iElem.GE.NELST.AND.iElem.LE.NELEND) then
            iEL_proc = iElem - NELST + 1
            SVARS1((iEL_proc-1)*NSVARS*NGAUSS+305: &                   ! hardening params (GA) for HCP and trans slip families.
                   (iEL_proc-1)*NSVARS*NGAUSS+351)=hard(1:47)          ! See STATEVAR documentation.
            
            ! later on, consider moving GA's to material props vector.
         endif
      ENDDO
      CLOSE(800)
      END SUBROUTINE
      
!**********************************************************
      FUNCTION DOTPRD(V1,V2,NBC1,IEQNO)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION V1(1),V2(1),NBC1(1),NFAC(1),IEQNO(1)
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
       COMMON/DOF/NDOF

      NDF=NDOF
!
      DOTPRD=FDOT(V1,V2)
      NDF1=NDF+1
      IF(NBOUD1.EQ.0)GO TO 15
      DO 10 N=1,NBOUD1
      NODE=NBC1(NDF1*(N-1)+1)
      DO 12 I=1,NDF
      INDEX=NBC1(NDF1*(N-1)+1+I)
      IF(INDEX.NE.1)GO TO 12
      NN=NDF*(NODE-1)+I
      NN=IEQNO(NN)
      DOTPRD=DOTPRD-V1(NN)*V2(NN)
12    CONTINUE
10    CONTINUE
15    CONTINUE
      RETURN
      END FUNCTION

!***************************************************************
      ! calculates the stiffness matrix and load vector of an element
      ! for a given total deformation field
      ! deniz - modified - added output variable: 
      !         FE_appl(:) - external applied traction/body forces on the element nodes
      SUBROUTINE ELST1(JELEM,NNODE,XEUN,SKE,SKE_Patch,FE,FE_appl,TIME,DELTME, & 
          stressElement,ddsddeElement,FtauElement,dpElement, &
          volumeElement,volumeElement0, &
          PROPS1,NPROPS,NSVARS,NDOFEL,MDLOAD,MLVARX,NPREDF,MCRD, & 
          VTAU_E,VT_E,KINC,IFLAG,PNEWDT,umat_Err)
          
      use options, only: pressureDirectionZ, DEBUG_UMAT,DEBUG_RESID, &
                         EXPORT_stiffnessMatrix,SymmetrizeStiffnessMatrix, &
                         StiffnessMatrixFormulation,SM_delPdelF,SM_delSdelLnU
      use fBarPatch

      IMPLICIT NONE
      INCLUDE 'PARDIS.H'
      
      ! JELEM : the element, for which this ELST1 call is made
      INTEGER::  NTENS,NSTR
      PARAMETER (NTENS=6,NSTR=6)
      real(8), intent(out) :: stressElement(6),ddsddeElement(6,6),FtauElement(3,3),dpElement(6)
      real(8), intent(out) :: volumeElement,volumeElement0
      REAL(8), intent(out) :: FE(MDOFE),FE_appl(MDOFE),SKE(MDOFE,MDOFE)
      real(8), intent(out) :: SKE_Patch(MDOFE,MDOFE,MaxElemInPatch) ! FBar: non-local element couplings
      REAL(8)::  XEUN(MAXDIM,MAXELN)
      REAL(8)::  SVARS(MSV)  ! SVARS : state variables of the elements handled by this processor only
      REAL(8)::  U_TAU(MDOFE),U_T(MDOFE,1),V(MDOFE),A(MDOFE)
      REAL(8)::  TIME(2),PARAMS(1)
      REAL(8)::  PROPS(NPROPS),PROPS1(MAXPROPS)
      REAL(8)::  SVARS1(MAXSV),COORDS(MCRD,NNODE),TIME1(2)
      REAL(8)::  SVARS2(MAXSV),AMATRX(MDOFE,MDOFE)
      REAL(8)::  VTAU_E(MDOFE),VT_E(MDOFE)
!--------------- deniz - thermal -------------------------
      REAL(8):: TEMP(2), TEMP_0
!-----------------------------------------------------------------------
      REAL(8):: STRESS(6),DDSDDE(NTENS,NTENS)

      REAL(8):: UE_TAU(MCRD,NNODE)
      real(8):: FTau(3,3),Ft(3,3),FTau_FBar(3,3),Ft_FBar(3,3)
      REAL(8):: UE_T(MCRD,NNODE),DFGRD1_INV(3,3)
      REAL(8):: F1(MDOFE),SIGMA(3,3)
      
      REAL(8):: XE_T(MCRD,NNODE),XE_TAU(MCRD,NNODE)         !OK
      REAL(8):: SHP0(3,4),SHPTAU(3,4)                       !OK
      REAL(8):: SHP0_ST(3,4,MAXNP),SHPTAU_ST(3,4,MAXNP)     !OK
      REAL(8):: DVOL0_ST(MAXNP),DVOLTAU_ST(MAXNP),XB_AV(12) !OK
      
      REAL(8)::  DFGRD1_ST(3,3,MAXNP)
      REAL(8)::  DFGRD0_ST(3,3,MAXNP)
      REAL(8)::  XBAR_ST(2*MAXDIM,MDOFE,MAXNP)
      REAL(8)::  XC0(MAXDIM)
      REAL(8)::  XCTAU(MAXDIM)
      REAL(8)::  XJAC_ST(MAXNP)
      REAL(8)::  XBTR_ST(MDOFE,MAXNP)
      REAL(8)::  SKE1(MDOFE,MDOFE)

      REAL(8)::  SKE3(MDOFE,MDOFE)
      REAL(8)::    SKE2(MDOFE,MDOFE), STMAT(MAXDIM**2,MAXDIM**2)
      REAL(8)::  STMAT1(2*MAXDIM,2*MAXDIM),XMTEMP1(MAXDIM**2,MDOFE)
      REAL(8)::    RT_TAU(3,3),UT_TAU(3,3),DFGRDT(3,3)
      REAL(8):: CROT(2*MAXDIM,2*MAXDIM), FTEMP(MDOFE)
      REAL(8)::    XTM(2*MAXDIM,MDOFE),DROTM(2*MAXDIM,2*MAXDIM)
      REAL(8)::    XG_ST(MAXDIM**2,MDOFE,MAXNP),DFGRD0_INV(3,3)
      REAL(8)::    XBAR(NTENS,MDOFE)
      REAL(8)::  STATEV(MAXSTATEV),PRESSB(MAXEL),HARD(47)

      INTEGER::  IBELEM(MAXEL),ILF(MAXEL),IFACE_TET(12)
      INTEGER::  IFACE_BRICK(24),IBNODE(4)
      REAL(8):: XCG(3),UV1(3),UV2(3),UV3(3),UV4(3)
      REAL(8):: FLOAD(12)
      REAL(8):: SKE_LOAD(12,12)


      REAL(8):: XCYC(4,4,3),XJAC_TET(3,3),XJAC_TET_INV(3,3)
      REAL(8):: GRAD_TET_B(3,4),GRAD_TET(9,12),DUGRAD(9)
      REAL(8):: DFGRD1_TET(3,3),DFGRD0_TET(3,3),XJAC_TET0(3,3)
      REAL(8)::    XJAC_TET0_INV(3,3),GRAD_TET0_B(3,4)
      REAL(8)::  XCYC0(4,4,3),XBAR_TET(6,12)
      REAL(8):: TEMP2D(3,4),GRAD_TET0(9,12)

      REAL(8)::  X1(3),X2(3),X3(3),W1(3,3),W2(3,3),W3(3,3)

      INTEGER::  I,J,JELEM,NNODE,IELEM,IEL,iPhase

      REAL(8):: DELTME,PNEWDT,DET,XULR1,XULR2,XULR3
      INTEGER:: NPROPS,NSVARS,NDOFEL,MDLOAD
      INTEGER:: MLVARX,NPREDF,MCRD,KINC,IFLAG
      INTEGER:: JN, ISTA,IINP
      INTEGER:: IGAUSS_INIT,IFS,ICONV1,II
      REAL(8):: DVOLTAU
      REAL(8):: DVOL0,DTIME,DVOL0_TET,DVOL_TET
      INTEGER:: IST1,IPT,INUMJ,IN
      INTEGER:: IST,IS,IST2,IWR,ISTV_START
      INTEGER:: NELST,K,KSTEP,MCRD1,MDLOAD1
      INTEGER:: NDOFEL1,NPREDF1,NELEND,NELNOX,NEQ
      INTEGER:: NELX,NODE,NPLANE,NDF,NGAUSS
      INTEGER:: NSTATV,NSV,NSVARS1
      INTEGER:: NX
      REAL(8):: PR,TEMPVEC,TVOL0,TVOLTAU
      REAL(8):: XJBAR

!------------ UMAT ------------------------------------------

      REAL(8):: DTEMP, CELENT, DRPLDT,SSE
      REAL(8):: LAYER,  SCD, RPL, SPD
      INTEGER:: NDI, NSHR, KSPT

      CHARACTER*80 CMNAME
!----------------------------------------------------------------------
      COMMON/PRESS/PRESSB,IBELEM,ILF,IFACE_TET,IFACE_BRICK
      COMMON/ABQ/TIME1,NSVARS1,NDOFEL1,MDLOAD1,NPREDF1,SVARS1,MCRD1
      COMMON/STATEUPDATE/SVARS2
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/EL_PAR/NELST,NELEND,NELNOX
      ! deniz - thermal
      COMMON/THERM/TEMP,TEMP_0
      
      ! deniz - FBar      
      real(8) :: ksi_FBar(2)
      real(8) :: G_patch(3,4,MaxElemInPatch) 
      real(8) :: VolTau_patch(MaxElemInPatch)
      real(8) :: DetTau_patch(MaxElemInPatch)
      integer :: NElementsInPatch_e,ElementsInPatch_e(MaxElemInPatch)
      integer :: patchElemIdx,patchElemIdx_J,patchElemID,patchID,elemIdx
      logical :: foundElementInPatch
      
      ! deniz - FBar - delPdelF and Stiffness Matrix calculation
      real(8) :: delPdelF(3,3,3,3)
      real(8) :: matrixStiffnessTangent(maxdim*maxeln,maxdim*maxeln)
      real(8) :: matrixStiffnessTangentNonlocal(maxdim*maxeln,maxdim*maxeln,MaxElemInPatch)
      
      logical :: calcStiffnessMatrix
      
      ! error triggers
      logical :: listLocalUmatErrorsOccurred(7)
      logical :: trigger_umatError
      logical :: trigger_BroyDidNotConv
      logical :: trigger_cutBack
      logical :: trigger_reformK
      COMMON/triggers/trigger_umatError,listLocalUmatErrorsOccurred, &
                      trigger_BroyDidNotConv,trigger_cutBack,trigger_reformK

      ! processor info
      INTEGER::IPROC,NPROCS,NPROW,NPCOL,MPIINFSIZE
      INTEGER::MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      INTEGER::ISENDORRECV(MAXPROC)

      COMMON/PROCESSORINFO/IPROC,NPROCS,MPIINFSIZE,MPIINF, & 
                     IPROCINFO,ISENDORRECV
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT

      
      ! deniz - umat error flag
      integer, intent(out) :: umat_Err

      ! deniz
      integer :: iGrain

!------------------------------- INITIALIZE
      umat_Err = 0
      
      stressElement = 0.d0
      dpElement = 0.d0
      ddsddeElement = 0.d0
      volumeElement = 0.d0
      volumeElement0 = 0.d0

      DFGRD1_ST(:,:,:) = 0.D0
      DFGRD0_ST(:,:,:) = 0.D0
      XBAR_ST(:,:,:) = 0.D0
      XC0(:) = 0.D0
      XCTAU(:) = 0.D0
      XJAC_ST(:) = 0.D0
      XBTR_ST(:,:) = 0.D0
      SKE1(:,:) = 0.D0

      SKE3(:,:) = 0.D0
      SKE2(:,:) = 0.D0
      STMAT(:,:) = 0.D0
      STMAT1(:,:) = 0.D0
      XMTEMP1(:,:) = 0.D0
      RT_TAU(:,:) = 0.D0
      UT_TAU(:,:) = 0.D0
      DFGRDT(:,:) = 0.D0
      CROT(:,:) = 0.D0
      FTEMP(:) = 0.D0
      XTM(:,:) = 0.D0
      DROTM(:,:) = 0.D0
      XG_ST(:,:,:) = 0.D0
      DFGRD0_INV(:,:) = 0.D0
      XBAR(:,:) = 0.D0
      STATEV(:) = 0.D0
      HARD(:) = 0.D0

      XCG(:)=0.D0
      UV1(:)=0.D0
      UV2(:)=0.D0
      UV3(:)=0.D0
      UV4(:)=0.D0
      SKE_LOAD(:,:)=0.D0

      XCYC(:,:,:)=0.D0
      XJAC_TET(:,:)=0.D0
      XJAC_TET_INV(:,:)=0.D0
      GRAD_TET_B(:,:)=0.D0
      GRAD_TET(:,:)=0.D0
      DUGRAD(:)=0.D0
      DFGRD1_TET(:,:)=0.D0
      DFGRD0_TET(:,:)=0.D0
      XJAC_TET0(:,:)=0.D0
      XJAC_TET0_INV(:,:)=0.D0
      GRAD_TET0_B(:,:)=0.D0
      XCYC0(:,:,:)=0.D0
      XBAR_TET(:,:)=0.D0
      TEMP2D(:,:)=0.D0
      GRAD_TET0(:,:)=0.D0
      
!--------------------------------

      FE_appl(:) = 0.D0


      DO I=1,NPROPS
         PROPS(I)=PROPS1(I)
      ENDDO

      NSV=NSVARS*NGAUSS*(JELEM-1)-NSVARS*NGAUSS*(NELST-1)
      DO I=1,NSVARS*NGAUSS
         SVARS(I)=SVARS1(NSV+I)  ! SVARS : state variables of the element, at time t.
      ENDDO

      ! COORDS: stores the (undeformed) node positions of the element
      DO I=1,MCRD
         DO J=1,NNODE
            COORDS(I,J)=XEUN(I,J)
         ENDDO
      ENDDO

      DTIME=DELTME

!********************************
      DO I=1,MDOFE
         U_T(I,1)=0.0
      ENDDO
!*******************************

      DO I=1,MDOFE
         U_TAU(I)=VTAU_E(I)
         U_T(I,1)=VT_E(I)
      ENDDO

      KSTEP=1
      PNEWDT=1.0

      INUMJ=0


      UE_TAU(1:MCRD,1:NNODE)=RESHAPE(U_TAU,(/MCRD,NNODE/))
      UE_T(1:MCRD,1:NNODE)=RESHAPE(U_T(1:MCRD*NNODE,1),(/MCRD,NNODE/))

      XE_TAU(1:MCRD,1:NNODE)=COORDS(1:MCRD,1:NNODE) +  & 
                             UE_TAU(1:MCRD,1:NNODE)
      XE_T(1:MCRD,1:NNODE)=COORDS(1:MCRD,1:NNODE) +  & 
                           UE_T(1:MCRD,1:NNODE)

      TVOL0=0.D0
      TVOLTAU=0.D0


      SKE(1:MDOFE,1:MDOFE)=0.D0
      F1(1:MDOFE)=0.D0
      SKE_Patch = 0.d0

      DO IPT=1,NGAUSS
         CALL SHAPE_TET4(COORDS,IPT,SHP0,DVOL0,XC0)
         CALL SHAPE_TET4(XE_TAU,IPT,SHPTAU,DVOLTAU,XCTAU)


         DVOL0_ST(IPT)=DVOL0
         DVOLTAU_ST(IPT)=DVOLTAU

         TVOL0=TVOL0+DVOL0
         TVOLTAU=TVOLTAU+DVOLTAU

         SHP0_ST(1:MCRD,1:NNODE,IPT)=SHP0(1:MCRD,1:NNODE)
         SHPTAU_ST(1:MCRD,1:NNODE,IPT)=SHPTAU(1:MCRD,1:NNODE)

      ENDDO



      DO I=1,3
         DO J=1,4
            DO K=1,4
          XCYC(J,K,I)=XE_TAU(I,J)-XE_TAU(I,K)
          XCYC0(J,K,I)=COORDS(I,J)-COORDS(I,K)
            ENDDO
         ENDDO
      ENDDO


      XJAC_TET=XCYC(1:3,4,1:3)


      XJAC_TET0=XCYC0(1:3,4,1:3)

      CALL MATINV3_UEL(XJAC_TET,XJAC_TET_INV,DVOL_TET)

      CALL MATINV3_UEL(XJAC_TET0,XJAC_TET0_INV,DVOL0_TET)



      TEMP2D=0

      TEMP2D(1,1)=1.D0
      TEMP2D(2,2)=1.D0
      TEMP2D(3,3)=1.D0
      TEMP2D(1:3,4)=-1.D0

      GRAD_TET_B=MATMUL(XJAC_TET_INV,TEMP2D)
      GRAD_TET0_B=MATMUL(XJAC_TET0_INV,TEMP2D)

      !added by deniz
      !must initialize GRAD_TETs!
      GRAD_TET(:,:)=0.D0
      GRAD_TET0(:,:)=0.D0
      
      DO I=1,3
         IST1=(I-1)*3
         DO J=1,4
            IST2=(J-1)*3+I
            GRAD_TET(IST1+1:IST1+3,IST2)=GRAD_TET_B(1:3,J)
            GRAD_TET0(IST1+1:IST1+3,IST2)=GRAD_TET0_B(1:3,J)
         ENDDO
      ENDDO

      DUGRAD=MATMUL(GRAD_TET0,U_TAU(1:12))

      DFGRD1_TET=RESHAPE(DUGRAD,(/3,3/))

      DUGRAD=MATMUL(GRAD_TET0,U_T(1:12,1))

      DFGRD0_TET=RESHAPE(DUGRAD,(/3,3/))

      DO I=1,3
         DFGRD0_TET(I,I)=1.D0+DFGRD0_TET(I,I)
         DFGRD1_TET(I,I)=1.D0+DFGRD1_TET(I,I)
      ENDDO

      XBAR_TET(1,1:12)=GRAD_TET(1,1:12)
      XBAR_TET(2,1:12)=GRAD_TET(5,1:12)
      XBAR_TET(3,1:12)=GRAD_TET(9,1:12)
      XBAR_TET(4,1:12)=GRAD_TET(2,1:12)+GRAD_TET(4,1:12)
      XBAR_TET(5,1:12)=GRAD_TET(3,1:12)+GRAD_TET(7,1:12)
      XBAR_TET(6,1:12)=GRAD_TET(6,1:12)+GRAD_TET(8,1:12)

      CALL CALC_DFG_TET4(UE_T,NNODE,NDOFEL,MCRD,DFGRD0_ST,DVOL0_ST, & 
      XJAC_ST, SHP0_ST,XJBAR,1,NGAUSS)

      CALL CALC_DFG_TET4(UE_TAU,NNODE,NDOFEL,MCRD,DFGRD1_ST,DVOL0_ST, & 
       XJAC_ST,SHP0_ST,XJBAR,1,NGAUSS)


      CALL MAKEGRAD_TET4(XBAR_ST,XBTR_ST,XG_ST,SHPTAU_ST,DVOLTAU_ST, & 
       NNODE,MCRD,NDOFEL,NGAUSS)

      XBAR_ST(1:6,1:12,1)=XBAR_TET(1:6,1:12)

      F1(1:NDOFEL)=0.D0

      NSTATV=NSVARS
      DO IPT=1,NGAUSS
      
         Ftau(1:3,1:3)=DFGRD1_ST(1:3,1:3,IPT)
         Ft(1:3,1:3)=DFGRD0_ST(1:3,1:3,IPT)
         FtauElement = Ftau
         ! apply FBar stabilization
         elemIdx = JELEM-NELST+1
         patchID = patchIDfromElemID(JELEM)
         ksi_FBar(:)     = ksi_local(:,elemIdx)
         G_patch(:,:,:)  = G_local(:,:,:,elemIdx)
         VolTau_patch(:) = VolTau_local(:,elemIdx)
         NElementsInPatch_e = NElementsInPatch(patchID)
         ElementsInPatch_e(:) = ElementsInPatch(patchID,:)
         ! FBar patches for the element
         FTau_FBar = ksi_FBar(timeTau)*Ftau
         Ft_FBar   = ksi_FBar(timeT)*Ft

         KSTEP=1

         ISTV_START=NSVARS*(IPT-1)
         DO ISTA=1,NSTATV
            STATEV(ISTA)=SVARS(ISTV_START+ISTA)
         ENDDO

         ! deniz - modified
         ! OLD: Euler angles, phase IDs, initial hardening params.
         !      were read from disk into STATEV, here.
         ! NEW: moved back to loadvc, read into SVARS1
         
         ! UMAT DEBUG BLOCK
         !-----------------
         !Ft = 0.d0
         !Ft(1,1) = 1.d0
         !Ft(2,2) = 1.d0
         !Ft(3,3) = 1.d0
         !Ftau = 0.d0
         !Ftau(1,1) = 0.995
         !Ftau(2,2) = 0.995
         !Ftau(3,3) = 1.010
         !DTIME = 0.1d0
         !write(*,*) time(:), dtime 
         
         ! Ft: deform gradient at time t
         ! Ftau: trial deform gradient at time tau
         calcStiffnessMatrix = (IFLAG == 0)
         CALL UMAT_CPFEM_THERM(STRESS,STATEV,PROPS,DDSDDE,NTENS,delPdelF,NSTATV, & 
         NPROPS,PNEWDT,Ft_FBar,Ftau_FBar,TEMP(:),TEMP_0,KINC,DTIME,JELEM, &
         calcStiffnessMatrix,umat_Err) ! Deniz - thermal
         
         ! UMAT DEBUG BLOCK
         !-----------------
         !if (umat_Err/=0) write(*,*) 'umat error!'
         !write(*,*) JELEM, stress(:), pnewdt
         !write(*,*) statev(:)
         !open(5001,file='svars.txt')
         !write(5001,*) statev(1:nsvars)
         !close(5001)
         !call mpi_FINALIZE()
         
         if(DEBUG_UMAT)then  ! export UMAT debug info
            write(543,*) KINC, JELEM, intData(10), intData(11)
         endif
         
         ! UMAT_CPFEM solves the material model at the gauss point.
         ! updated state variables (at time tau) are stored in STATEV now.
         CALL assertBool((umat_Err.NE.0.AND.PNEWDT.LT.1.D0).OR. &
                         (umat_Err.EQ.0.AND.PNEWDT.GE.1.D0),.TRUE., &
                         'unhandled umat error state')
                    
         if(umat_Err.NE.0) return
         ! IF(PNEWDT.LT.1.D0)RETURN
         
         volumeElement = volumeElement + DVOLTAU_ST(IPT)
         volumeElement0 = volumeElement0 + DVOL0_ST(IPT)
         stressElement = stressElement + stress*DVOLTAU_ST(IPT)  ! will be normalized below outside the GaussPt loop
         ddsddeElement = ddsddeElement + ddsdde*DVOLTAU_ST(IPT)  ! will be normalized below outside the GaussPt loop
         dpElement = statev(327:332)

         STATEV(300)=DVOLTAU_ST(IPT)

         ! save back to the list of state variables for the element SVARS
         DO ISTA=1,NSTATV
            SVARS(ISTV_START+ISTA)=STATEV(ISTA)
         ENDDO


         ! **************************
         ! add internal forces to F1:
         ! **************************
         XBAR(1:NSTR,1:NDOFEL)=XBAR_ST(1:NSTR,1:NDOFEL,IPT) & 
           *DVOLTAU_ST(IPT)
         FTEMP=MATMUL(TRANSPOSE(XBAR(1:NSTR,1:NDOFEL)),STRESS(1:NSTR))
         F1(1:NDOFEL)=F1(1:NDOFEL)+FTEMP(1:NDOFEL)

         ! FBar - find the index of the element within its patch
         foundElementInPatch = .false.
         do patchElemIdx=1,NElementsInPatch_e
            if (jelem == ElementsInPatch_e(patchElemIdx)) then
               foundElementInPatch = .true.
               exit
            endif
         enddo
         if (.not.foundElementInPatch) then
            write(*,*) 'corrupt patch data structure'
            write(*,*) 'umat called for elemID ',jelem
            write(*,*) 'not found in patch of elements',ElementsInPatch_e(1:NElementsInPatch_e)
            CALL MPI_ABORT()
         endif

         ! STIFFNESS MATRIX CALCULATION
         !-----------------------------
         if (calcStiffnessMatrix) then
            if (StiffnessMatrixFormulation == SM_delSdelLnU) then  ! old way of calculating stiffness matrix ( ddsdde=delS/del(lnU) )

               IF(INUMJ.EQ.0)THEN
                  DO I=1,MCRD
                     SIGMA(I,I)=STRESS(I)
                  ENDDO
                  IF(MCRD.EQ.2)THEN
                     SIGMA(1,2)=STRESS(3)
                     SIGMA(2,1)=STRESS(3)
                  ENDIF
                  IF(MCRD.EQ.3)THEN
                     SIGMA(1,2)=STRESS(4)
                     SIGMA(1,3)=STRESS(5)
                     SIGMA(2,3)=STRESS(6)
                     DO I=1,3
                        DO J=I,3
                           SIGMA(J,I)=SIGMA(I,J)
                        ENDDO
                     ENDDO
                  ENDIF

                  STMAT(1:MCRD**2,1:MCRD**2)=0.D0

                  DO II=1,MCRD
                     IST=(II-1)*MCRD
                     DO I=1,MCRD
                        DO J=1,MCRD
                           STMAT(IST+I,IST+J)=SIGMA(I,J)
                        ENDDO
                     ENDDO
                  ENDDO


                  DO I=1,MCRD
                     STMAT1(I,I)=2*SIGMA(I,I)
                  ENDDO

                  IF(MCRD.EQ.3)THEN
                   STMAT1(1,4)=SIGMA(1,2)
                   STMAT1(4,1)=SIGMA(1,2)
                   STMAT1(1,5)=SIGMA(3,1)
                   STMAT1(5,1)=SIGMA(3,1)
                   STMAT1(2,4)=SIGMA(1,2)
                   STMAT1(4,2)=SIGMA(1,2)
                   STMAT1(2,6)=SIGMA(2,3)
                   STMAT1(6,2)=SIGMA(2,3)
                   STMAT1(3,6)=SIGMA(2,3)
                   STMAT1(6,3)=SIGMA(2,3)
                   STMAT1(3,5)=SIGMA(1,3)
                   STMAT1(5,3)=SIGMA(1,3)
                  ELSE
                   STMAT1(1,3)=SIGMA(1,2)
                   STMAT1(2,3)=SIGMA(2,1)
                   STMAT1(3,3)=0.5D0*(SIGMA(1,1)+SIGMA(2,2))
                  ENDIF
                  IF(MCRD.EQ.3)THEN
                   STMAT1(4,4)=0.5D0*(SIGMA(1,1)+SIGMA(2,2))
                   STMAT1(6,6)=0.5D0*(SIGMA(2,2)+SIGMA(3,3))
                   STMAT1(5,5)=0.5D0*(SIGMA(1,1)+SIGMA(3,3))
                   STMAT1(4,6)=0.5D0*SIGMA(3,1)
                   STMAT1(6,4)=0.5D0*SIGMA(3,1)
                   STMAT1(4,5)=0.5D0*SIGMA(3,2)
                   STMAT1(5,4)=0.5D0*SIGMA(3,2)
                   STMAT1(6,5)=0.5D0*SIGMA(1,2)
                   STMAT1(5,6)=0.5D0*SIGMA(1,2)
                  ENDIF
                  DO I=1,NDOFEL
                     DO J=1,NDOFEL
                        SKE1(I,J)=FTEMP(I)*XBTR_ST(J,IPT)
                     ENDDO
                  ENDDO

                  XMTEMP1(1:MCRD**2,1:NDOFEL)= & 
                        MATMUL(STMAT(1:MCRD**2,1:MCRD**2),XG_ST(1:MCRD**2,1:NDOFEL,IPT))
                  XMTEMP1(1:MCRD**2,1:NDOFEL)= & 
                        XMTEMP1(1:MCRD**2,1:NDOFEL)*DVOLTAU_ST(IPT)

                  SKE2(1:NDOFEL,1:NDOFEL)=MATMUL(TRANSPOSE(XG_ST(1:MCRD**2,1: & 
                        NDOFEL,IPT)),XMTEMP1(1:MCRD**2,1:NDOFEL))

                  CALL MATINV3_UEL(Ft,DFGRD0_INV,DET)

                  CALL MAT33(DFGRDT,Ftau,DFGRD0_INV,3)

                  CALL RUDCMP_UEL(DFGRDT,RT_TAU,UT_TAU)

                  CALL DROTMAT(RT_TAU,DROTM)

                  IF(MCRD.EQ.2)THEN
                  DO I=1,NSTR-1
                     DROTM(3,I)=DROTM(4,I)
                     DROTM(I,3)=DROTM(I,4)
                  ENDDO
                     DROTM(3,3)=DROTM(4,4)
                  ENDIF


                  CROT(1:NSTR,1:NSTR)=MATMUL(DDSDDE(1:NSTR,1:NSTR), & 
                    DROTM(1:NSTR,1:NSTR))
                  CROT(1:NSTR,1:NSTR)= & 
                    (CROT(1:NSTR,1:NSTR)-STMAT1(1:NSTR,1:NSTR))*DVOLTAU_ST(IPT)

                  XTM(1:NSTR,1:NDOFEL)=MATMUL(CROT(1:NSTR,1:NSTR), & 
                     XBAR_ST(1:NSTR,1:NDOFEL,IPT))

                  SKE3(1:NDOFEL,1:NDOFEL)= & 
                     MATMUL(TRANSPOSE(XBAR_ST(1:NSTR,1:NDOFEL,IPT)), & 
                     XTM(1:NSTR,1:NDOFEL))

                  SKE(1:NDOFEL,1:NDOFEL)=SKE(1:NDOFEL,1:NDOFEL)+ & 
                     SKE3(1:NDOFEL,1:NDOFEL)+SKE2(1:NDOFEL,1:NDOFEL)

               ENDIF
            
            elseif (StiffnessMatrixFormulation == SM_delPdelF) then ! the new way (delPdelF)               
            
               DetTau_patch = VolTau_patch * 6.d0
            
               CALL computeMatrixStiffnessTangent     &
               (                                            &
                  stress,                     &
                  delPdelF,                    &
                  Ftau,                       &
                  G_patch(:,:,:),            &
                  DetTau_patch,                 &
                  patchElemIdx,                                 &
                  NElementsInPatch_e,                         &
                  matrixStiffnessTangent,                   &
                  matrixStiffnessTangentNonlocal            &
               )   
               
               SKE(:,:) = matrixStiffnessTangent(:,:)
              
               if(DEBUG_RESID) THEN
                  
                  CALL printMatrix(matrixStiffnessTangent,maxdim*maxeln,maxdim*maxeln, &
                     'copied element stiffness matrix',970)
                     
               endif
              
               do patchElemIdx_J = 1,NElementsInPatch_e
                  SKE_Patch(:,:,patchElemIdx_J) = matrixStiffnessTangentNonlocal(:,:,patchElemIdx_J)
               enddo

               if (EXPORT_stiffnessMatrix .and. jelem == 1) then
                  write(972,*) 'element 1 - (local) stiffness matrix delPdelF'
                  do i = 1,ndofel
                     write(972,*) SKE(i,1:ndofel)
                  enddo
               endif
         
            endif ! stiffness matrix type
            
         endif
            
      enddo ! gauss points
      
      stressElement = stressElement / volumeElement
      ddsddeElement = ddsddeElement / volumeElement

      IF(IBELEM(JELEM).NE.0)THEN    ! an external applied surface traction was specified on this element

         IFS=IBELEM(JELEM) ! traction defined of the face IFS. see fem_uel for the numbering of faces
         IST=3*(IFS-1)
         IBNODE(1:3)=IFACE_TET(IST+1:IST+3)  ! IFACE_TET stores the local node IDs on each face of a TET elem

         XCG(1:3)=1.D0/4.D0*SUM(XE_TAU(1:3,1:4),2)
         UV1(1:3)=XCG-XE_TAU(1:3,IBNODE(1))
         UV2(1:3)=XE_TAU(1:3,IBNODE(2))-XE_TAU(1:3,IBNODE(1))
         UV3(1:3)=XE_TAU(1:3,IBNODE(3))-XE_TAU(1:3,IBNODE(1))
         
         CALL CROSS(UV4,UV3,UV2)
         CALL DOT3(TEMPVEC,UV1,UV4) !UV4 : face normal vector*2
         IS=1
         IF(TEMPVEC.LE.0)IS=-1
         IF(ILF(JELEM).EQ.1)THEN ! flag=1: loading applied by dload()
            CALL DLOAD(JELEM,PR,TIME)
         ENDIF

         PR=IS*PR
         FLOAD(:)=0.D0
         if(pressureDirectionZ) then  ! apply traction along the z-direction only
            FLOAD(3:3)=UV4(3)*PR/3.D0     ! nodal force vector for node 1, in the direction of face normal
            FLOAD(6:6)=UV4(3)*PR/3.D0     ! nodal force vector for node 2, in the direction of face normal
            FLOAD(9:9)=UV4(3)*PR/3.D0     ! nodal force vector for node 3, in the direction of face normal
         else                             ! apply traction along the surface normal (pressure)
            FLOAD(1:3)=UV4*PR/3.D0     ! nodal force vector for node 1, in the direction of face normal
            FLOAD(4:6)=UV4*PR/3.D0     ! nodal force vector for node 2, in the direction of face normal
            FLOAD(7:9)=UV4*PR/3.D0     ! nodal force vector for node 3, in the direction of face normal
         endif
         FLOAD(1:9)=0.5D0*FLOAD(1:9)! the face is triangle. divide by 2
         SKE_LOAD=0.D0

         X1=XE_TAU(:,IBNODE(1))  ! these are the positions of the three face nodes
         X2=XE_TAU(:,IBNODE(2))
         X3=XE_TAU(:,IBNODE(3))


         CALL MAKESPIN(W1,X3-X2)

         CALL MAKESPIN(W2,X1-X3)

         CALL MAKESPIN(W3,X2-X1)


         SKE_LOAD(1:9,1:9)=0.D0

         DO I=0,2
            SKE_LOAD(3*I+1:3*I+3,1:3)=W1
            SKE_LOAD(3*I+1:3*I+3,4:6)=W2
            SKE_LOAD(3*I+1:3*I+3,7:9)=W3
         ENDDO

         SKE_LOAD(1:9,1:9)=0.5D0*SKE_LOAD(1:9,1:9)*PR/3.D0

         ! if loads are fixed in a direction, then
         ! there's no contribution to stiffness matrix due to rotating/'follower' forces
         if (pressureDirectionZ) SKE_LOAD(:,:) = 0.D0         

         ! **************************************************
         ! save applied forces in a separate vector : FE_appl
         ! **************************************************

         DO I=1,3
            IN=IBNODE(I)
            F1(3*(IN-1)+1:3*(IN-1)+3) = F1(3*(IN-1)+1:3*(IN-1)+3) &
                           -FLOAD(3*I-2:3*I)      ! assemble the face nodal forces to the global force vector
            ! applied forces vector:
            FE_appl(3*IN-2:3*IN) = FLOAD(3*I-2:3*I)
            DO J=1,3
          JN=IBNODE(J)
          SKE(3*IN-2:3*IN,3*JN-2:3*JN)=SKE(3*IN-2:3*IN,3*JN-2:3*JN) &
                             -SKE_LOAD(3*I-2:3*I,3*J-2:3*J)
            ENDDO
         ENDDO

      ENDIF

      if (SymmetrizeStiffnessMatrix) then
         SKE(1:NDOFEL,1:NDOFEL)=0.5D0*(SKE(1:NDOFEL,1:NDOFEL)+ & 
            TRANSPOSE(SKE(1:NDOFEL,1:NDOFEL)))
      endif

      FE(1:NDOFEL)=-F1(1:NDOFEL) ! convention: Resiudal=out of balance forces. K*du = F results in G=F*du > 0

      IF(PNEWDT.LT.1.D0)RETURN

!****************************************************************************

      ! SVARS : state variables of the element, at time t.
      ! update the state variables for time Tau, (SVARS2 in common block)
      DO I=1,NSVARS*NGAUSS
         SVARS2(NSV+I)=SVARS(I)
      ENDDO

      END SUBROUTINE

!****************************************************
      FUNCTION DOT(A,B,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
      DOT=0.D0
      DO 100 I=1,N
  100 DOT=DOT+A(I)*B(I)
      RETURN
      END FUNCTION
!****************************************************
      FUNCTION FDOT(A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      FDOT=0.D0
      DO 100 I=1,NEQ
  100 FDOT=FDOT+A(I)*B(I)
      RETURN
      END FUNCTION
!****************************************************
      SUBROUTINE GET_SUBSCRIPT(IPROC,FNAME)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      CHARACTER*10 FNAME
      DIMENSION IN(10),INTMP(10)                                          !  <---------------------

      IPROC1=IPROC
      I=0

      DO WHILE(IPROC1.GT.9)
         I1=IPROC1/10
         I2=IPROC1-I1*10
         I=I+1
         IPROC1=I1
         IN(I)=I2
      ENDDO

      I=I+1
      IN(I)=IPROC1

      DO J=1,I
         INTMP(I-J+1)=IN(J)
      ENDDO

      DO J=1,I
         IN(J)=INTMP(J)
      ENDDO

      IF(I.EQ.1)THEN
         FNAME='_'//ACHAR(48+IN(1))//'.OUT'
      ELSEIF(I.EQ.2)THEN
         FNAME='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//'.OUT'
      ELSE
         FNAME='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//ACHAR(48+IN(3)) & 
         //'.OUT'
      ENDIF


      RETURN
      END SUBROUTINE
!---------------------------------------------------------------
! RETURNS TIME BASED ON DWELL LOAD SUCH THAT A TIME POINT IS PRESENT
! AT THE CHANGE POINTS
      SUBROUTINE TDWELL(FINTME,DELTME,XIDT,XMINDT,IWRITE_FLAG,NSEG,TSEG)
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION TSEG(NSEG+1)

      COMMON/LOAD_DWELL/TPERIOD,TRAMP,T_NODWELL,T_DWELL,SAMP,SMIN,UAMP, & 
             UMIN
      COMMON/LOAD_MONO/SLOAD,ST_RATE
      COMMON/ICYC_FLAG/ICYC


      XTOL=1D-6


      IF(TRAMP-FINTME.GT.XTOL)THEN

         FINTME1=FINTME+DELTME+XMINDT
         IF(FINTME1.GE.TRAMP)THEN
            DELTME=TRAMP-FINTME
            IWRITE_FLAG=1
         ENDIF



      ELSE

         FINTME2=FINTME-TRAMP

         REM=DMOD(FINTME2,TPERIOD)

         IF(DABS(REM).LT.XTOL.OR.DABS(REM-TPERIOD).LT.XTOL)THEN
            XNCYC=FINTME/TPERIOD
         ELSE
            XNCYC=(FINTME2-REM)/TPERIOD
         ENDIF

         T=FINTME2-XNCYC*TPERIOD


         DO I=1,NSEG

            TST=TSEG(I)
            TEND=TSEG(I+1)

            TEQ=T-TST


         IF((TEQ.GT.XTOL.OR.DABS(TEQ).LT.XTOL).AND.(TEND-T).GT.XTOL)THEN


               T1=T+DELTME+XMINDT
               T1EQ=T1-TEND

               IF(T1EQ.GT.XTOL.OR.DABS(T1EQ).LT.XTOL)DELTME=TEND-T
               IF(DABS(TEQ).LT.XTOL)DELTME=XIDT


            ENDIF


         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------------------
      ! WAR: stacked list of all PROC-local DVTAUs (gamma=u(k+1)-u(k) in Strang's notation)
      ! from the beginning of this time increment
      ! V: 
      SUBROUTINE BROY_PARALLEL(DVTAU,VAR,WAR,NUMUPD,NMPD1,MAXREF,IN_ST, & 
                   N_UPDATE,NEQ,STEP)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      INCLUDE 'mpif.h'
      DIMENSION DVTAU(NEQ),VAR(*),WAR(*)
      allocatable::V(:),W(:),V_RECV(:),W_RECV(:)

      ALLOCATE(V(NEQ),W(NEQ),V_RECV(NEQ),W_RECV(NEQ))
      V(1:NEQ)=0.D0
      W(1:NEQ)=0.D0
      V_RECV(1:NEQ)=0.D0
      W_RECV(1:NEQ)=0.D0

      STEP1=1.D0/STEP


      IF(NUMUPD.EQ.0)GOTO 30

      DO I=1,NUMUPD
         NN=(I-1)*N_UPDATE

         V(1:NEQ)=0.D0
         W(1:NEQ)=0.D0

         V(IN_ST+1:IN_ST+N_UPDATE)=VAR(NN+1:NN+N_UPDATE)
         W(IN_ST+1:IN_ST+N_UPDATE)=WAR(NN+1:NN+N_UPDATE)

         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         CALL MPI_ALLREDUCE(V,V_RECV,NEQ, & 
          MPI_DOUBLE_PRECISION, & 
          MPI_SUM,MPI_COMM_WORLD,IERROR)
         CALL MPI_ALLREDUCE(W,W_RECV,NEQ, & 
          MPI_DOUBLE_PRECISION, & 
          MPI_SUM,MPI_COMM_WORLD,IERROR)

         V(1:NEQ)=V_RECV(1:NEQ)
         W(1:NEQ)=W_RECV(1:NEQ)

         DW=DOT_PRODUCT(DVTAU(1:NEQ),W)

         DVTAU(1:NEQ)=DVTAU(1:NEQ)+DW*V(1:NEQ)
      ENDDO
 30   CONTINUE

      NUMUPD=NUMUPD+1

      W(1:NEQ)=0.D0
      W(IN_ST+1:IN_ST+N_UPDATE)=WAR(NMPD1+1:NMPD1+N_UPDATE)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(W,W_RECV,NEQ, & 
       MPI_DOUBLE_PRECISION, & 
       MPI_SUM,MPI_COMM_WORLD,IERROR)

      W(1:NEQ)=W_RECV(1:NEQ)

      V(1:NEQ)=STEP1*W(1:NEQ)-DVTAU(1:NEQ)

      FAC1=DOT_PRODUCT(V,W)

      FAC=1.D0/FAC1


      V(1:NEQ)=(W(1:NEQ)-V(1:NEQ))*FAC

      VAR(NMPD1+1:NMPD1+N_UPDATE)=V(IN_ST+1:IN_ST+N_UPDATE)

      DW=DOT_PRODUCT(DVTAU(1:NEQ),W)

      DVTAU(1:NEQ)=DVTAU(1:NEQ)+DW*V(1:NEQ)

      DEALLOCATE(V,W,V_RECV,W_RECV)

      RETURN
      END SUBROUTINE

      SUBROUTINE FROB_NORM(X,N,XNORM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N,N)

      XNORM=0.D0
      DO I=1,N
         DO J=1,N
            XNORM=XNORM+X(I,J)**2.D0
         ENDDO
      ENDDO

      XNORM=DSQRT(XNORM)
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------------------
! call this if the loading involves a ramping time.
! sets the time step to the end of the ramping time, if necessary
      SUBROUTINE T_CREEP_CSTRAINRATE(FINTME,DELTME,XMINDT,XMAXDT)
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON/LOAD_DWELL/TPERIOD,TRAMP,T_NODWELL,T_DWELL,SAMP,SMIN,UAMP,UMIN

      XTOL=1D-6

      IF(TRAMP-FINTME.GT.XTOL)THEN

         FINTME1=FINTME+DELTME+XMINDT
         IF(FINTME1.GE.TRAMP)THEN
            DELTME=TRAMP-FINTME
         ENDIF

      ENDIF

      RETURN
      END

!*******************************************************************
!     CREATES A SEQUENCE OF PROCESSOR PAIRS COMMUNICATING(SEND/RECV)
!     OBJECTIVE: REDUCES WAITING TIME FOR COMMUNICATING PROCESSOR
!*******************************************************************
      SUBROUTINE PROC_COMB(NPROC,MPIINF,MPIINFSIZE)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'PARDIS.H'

      DIMENSION MPIINF(MPIINFSIZEMAX)
      DIMENSION IPARR1(MAXCOMM,2),IPARR(2)

      allocatable::IGRPARR(:,:),IGRPNDX(:)

!*******************************************************************
!     NPROC(I)-NUMBER OF PROCESSORS
!     MPIINF(O)-PROCESSOR NUMBERS IN PAIRS. EX - (0123..), MEANS
!     01,23,.. ARE DIFFERENT COMMUNICATING PAIRS
!     MPIINFSIZE(O)-LENGTH OF ARRAY MPIINF
!*******************************************************************

      NSIZE=NPROC*(NPROC-1)

      IF(NSIZE.GT.MAXCOMM)THEN
        WRITE(*,*)'INCRESE MAXCOMM TO',NSIZE
        STOP
      ENDIF


      ILEVEL=0
      I2=NPROC
      DO WHILE(I2.GT.1)
         ILEVEL=ILEVEL+1
         I1=I2/2
         I2=I2-I1
      ENDDO


      ALLOCATE(IGRPARR(MAXCOMM,ILEVEL),IGRPNDX(ILEVEL))

      DO I=1,ILEVEL

         IF(I.EQ.1)THEN
           I1=NPROC/2
           I2=NPROC-I1
           IGRPNDX(I)=2
           IGRPARR(1,I)=I2
           IGRPARR(2,I)=I1

         ELSE
           ILOOP=IGRPNDX(I-1)
           ICNT=0

           DO J=1,ILOOP
              IP1=IGRPARR(J,I-1)
              IF(IP1.GT.1)THEN
                 I1=IP1/2
                 I2=IP1-I1
                 ICNT=ICNT+1
                 IGRPARR(ICNT,I)=I2
                 ICNT=ICNT+1
                 IGRPARR(ICNT,I)=I1
              ELSE
                 ICNT=ICNT+1
                 IGRPARR(ICNT,I)=0
              ENDIF
           ENDDO

           IGRPNDX(I)=ICNT

         ENDIF

      ENDDO

      IGNDX=0
      DO I=1,ILEVEL
         ITEMP=1
         IPROCBEGIN=0

         DO J=1,IGRPNDX(I)

            IF(ITEMP.LE.2)THEN
               IP1=IGRPARR(J,I)
               IF(IP1.EQ.0) THEN
                  IPROCBEGIN=IPROCBEGIN+1
                  GOTO 130
               ELSE
                  IPARR(ITEMP)=IGRPARR(J,I)
                  DO K=0,IPARR(ITEMP)-1
                     IPARR1(K+1,ITEMP)=IPROCBEGIN
                     IPROCBEGIN=IPROCBEGIN+1
                  ENDDO
                  ITEMP=ITEMP+1
               ENDIF
 130        ENDIF

            IF(ITEMP.GT.2)THEN
               ITEMP=1
               INDX1=1
               INDX2=2
               IF(IPARR(1).GT.IPARR(2)) THEN
                  INDX1=2
                  INDX2=1
               ENDIF

               DO K=0,IPARR(INDX2)-1
                  DO L=1,IPARR(INDX1)
                     IGNDX=IGNDX+1
                     MPIINF(IGNDX)=IPARR1(L,INDX1)
                     IGNDX=IGNDX+1
                     I1=IPARR(INDX2)
                     I2=L+K
                     I3=(I2-1)/I1
                     ILOC=I2-I3*I1
                     MPIINF(IGNDX)=IPARR1(ILOC,INDX2)
                  ENDDO
               ENDDO

            ENDIF

         ENDDO
      ENDDO

      MPIINFSIZE=IGNDX

      IF(MPIINFSIZE.GT.MPIINFSIZEMAX)THEN
         WRITE(*,*)'MPIINFSIZE',MPIINFSIZE,MPIINFSIZEMAX
         STOP
      ENDIF

      DEALLOCATE(IGRPARR,IGRPNDX)

      END SUBROUTINE
!***************************END*****************************************
!***********************************************************************
!     EACH PROCESSOR IDENTIFIES ITS SEQUENCE OF COMMUNICATION WITH
!     OTHER PROCESSORS FORM MPIINF ARRAY
!***********************************************************************
      SUBROUTINE ARRANGE_PROC(IPROCINFO,ISENDORRECV,MPIINF, & 
      			MPIINFSIZE,IPROC)

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION  MPIINF(MPIINFSIZEMAX),IPROCINFO(MAXPROC)
      DIMENSION ISENDORRECV(MAXPROC)
!***********************************************************************
!     IPROCINFO(O)-ARRAY CONTAINING PROCESSOR NUMBERS ARRANGED IN BEST
!     COMMUNICATION SEQUENCE
!     ISENDORRECV(O)-ARRAY CONTAING TAG WHICH DECIDES WHETHER PROC A WOULD
!     SEND AND THEN RECEIVE FROM B, OR VICE VERSA
!     TAG 1 - SEND FIRST
!     TAG 2 - RECEIVE FIRST
!     MPIINF(I)
!      MPIINFSIZE(I)
!     IPROC(I)
!***********************************************************************
      INDX=0

      DO I=1,MPIINFSIZE,2
	 IP1=MPIINF(I)
	 IP2=MPIINF(I+1)
	 IF(IP1.EQ.IPROC)THEN
	    INDX=INDX+1
	    IPROCINFO(INDX)=IP2
	    ISENDORRECV(INDX)=1
         ELSEIF(IP2.EQ.IPROC)THEN
            INDX=INDX+1
            IPROCINFO(INDX)=IP1
            ISENDORRECV(INDX)=2
	 ENDIF
      ENDDO

      END SUBROUTINE
!***********************END*********************************************

!***********************************************************************
!C     CREATES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE
!        modified by Kourosh for faster calculation
!***********************************************************************
      SUBROUTINE NODE_ELEM_CONNEC(IJK)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/NODECONNEC/IELEMNO,NINDX
      DIMENSION IJK(MNELX),IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************
      integer, allocatable :: IELEMNO_TMP(:,:),NINDX_TMP(:)
      integer :: nodeID

! ----------- slow way of doing this
!      ITEMP=1
!      DO I=1,NX
!         DO J=1,NELX
!            DO K=1,NODE
!               IF(IJK((J-1)*NODE+K).EQ.I)THEN
!                  IELEMNO(ITEMP)=J
!                  ITEMP=ITEMP+1
!                  GOTO 20
!               ENDIF
!            ENDDO
!   20    ENDDO
!         NINDX(I)=ITEMP
!      ENDDO
! ----------- faster way of doing this:
      ALLOCATE(IELEMNO_TMP(NX,MAXNODELEM),NINDX_TMP(NX))
      NINDX_TMP=0
      DO I=1,NELX
         DO J=1,NODE
            nodeID = IJK((I-1)*NODE+J)
            NINDX_TMP(nodeID)=NINDX_TMP(nodeID)+1
            IELEMNO_TMP(nodeID,NINDX_TMP(nodeID))=I
         ENDDO
      ENDDO
      
      KOUNT=1
      DO I=1,NX
         DO J=1,NINDX_TMP(I)
            IELEMNO(KOUNT)=IELEMNO_TMP(I,J)
            KOUNT=KOUNT+1
         ENDDO
         NINDX(I)=KOUNT
      ENDDO
      
      DEALLOCATE(IELEMNO_TMP,NINDX_TMP)
!------------------------------------------
      
      
      END SUBROUTINE
      
!***********************************************************************
!C     READS LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE readNODE_ELEM_CONNEC(strConnectFile,success)
      implicit none
      INCLUDE 'PARDIS.H'
      character(len=*), intent(in) :: strConnectFile      
      logical, intent(out) :: success
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/NODECONNEC/IELEMNO,NINDX
      integer :: NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      integer :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, nNodes, I, lastIndex

      success = .false.
      
      inquire(file=trim(strConnectFile),exist=fileExists)
      if(.NOT.fileExists) return
      
      open(716,file=trim(strConnectFile))
      read(716,*,iostat=error) nNodes
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      lastIndex = 1
      do iNode=1,nNodes
         read(716,*,iostat=error) nElem
         read(716,*,iostat=error) (IELEMNO(I),I=lastIndex,lastIndex+nElem-1)
         lastIndex = lastIndex + nElem
         NINDX(iNode) = lastIndex
      enddo
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      
      success = .true.
      close(716)
      
      END SUBROUTINE
      
!***********************************************************************
!C     WRITES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE writeNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX, &
                                       NELX,NX,MAXNODELEM,MAXNODE)
                                       
      implicit none
      character(len=*), intent(in) :: strConnectFile      

      integer, intent(in) :: NX,NELX,MAXNODELEM,MAXNODE
      integer, intent(in) :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, I, lastIndex
      integer :: elemStart, elemEnd

      
      open(716,file=trim(strConnectFile))
      write(716,*,iostat=error) NX
      lastIndex = 1
      do iNode=1,NX
         elemStart=1
         if(iNode.ne.1) elemStart=NINDX(iNode-1)
         elemEnd = NINDX(iNode)-1
         nElem = elemEnd - elemStart + 1
         
         write(716,*,iostat=error) nElem
         write(716,*,iostat=error) (IELEMNO(iElem),iElem=elemStart,elemEnd)
         lastIndex = lastIndex + nElem
      enddo
      
      close(716)
      
      END SUBROUTINE


      SUBROUTINE readElementNonLocalList(strNonlocalFile,elemNonLocalList,elemNonLocalCount, &
                                          elemNonLocalWeights,damageCharacteristicLength, &
                                          MAXNONLOCALELEMS,NELX,success)
      implicit none

      character(len=*), intent(in) :: strNonlocalFile      
      integer, intent(out):: elemNonLocalList(MAXNONLOCALELEMS,NELX)
      integer, intent(out):: elemNonLocalCount(NELX)
      real(8), intent(out):: elemNonLocalWeights(MAXNONLOCALELEMS,NELX)
      real(8), intent(out):: damageCharacteristicLength
      integer, intent(in) :: NELX,MAXNONLOCALELEMS
      logical, intent(out) :: success

      !Local
      logical :: fileExists
      integer :: error, iElem, jElem, nElem, nElemCount, iNode, nElems, I, lastIndex

      success = .false.
      
      elemNonLocalList = 0
      elemNonLocalCount = 0
      elemNonLocalWeights = 0.d0
            
      inquire(file=trim(strNonlocalFile),exist=fileExists)
      if(.NOT.fileExists) then
         return
      endif
      
      open(716,file=trim(strNonlocalFile))
      read(716,*) damageCharacteristicLength
      read(716,*,iostat=error) nElems
      if (error.ne.0 .or. nElems.NE.NELX) then
         close(716)
         return
      endif
      
      do iElem=1,nElems
      
         read(716,*,iostat=error) nElemCount
         if (nElemCount > MAXNONLOCALELEMS) then
            write(*,*) 'INCREASE MAXNONLOCALELEMS to ',nElemCount
            STOP
         endif
         elemNonLocalCount(iElem) = nElemCount
         read(716,*,iostat=error) (elemNonLocalList(I,iElem),I=1,nElemCount)
         read(716,*,iostat=error) (elemNonLocalWeights(I,iElem),I=1,nElemCount)
                                  
      enddo
      
      if (error.ne.0) then
         close(716)
         return
      endif
      
      success = .true.
      close(716)
      
      END SUBROUTINE


!*************************END*******************************************
!***********************************************************************
!C     INITIALIZES SIZE OF MODIFIED SPARSE ROW FORMAT
!***********************************************************************
      SUBROUTINE INIT_SIZE_SLU(IJK,IPROC,N_UPDATE,IN_ST,NNZ_LOC)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/NODECONNEC/IELEMNO,NINDX
       DIMENSION IJK(MNELX),IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
       DIMENSION INCOMP(MAXNODE)
!***********************************************************************
! C     IJK(I)
! C     N_UPDATE(I)
! C     NDF(I)
!***********************************************************************

       INDX=1

       DO I=1,N_UPDATE

          NNO=(IN_ST+I-1)/NDF+1
          ICOMP=0

          IF(NNO.EQ.1)THEN
            ISTART=1
          ELSE
            ISTART=NINDX(NNO-1)
          ENDIF

          IEND=NINDX(NNO)-1
          DO  J=ISTART,IEND

              IELNO=IELEMNO(J)

              DO K=1,NODE

                 IMATCH=0
                 NN=IJK((IELNO-1)*NODE+K)

                 DO LL=1,ICOMP
                    IF(NN.EQ.INCOMP(LL))THEN
                       IMATCH=1
                       GOTO 30
                    ENDIF
                 ENDDO

  30             IF(IMATCH.EQ.0)THEN

                    ICOMP=ICOMP+1
                    INCOMP(ICOMP)=NN

                 ENDIF

             ENDDO
          ENDDO



          DO J=1,NX
             DO K=1,ICOMP
                IF(J.EQ.INCOMP(K))THEN
                   DO L=1,NDF
                      INDX=INDX+1
                   ENDDO
                ENDIF
             ENDDO
          ENDDO


      ENDDO

      NNZ_LOC=INDX-1

      END SUBROUTINE
!*************************END*******************************************
!***********************************************************************
!C     INITIALIZES MODIFIED SPARSE ROW FORMAT
!***********************************************************************
      SUBROUTINE INIT_MSR_SLU(IJK,IPROC,N_UPDATE,IN_ST,IROWPTR, & 
                   IBINDX,NNZ_LOC)
!      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT NONE
      INCLUDE 'PARDIS.H'

      INTEGER:: IJK(MNELX),IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      INTEGER:: INCOMP(MAXNODE)
!***********************************************************************
! C     IJK(I)
! C     N_UPDATE(I)
! C     IUPDATE(I)-CONTAINS ROWS TO BE UPDATED BY THIS PROCESSOR
! C     IROWPTR(O) - stores the data-index of first elements of each row. 
! C                  IROWPTR(:) = index of 1st el ol 1st row, index of 1st el of 2nd row... 
! C                  usage: k = IROWPTR(i) to IROWPTR(i+1)-1, VAL(k) traverses elements on row i
! C     IBINDX(O)-ARRAY CONTAINING NON-ZERO COLUMN INDICES FOR EACH ROW
! C     CONSISTENT WITH AZTEC
! C     VAL(O)-ARRAY WITH ELEMENTS INITIALIZED TO ZERO
      ! IBINDX has all the non-zero matrix coefficients for the equations of this CPU.
      ! IBINDX is ordered. (all DOFS coupled to the 1st DOF of this CPU, all DOFS coupled to the 2nd DOF of this CPU, all DOFS coupled to the 3rd DOF of this CPU, ... )
      ! INDX: total # of couplings from this CPU, size of IBINDX
!***********************************************************************
      INTEGER:: IBINDX(NNZ_LOC),IROWPTR(N_UPDATE+1)
      INTEGER:: N_UPDATE,IMATCH,IEND,IELNO
      INTEGER:: ICOMP,NNO,K
      INTEGER::ISTART, NN,INDX,LL
      INTEGER:: IN_ST,NDF,I,NODE
      INTEGER:: J,NX,L,IPROC
      INTEGER:: NNZ_LOC,NELX,NEQ,NPLANE,NGAUSS

      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/NODECONNEC/IELEMNO,NINDX

      INDX=1

      DO I=1,N_UPDATE         ! traverse all the nodes/DOFs handled by this processor

         NNO=(IN_ST+I-1)/NDF+1
         ICOMP=0
         IROWPTR(I)=INDX         !INDX = index of the coefficient in the 1-D data array

          IF(NNO.EQ.1)THEN
            ISTART=1
          ELSE
            ISTART=NINDX(NNO-1)  ! NINDX(i): the index of connectivity data of node a, in IELEMNO. IELEMNO(NINDX(i) + ...) = connected elemIDs
          ENDIF
          IEND=NINDX(NNO)-1

          DO  J=ISTART,IEND   ! traverse all elements connected to the node NNO. DOFs (I) of the node NNO are handled by this processor.
                              ! will now the node coupling information to determine positions of non-zero coefficients

              IELNO=IELEMNO(J) 

              DO K=1,NODE

                 IMATCH=0
                 NN=IJK((IELNO-1)*NODE+K)    ! found coupling: node NN is coupled to node NNO, (through element IELNO)

                 DO LL=1,ICOMP               ! check if already added to the list of couplings for this DOF (IN_ST+I)
                    IF(NN.EQ.INCOMP(LL))THEN
                       IMATCH=1
                       GOTO 30
                    ENDIF
                 ENDDO

  30             IF(IMATCH.EQ.0)THEN

                    ICOMP=ICOMP+1
                    INCOMP(ICOMP)=NN         ! add the coupled node to the list of couplings for this DOF (IN_ST+I)

                 ENDIF

             ENDDO
          ENDDO
          
          ! ICOMP is the total # of nodes coupled to this DOF (IN_ST+I)
          ! INCOMP(:) is the list of nodes coupled to this DOF (IN_ST+I)
          
          ! traverse all the DOFs (NX -> (J-1)*NDF+L) and for the ones that are coupled to this DOF (INCOMP), mark the indices of non-zero coefficients in K(:,:) on IBINDX(:)
          DO J=1,NX
             DO K=1,ICOMP
                IF(J.EQ.INCOMP(K))THEN ! coupled node. mark all the DOFs of this node. store their EQN #s in order into IBINDX
                   DO L=1,NDF
                      IBINDX(INDX)=(J-1)*NDF+L
                      INDX=INDX+1
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          
          ! now all DOFs coupled to the DOF (IN_ST+I-1), are added to the end of IBINDX

      ENDDO
      
      IROWPTR(N_UPDATE+1)=INDX
      ! deniz - modified this line for a more conventional indexing structure
      ! OLD: IROWPTR(N_UPDATE+1)=INDX-1
      ! NEW: IROWPTR(N_UPDATE+1)=INDX
      ! now we can simply use:
      ! IS=IROWPTR(I)
      ! IE=IROWPTR(I+1)-1
      ! without the correction: IF(I.EQ.N_UPDATE)IE=IROWPTR(I+1)..
      
      END SUBROUTINE
      
   ! deniz - initializes modified sparse row format for the constraint/support reaction stiffness  matrix
      SUBROUTINE INIT_SIZE_RXN(N_UPDATE_rxn,ISIZE_rxn,IN_ST,N_UPDATE, & 
                               IROWPTR,IBC)
      implicit none
      include 'PARDIS.H'
      
      integer, intent(out):: N_UPDATE_rxn, ISIZE_rxn
      integer, intent(in) :: IN_ST, N_UPDATE
      integer, intent(in) :: IROWPTR(N_UPDATE+1)
      integer, intent(in) :: IBC(MDOFX)
      
      integer :: j_row

      N_UPDATE_rxn = 0
      ISIZE_rxn = 0
      
      do j_row=IN_ST+1,IN_ST+N_UPDATE
         if(IBC(j_row).EQ.1) then
            N_UPDATE_rxn = N_UPDATE_rxn + 1
            ISIZE_rxn = ISIZE_rxn  & 
                  + IROWPTR(j_row-IN_ST+1)-IROWPTR(j_row-IN_ST)
         endif
      enddo
      
      END SUBROUTINE
      
      SUBROUTINE INIT_MSC_rxn(ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn, & 
                              ISIZE_rxn,IN_ST,N_UPDATE,NNZ_LOC, & 
                              IROWPTR,IBINDX,IBC)
      implicit none
      include 'PARDIS.H'

      integer, intent(out):: ICOLPTR_rxn(N_UPDATE_rxn+1)
      integer, intent(out):: IBINDX_rxn(ISIZE_rxn)      
      integer, intent(in) :: N_UPDATE_rxn, ISIZE_rxn
      integer, intent(in) :: IN_ST, N_UPDATE, NNZ_LOC
      integer, intent(in) :: IROWPTR(N_UPDATE+1)
      integer, intent(in) :: IBINDX(NNZ_LOC)
      integer, intent(in) :: IBC(MDOFX)

      integer :: k_rxn, j_row, i_col, n_row_entries
      integer :: k_sta, k_end
      
      IBINDX_rxn(:) = 0
      ICOLPTR_rxn(:) = 0
      
      i_col = 1
      k_rxn = 1
      
      ICOLPTR_rxn(i_col) = 1
      
      do j_row=IN_ST+1,IN_ST+N_UPDATE
         if(IBC(j_row).EQ.1) then
         
            k_sta=IROWPTR(j_row-IN_ST)
            k_end=IROWPTR(j_row-IN_ST+1)-1
            n_row_entries = k_end - k_sta + 1

            ! move the row structure from K as a column structure to K_rxn
            i_col = i_col + 1
            ICOLPTR_rxn(i_col) = ICOLPTR_rxn(i_col-1) + n_row_entries
            IBINDX_rxn(k_rxn:k_rxn+n_row_entries-1) = IBINDX(k_sta:k_end)
            k_rxn = k_rxn + n_row_entries
         endif
      enddo
      
      ! ICOLPTR_rxn has the proper indexing structure, (unlike IROWPTR, which has to be fixed)
      ! so I'm not going to assign ICOLPTR_rxn(i_col) = k_rxn - 1
      
      END SUBROUTINE
      
      ! extracts the list of constrained DOF/EQN numbers on this proc --> EQ_uBC_proc.
      SUBROUTINE INIT_uBC_proc(EQ_uBC_proc,N_UPDATE_rxn,IBC,   &
                               IN_ST,N_UPDATE)
      implicit none
      include 'PARDIS.H'

      integer, intent(out):: EQ_uBC_proc(1:N_UPDATE_rxn)
      integer, intent(in) :: N_UPDATE_rxn,IN_ST,N_UPDATE
      integer, intent(in) :: IBC(MDOFX)

      integer :: iDOF, i_uBC_proc
      
      i_uBC_proc=0
  
      do iDOF=IN_ST+1,IN_ST+N_UPDATE
         if(IBC(iDOF).EQ.1) then
            i_uBC_proc=i_uBC_proc+1
            EQ_uBC_proc(i_uBC_proc) = iDOF
         endif
      enddo
      
      CALL assertInt(i_uBC_proc,N_UPDATE_rxn,   &
      'assertion error: i_uBC_proc!=N_UPDATE_rxn at INIT_uBC_proc')
      END SUBROUTINE
      
      ! extracts the imposed values of constrained DOF/EQN on this proc --> uBC_proc.
      SUBROUTINE extract_uBC_proc(dU_BCs,uBC_proc,EQ_uBC_proc,N_UPDATE_rxn,NEQ)
      implicit none
      include 'PARDIS.H'

      real(8), intent(in) :: dU_BCs(NEQ)
      real(8), intent(out):: uBC_proc(1:N_UPDATE_rxn)
      integer, intent(in) :: EQ_uBC_proc(1:N_UPDATE_rxn)
      integer, intent(in) :: N_UPDATE_rxn, NEQ

      integer :: iDOF, i_uBC_proc
  
      do i_uBC_proc=1,N_UPDATE_rxn
         iDOF = EQ_uBC_proc(i_uBC_proc)
         uBC_proc(i_uBC_proc) = dU_BCs(iDOF)
      enddo
      END SUBROUTINE
      
      
      SUBROUTINE MATMUL_ColSparse(MAT,vec,res,IBINDX,ICOLPTR, & 
                                  n_cols,n_rows,n_entries)
      implicit none
      real(8), intent(in) :: MAT(n_entries)
      real(8), intent(in) :: vec(n_cols)
      real(8), intent(out):: res(n_rows)
      integer, intent(in) :: IBINDX(n_entries)      
      integer, intent(in) :: ICOLPTR(n_cols+1)
      integer, intent(in) :: n_cols, n_rows, n_entries
      
      integer :: i_col, i_row, k, k_sta, k_end
      res = 0.D0
      
      do i_col = 1,n_cols
      
         k_sta = ICOLPTR(i_col)
         k_end = ICOLPTR(i_col+1)-1 ! COLPTR has the proper index at ICOLPTR(n_cols+1)
         
         do k = k_sta, k_end
            i_row = IBINDX(k)
            res(i_row) = res(i_row) + MAT(k)*vec(i_col)
         enddo
               
      enddo      
      
      END SUBROUTINE

      ! converts IROWPTR and IBINDX sparse matrix structure vectors
      ! from Fortran (1-indexed) to C (0-indexed) format
      ! SuperLU is a C library. C is zero-indexed.
      SUBROUTINE convertSPARSE_f2C(IROWPTR_F,IBINDX_F,IROWPTR_C,IBINDX_C, &
                                   N_UPDATE,NNZ_LOC)
      implicit none
      integer, intent(in) :: IROWPTR_F(N_UPDATE+1)
      integer, intent(in) :: IBINDX_F(NNZ_LOC)      
      integer, intent(out):: IROWPTR_C(N_UPDATE+1)
      integer, intent(out):: IBINDX_C(NNZ_LOC)      
      integer, intent(in) :: N_UPDATE, NNZ_LOC
      
      IROWPTR_C(1:N_UPDATE+1) = IROWPTR_F(1:N_UPDATE+1) - 1
      IBINDX_C(1:NNZ_LOC) = IBINDX_F(1:NNZ_LOC) - 1
      
      END SUBROUTINE


!*************************END*******************************************
!***********************************************************************
!C     FINDS PROC NO ASSOCIATED WITH A NODE CONSIDERING LINEAR
!C     PARTITIONING OF grid POINTS
!***********************************************************************
      SUBROUTINE FINDPROC(NNO,NPROCS,NODEINPROC)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
!***********************************************************************
!C     NNO(I)-NODE NUMBER
!C     NPROCS-NUMBER OF PROCESSORS
!C     NODEINPROC(O)-PROCESSOR NUMBER TO WHICH NODE IS ASSOCIATED
!***********************************************************************

      I1=NX/NPROCS+1
      I2=NX/NPROCS

      IA=NX-I2*NPROCS
      IB=NPROCS-IA


      IF(NNO.LE.I1*IA)THEN
         NODEINPROC=(NNO-1)/I1
      ELSE
         NNO1=NNO-I1*IA
         NODEINPROC=IA+(NNO1-1)/I2
      ENDIF


      END SUBROUTINE

! added by deniz
! *******************************************************************
! ***************** LINEAR PARTITIONING OF ELEMENTS *****************

      SUBROUTINE partitionElements(NELX,NPROCS,PROCID, & 
                                   NELST,NELEND,NELNOX)
      implicit none
      integer, intent(in) :: NELX, NPROCS, PROCID
      integer, intent(out) :: NELST,NELEND,NELNOX
      
      integer :: I1, I2, IA

      if (NELX < NPROCS) then ! each proc gets one.
         if (PROCID+1 <= NELX) then
            NELST = PROCID+1
            NELEND = PROCID+1
            NELNOX = 1
         else
            NELST = 0
            NELEND = 0
            NELNOX = 0
         endif
         return
      endif
      
      I1=NELX/NPROCS+1
      I2=NELX/NPROCS

      IA=NELX-I2*NPROCS

      IF((PROCID+1).LE.IA)THEN
         NELST=PROCID*I1+1
         NELNOX=I1
      ELSE
         NELST=IA*I1+(PROCID-IA)*I2+1
         NELNOX=I2
      ENDIF

      ! range of the elements handled by each processor: NELST -- NELEND
      NELEND=NELST+NELNOX-1
      END SUBROUTINE
!*************************END*******************************************
      SUBROUTINE UPDATE_LOCAL_ARRAY(IJK, NROWS,RECVVAL,RECVRHS,IROWARR, & 
         IELARR,IRECVRINDX,IPROC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL, & 
         RESID,NNZ_LOC,IFLAG)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION IJK(MNELX)
      DIMENSION IROWARR(MDOFPROC),IELARR(MDOFPROC), & 
                RECVVAL(MDOFPROC*MAXDIM*8),RECVRHS(MDOFPROC), & 
                IRECVRINDX(MDOFPROC)
      DIMENSION ICOLMARR(MDOFE)
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      DIMENSION IBINDX(NNZ_LOC),VAL(NNZ_LOC),RESID(N_UPDATE), & 
       IROWPTR(N_UPDATE+1)



      DO I=1,NROWS
         IROWNO=IROWARR(I)
         IELNO=IELARR(I)


         DO J=1,NODE
            NNO=IJK((IELNO-1)*NODE+J)
            DO K=1,NDF
               ITEMP=(NNO-1)*NDF+K
               ICOLMARR((J-1)*NDF+K)=ITEMP
            ENDDO
         ENDDO

         IROW1=IROWNO-IN_ST
         IS=IROWPTR(IROW1)
         IE=IROWPTR(IROW1+1)-1

         ! IF(IROW1.EQ.N_UPDATE)IE=IROWPTR(N_UPDATE+1) -- deniz modified into conventional form

         IF(IFLAG.EQ.0)THEN

            DO K=IS,IE
               ICOLM=IBINDX(K)
               DO KK=1,NODE*NDF
                  IF(ICOLM.EQ.ICOLMARR(KK))THEN
                     VAL(K)=VAL(K)+RECVVAL((I-1)*NDF*NODE+KK)
                     GOTO 120
                  ENDIF
               ENDDO
 120        ENDDO
         ENDIF


         RESID(IROW1)=RESID(IROW1)+RECVRHS(I)

      ENDDO



      END SUBROUTINE
      
      SUBROUTINE UPDATE_couplings_FBar(IJK,NROWS,RECVVAL_Patch,RECVRHS,IROWARR, & 
         IELARR,IRECVRINDX,IPROC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL, & 
         NNZ_LOC,IFLAG)
      
      use fBarPatch
      
      implicit none
      INCLUDE 'PARDIS.H'
      integer, intent(in) :: NROWS,IPROC,N_UPDATE,IN_ST,NNZ_LOC,IFLAG
      integer, intent(in) :: IJK(MNELX)
      real(8), intent(inout) :: VAL(NNZ_LOC)
      integer, intent(in) :: IROWARR(MDOFPROC),IELARR(MDOFPROC)
      real(8), intent(in) :: RECVVAL_Patch(MDOFPROC*MAXDIM*8*MaxElemInPatch)
      real(8), intent(in) :: RECVRHS(MDOFPROC)
      integer, intent(in) :: IRECVRINDX(MDOFPROC)
      integer, intent(in) :: IBINDX(NNZ_LOC),IROWPTR(N_UPDATE+1)

      ! locals
      integer :: patchID,patchElemID,patchElemIdx
      integer :: ICOLMARR(MDOFE)
      integer :: I,J,K,KK,IS,IE,IELNO,IROW1,IROWNO,ICOLM,NNO,ITEMP
      ! common variables
      integer :: NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      
      IF (IFLAG /= 0) return

      DO I=1,NROWS
         IROWNO=IROWARR(I)
         IELNO=IELARR(I)      ! this is 'my' element in the patch.
                              ! loop over the others

         patchID = patchIDfromElemID(IELNO)
         do patchElemIdx = 1,NElementsInPatch(patchID)
            patchElemID = ElementsInPatch(patchID,patchElemIdx)
            if (patchElemID == IELNO) cycle ! this is 'me' - skip.

            DO J=1,NODE
               NNO=IJK((patchElemID-1)*NODE+J)  ! nodes of the coupled patch element
               DO K=1,NDF
                  ITEMP=(NNO-1)*NDF+K
                  ICOLMARR((J-1)*NDF+K)=ITEMP
               ENDDO
            ENDDO

            IROW1=IROWNO-IN_ST
            IS=IROWPTR(IROW1)
            IE=IROWPTR(IROW1+1)-1

            DO K=IS,IE
               ICOLM=IBINDX(K)
               DO KK=1,NODE*NDF
                  IF(ICOLM.EQ.ICOLMARR(KK))THEN
                     VAL(K) = VAL(K) + &
                           RECVVAL_Patch((I-1)*NODE*NDF*MaxElemInPatch+ &
                                         (patchElemIdx-1)*NODE*NDF+KK)
                     GOTO 120
                  ENDIF
               ENDDO
120         ENDDO
         enddo
      ENDDO

      END SUBROUTINE

!*************************END*******************************************
! -------------------------------------------------------------------
      ! deniz.  modification!
      ! domainRange imported from common block
      ! it stores the dimensions of the domain
      subroutine disp(u,dir,nodePos,loadType,functionDisp,strain, &
                      macroF_tau,uCorner_tau,ksiCorner_init,domainLen0)
      
      implicit none
      
      real(8), intent(out):: u(3)
      integer, intent(in) :: dir ! direction of displacement DOF
      real(8), intent(in) :: nodePos(3) ! underformed position of the node, measured from the center of the domain
      integer, intent(in) :: loadType
      real(8), intent(in) :: functionDisp,strain
      real(8), intent(in) :: macroF_tau(3,3),uCorner_tau(3,8),ksiCorner_init(3,8)
      real(8), intent(in) :: domainLen0(3)

      real(8) :: du(3),ksi(3)
      real(8) :: shapeFuncHEX8(8)
      real(8), parameter :: Identity(3,3) = &
         reshape((/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /), (/3,3/))

      u(3)=0.0 ! these are not used. see caller of disp()
      u(2)=0.0

      if(loadType==1 .or.  &  ! --------- CREEP LOADING
         loadType==3 .or.  &  ! --------- traction controlled cyclic loading
         loadType==5 .or.  &  ! --------- traction-time prescribed
         loadType==7) then    ! --------- macroscopic stress-time prescribed
      
         ! not applicable to this DOF
         write(*,*) 'disp() called for traction controlled loading'

         stop

      elseif(loadType==2) then   ! --------- CSR loading

         u(1) = domainLen0(dir)*strain

      elseif(loadType==4 .or. &  ! --------- displacement-controlled cyclic loading
             loadType==6) then   ! --------- disp-time provided
             
         u(1) = functionDisp
         
      elseif(loadType==8 .or. &  ! --------- strain-time provided
             loadType==9 .or. &  ! --------- VM-strain controlled, biaxial stress loading
             loadType==10) then  ! --------- VM-strain controlled, uniaxial loading along arbitrary direction
      
         ! calculate displacement from macroscopic deformation gradient
         du = matmul(macroF_tau - Identity, nodePos) 
         
         u(1) = du(dir)
         
      elseif(loadType==18) then
      
         ! transform to natural coordinates
         ksi(1) = 2.d0*nodePos(1)/domainLen0(1)
         ksi(2) = 2.d0*nodePos(2)/domainLen0(2)
         ksi(3) = 2.d0*nodePos(3)/domainLen0(3)
         ! evaluate shape functions at node position
         call evalShapeFuncHEX8(shapeFuncHEX8,ksi,ksiCorner_init)
         ! interpolate displacement at the node (nodePos) 
         ! from corner node displacements
         du = matmul(shapeFuncHEX8,transpose(uCorner_tau))
         
         u(1) = du(dir)
         
      endif

      return
      end SUBROUTINE
      
      subroutine evalShapeFuncHEX8(shapeFuncHEX8,ksi,ksi_nodes)
      implicit none
      real(8), intent(out) :: shapeFuncHEX8(8)
      real(8), intent(in)  :: ksi(3),ksi_nodes(3,8)
!      real(8), parameter :: ksi_nodes(3,8) = &
!         reshape((/ -1.d0, -1.d0, -1.d0, &
!                    +1.d0, -1.d0, -1.d0, &
!                    -1.d0, +1.d0, -1.d0, &
!                    +1.d0, +1.d0, -1.d0, &
!                    -1.d0, -1.d0, +1.d0, &
!                    +1.d0, -1.d0, +1.d0, &
!                    -1.d0, +1.d0, +1.d0, &
!                    +1.d0, +1.d0, +1.d0  /), (/3,8/))
                    
      integer :: iNode
      
      do iNode = 1,8
         shapeFuncHEX8(iNode) = &
            1.d0/8.d0*(1.d0 + ksi(1)*ksi_nodes(1,iNode)) &
                     *(1.d0 + ksi(2)*ksi_nodes(2,iNode)) &
                     *(1.d0 + ksi(3)*ksi_nodes(3,iNode)) 
      enddo
      end subroutine
!------------------------------------------------------------------
      subroutine dload(jelem,pr,time)
      implicit none
      include 'PARDIS.H'   ! need this just because MAX_HISTORY_N
      
      integer, intent(in) :: jelem
      real(8), intent(out) :: pr
      real(8), intent(in) :: time(2)
      
      real(8) :: PI
      parameter(PI=3.14159265359D0)

      real(8) :: tperiod,tramp,t_nodwell,t_dwell,samp,smin, & 
                 uamp,umin
      real(8) :: sload,st_rate
      integer :: icyc
      
      common/load_dwell/tperiod,tramp,t_nodwell,t_dwell,samp,smin, & 
                 uamp,umin
      common/load_mono/sload,st_rate
      common/icyc_flag/icyc

      ! deniz
      integer :: loadType      ! 1: creep, 2: CSR, 3: stress cyclic, 4: displ cyclic. 5/6: load/displ. history
      integer(4) :: padding
      real(8) :: P_cyclic_max  ! in MPa
      real(8) :: P_cyclic_min  ! in MPa
      real(8) :: P_cyclic_period ! in seconds
      integer :: P_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: P_dwell_ramptime
      integer :: P_history_repeat
      real(8) :: P_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: P_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: P_history_N    ! number of load-time points

      real(8) :: U_cyclic_max  !
      real(8) :: U_cyclic_min  !
      real(8) :: U_cyclic_period ! in seconds
      integer :: U_dwell_load    
      real(8) :: U_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: U_history_repeat    ! 1: repeat provided u-time history, 0: keep u at the last provided value
      real(8) :: U_history_values(MAX_HISTORY_N) ! provide u-time points. 
      real(8) :: U_history_times(MAX_HISTORY_N)  ! program will interpolate at integration time steps
      integer :: U_history_N    ! number of U-time points
      
      integer :: load_axis,iStressLateral,iStressMain
      real(8) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6) 
      logical :: specifiedLoadFrame
      logical :: strainControlledSimulation
      real(8) :: biaxial_StressRatioAngle      
      real(8) :: biaxial_StrainRate
      real(8) :: deformMultiplierTensionCompression
      real(8) :: macroStrain_t(6),macroStrain_tau(6)
      
      integer :: strain_history_repeat    ! 1: repeat provided strain-time history, 0: keep u at the last provided value
      real(8) :: strain_history_values(6,MAX_HISTORY_N) ! provide strain-time points. 
      real(8) :: strain_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: strain_history_N    ! number of strain-time points
      
      real(8) :: ksiCorner_init(3,8)  ! initial positions of corner nodes
      real(8) :: uCorner_history_values(3,8,MAX_HISTORY_N) ! provide displacement-time histories of corners of the cubic domain
      real(8) :: uCorner_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: uCorner_history_N           ! number of displacement-time points
      integer :: uCorner_remove_gradient     ! remove gradient-modes (trapezoids) from the deformation

      integer :: mark_for_export(MAX_HISTORY_N) ! mark time step for data-export. 0:default behavior, 1: force simple export, 2: force export all.
      integer :: export_special                 ! mark the special time points of cyclic loadings for data export

      common/load_cond/loadType,P_history_N,U_history_N,strain_history_N,uCorner_history_N, & 
                       P_history_repeat,U_history_repeat,P_dwell_load,U_dwell_load, &
                       P_dwell_ramptime,U_dwell_ramptime,P_cyclic_max, & 
                       P_cyclic_min,P_cyclic_period,P_history_values, & 
                       P_history_times,U_cyclic_max,U_cyclic_min, & 
                       U_cyclic_period,U_history_values,U_history_times, &
                       macroStrain_t, macroStrain_tau, load_axis,strainControlledSimulation, &
                       biaxial_StressRatioAngle,biaxial_StrainRate,deformMultiplierTensionCompression, &
                       iStressLateral,iStressMain,specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv, &
                       strain_history_values,strain_history_times,strain_history_repeat, &
                       uCorner_history_values,uCorner_history_times,uCorner_remove_gradient, &
                       ksiCorner_init,mark_for_export,export_special
                       
      real(8) :: domainRange0(3,2), domainLen0(3)
      real(8) :: domainRange(3,2), domainLen(3)
      COMMON/DOMAINGEO/domainRange0,domainLen0,domainRange,domainLen
      
      real(8) :: P_cyclic_amp,P_cyclic_mean,t_cyc,t2_cyc,tramp_cyc,smax
      real(8) :: w1,w2
      logical :: found
      integer :: iTime, i_cyc

      if(loadType==1) then        ! --------- CREEP LOADING
                  
         if(time(1).le.tramp)then

            pr=sload*time(1)/tramp

         else

            pr=sload

         endif
         
      elseif(loadType==2) then    ! --------- CSR loading
                                    ! not applicable
!         write(*,*)                 
!     &   'dload() called for CSR loading'
  
         pr=0.0

      elseif(loadType==3) then    ! --------- traction controlled cyclic loading
      
         if(P_dwell_load.EQ.1) then ! dwell loading.

            t_cyc = DMOD(time(1),P_cyclic_period)
            
            if(t_cyc.LT.P_dwell_ramptime) then                       ! on ramp-up
               pr = P_cyclic_max*t_cyc/P_dwell_ramptime
            elseif(t_cyc.LT.P_cyclic_period-P_dwell_ramptime) then   ! on P_max
               pr = P_cyclic_max
            elseif(t_cyc.LT.P_cyclic_period) then                    ! on ramp-down
               pr = P_cyclic_max*(P_cyclic_period-t_cyc)/P_dwell_ramptime
            else
               pr = 0.D0
            endif
         
         else
            P_cyclic_amp = 0.5D0*(P_cyclic_max-P_cyclic_min)
            P_cyclic_mean = 0.5D0*(P_cyclic_max+P_cyclic_min)
            tramp_cyc = P_cyclic_mean /  & 
                        (P_cyclic_amp*2*PI/P_cyclic_period)
     
            if (time(1).LT.tramp_cyc) then
               pr = time(1)/tramp_cyc*P_cyclic_mean
            else
               t_cyc = time(1) - tramp_cyc
               pr = P_cyclic_mean +  & 
                        P_cyclic_amp*SIN(2*PI*t_cyc/P_cyclic_period)
            endif
         endif
      elseif(loadType==4) then    ! --------- displacement-controlled cyclic loading
                                    ! not applicable
!            write(*,*)                 
!     &   'dload() called for displacement controlled cyclic loading'

         pr=0.0         

      elseif(loadType==5) then    ! --------- traction-time provided
         t_cyc = time(1)
         if(P_history_repeat.EQ.1) then
            t_cyc = DMOD(time(1),P_history_times(P_history_N))
         endif
         
         found=.FALSE.
         do iTime=1,P_history_N
            if(t_cyc.LT.P_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               pr=P_history_times(1)
            else
               w1 = (P_history_times(iTime)-t_cyc) / & 
            (P_history_times(iTime)-P_history_times(iTime-1))
               w2 = (t_cyc-P_history_times(iTime-1)) / & 
            (P_history_times(iTime)-P_history_times(iTime-1))
               pr = w1*P_history_values(iTime-1) & 
                    + w2*P_history_values(iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            pr=P_history_values(P_history_N)
         endif   
         
      elseif(loadType==6) then    ! --------- disp-time provided
                                    ! not applicable
!            write(*,*) 'dload() called with disp-time history loading'
         
         pr = 0.0
         
      endif

      pr=-pr

      
      return
      end SUBROUTINE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

      SUBROUTINE extractRxnMatrix(IBC,N_UPDATE,IN_ST, & 
                                  IROWPTR,IBINDX,NNZ_LOC,VAL, & 
                                  ICOLPTR_rxn,IBINDX_rxn,VAL_rxn, & 
                                  N_UPDATE_rxn,ISIZE_rxn)
      IMPLICIT none
      INCLUDE 'PARDIS.H'
      
      integer, intent(in) :: IBC(MDOFX),N_UPDATE,IN_ST,NNZ_LOC
      integer, intent(in) :: IBINDX(NNZ_LOC),IROWPTR(N_UPDATE+1)
      real(8), intent(in) :: VAL(NNZ_LOC)
      real(8), intent(out):: VAL_rxn(ISIZE_rxn)
      integer, intent(in) :: ICOLPTR_rxn(N_UPDATE_rxn+1)
      integer, intent(in) :: IBINDX_rxn(ISIZE_rxn)
      integer, intent(in) :: N_UPDATE_rxn, ISIZE_rxn
      
      integer :: NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS

      integer :: k,k_rxn,k_sta,k_end,k_rxn_sta,k_rxn_end
      integer :: I,j_col,i_row,iDOF
      k_rxn = 1
      j_col = 1
      
      VAL_rxn(:) = 0.D0
      
      do I=1,N_UPDATE
         i_row=IN_ST+I           ! EQNs handled by this CPU
         k_sta=IROWPTR(I)        ! data index of the first non-zero coefficient in this row
         k_end=IROWPTR(I+1)-1    ! data index of the last non-zero coefficient in this row

         IF(IBC(i_row).EQ.1)THEN   ! a DISP BC was defined for this EQN (i_row), which is handled by this CPU.

               ! copy this row of stiff. matrix into a the column-ordered K_rxn
               ! K_rxn will be used to determine the 'settlement' (displ. BC) forces on the structure
               k_rxn_sta = ICOLPTR_rxn(j_col)
               k_rxn_end = ICOLPTR_rxn(j_col+1)-1
               
               k = k_sta
               do k_rxn = k_rxn_sta, k_rxn_end
                  iDOF = IBINDX_rxn(k_rxn)    ! DOF/EQN number of this coefficient
                  if (IBC(iDOF).EQ.0) then   
                     VAL_rxn(k_rxn) = VAL(k)  ! a free EQN. copy the coefficient as it is.
                  elseif (iDOF.EQ.i_row) then
                     VAL_rxn(k_rxn) = - VAL(k)! a diagonal, constrained DOF x constrained DOF.
                  else                       
                     VAL_rxn(k_rxn) = 0.D0    ! a non-diagonal, constrained DOF x constrained DOF.
                  endif
                  k = k + 1
               enddo

               j_col = j_col + 1
               
               k_rxn = k_rxn + (k_end-k_sta+1)
                                                
         ENDIF
      enddo
      
      END SUBROUTINE

      ! this subroutine modifies VAL_rxn to give settlement forces due to imposed displs.
      ! dont need this anymore because VAL_rxn is now modified after being extracted in: extractRxnMatrix()
      ! calculates the 'settlement' forces due to non-homogeneous BCs --> F_rxn
      SUBROUTINE calc_F_rxn(VAL_rxn,F_rxn,uBC_proc,EQ_uBC_proc,EQ_uBC, &
                            ICOLPTR_rxn,IBINDX_rxn,IBC,     &
                            ISIZE_rxn,N_UPDATE_rxn,NEQ,NEQ_uBC)
      implicit none
      INCLUDE 'PARDIS.H'
      
      real(8), intent(inout):: VAL_rxn(ISIZE_rxn)
      real(8), intent(inout):: F_rxn(NEQ)
      real(8), intent(in)   :: uBC_proc(N_UPDATE_rxn)
      integer, intent(in)   :: EQ_uBC_proc(N_UPDATE_rxn)
      integer, intent(in)   :: EQ_uBC(NEQ_uBC)
      integer, intent(in)   :: ICOLPTR_rxn(N_UPDATE_rxn+1)
      integer, intent(in)   :: IBINDX_rxn(ISIZE_rxn)
      integer, intent(in)   :: IBC(MDOFX)
      integer, intent(in)   :: ISIZE_rxn,N_UPDATE_rxn,NEQ,NEQ_uBC
      
      integer :: iDOF,jDOF_uBC,j_col,i_diag_uBC,k,k_sta,k_end
      
      do j_col=1,N_UPDATE_rxn
      
         jDOF_uBC = EQ_uBC_proc(j_col)   ! get the EQN # of this DOF (the diagonal element)

         k_sta = ICOLPTR_rxn(j_col)
         k_end = ICOLPTR_rxn(j_col+1)-1 ! COLPTR has the proper index at ICOLPTR_rxn(n_cols+1)

         do k = k_sta, k_end
            iDOF = IBINDX_rxn(k)

            if (iDOF.EQ.jDOF_uBC) then ! diagonal, constrained DOF
               VAL_rxn(k) = - VAL_rxn(k)
            elseif (IBC(iDOF).EQ.1) then  ! non-diagonal, constrained DOF
               VAL_rxn(k) = 0.D0
            endif
         enddo
      enddo

      ! calculate the 'settlement' forces due to imposed displacements
      CALL MATMUL_ColSparse(VAL_rxn,uBC_proc,F_rxn,   &
           IBINDX_rxn,ICOLPTR_rxn,NEQ,N_UPDATE_rxn,ISIZE_rxn)
           
      F_rxn = - F_rxn

      END SUBROUTINE
      
      ! do this part at the beginning of each q-Newton iteration:
      SUBROUTINE extractRxnForce(RESID,F_rxn_proc,IBC,N_UPDATE,IN_ST,NEQ)
      IMPLICIT none
      INCLUDE 'PARDIS.H'

      real(8), intent(in) :: RESID(N_UPDATE)
      real(8), intent(out):: F_rxn_proc(NEQ)
      integer, intent(in) :: N_UPDATE,IN_ST,IBC(MDOFX),NEQ
      
      integer :: I,i_row
      
      F_rxn_proc(:) = 0.D0
      
      do I=1,N_UPDATE
         i_row=IN_ST+I       ! EQNs handled by this CPU
         ! if a DISP BC is defined for this EQN,  save the constraint/reaction force
         if (IBC(i_row).EQ.1) F_rxn_proc(i_row)=RESID(I)
      enddo

      END SUBROUTINE
!***********************************************************************

   ! deniz modification
   !    instead of inserting 1 on the diagonal, the original coeff. is retained,
   !    while setting rest of the row to zero. See modified REMOVE_DISPDOF.
      SUBROUTINE REMOVE_DISPDOF(IBC,N_UPDATE,IN_ST, & 
                                IROWPTR,IBINDX,VAL, & 
                                RESID,NNZ_LOC,IFLAG)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION IBC(MDOFX)
      integer :: NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS

      DIMENSION IBINDX(NNZ_LOC),VAL(NNZ_LOC),RESID(N_UPDATE), & 
                  IROWPTR(N_UPDATE+1)

      DO I=1,N_UPDATE
         IROWNO=IN_ST+I    ! EQNs handled by this CPU

         IS=IROWPTR(I)        ! data index of the first non-zero coefficient in this row
         IE=IROWPTR(I+1)-1    ! data index of the last non-zero coefficient in this row

         IF(IBC(IROWNO).EQ.1)THEN   ! a DISP BC was defined for this EQN (IROWNO), which is handled by this CPU.

            ! do this part only if re-assembling K:
            IF(IFLAG.EQ.0)THEN
               DO J=IS,IE
                  IF(IBINDX(J).NE.IROWNO)THEN   ! set all coeffs of this equation to 0, to impose the displ. B.C.
                     VAL(J)=0.D0                 
                  else                          ! set the diag. coeff to 1. this forces u to come out as 0,
!                    VAL(J)=1.D0                ! since we also set corresp. F to 0 (see below).
                     VAL(J)=VAL(J)! ! ! ! ! ! ! ! deniz - modification ! RHS is also VAL(J). See extractRxnMatrix
                  ENDIF
               ENDDO
            ENDIF

            ! do this part at each q-Newton iteration:
            RESID(I)=0.D0                       ! set corresp. F to 0. this serves two purposes:
                                                ! 1) F is out balance forces, and is used to test convergence. forces at disp BCs are automatically balanced by the support.
                                                ! 2) F is also used to solve for dU in iterations. K-modification above, and RHS_i=0 forces constrained DOFS dU_i to come out 0.
                                                !                 or, dU_i = dU_i*, if non-zero displacements were prescribed. do this only at the first iteration to get a dU that satisfies non-hom. displ. BCs
                                                
            GOTO 140
         ENDIF

         IF(IFLAG.EQ.0)THEN
            DO J=IS,IE
               ICOLM=IBINDX(J)
               IF(IBC(ICOLM).EQ.1)THEN          ! also set columns of SK to zero, to keep the symmetry. this is 
                  VAL(J)=0.D0
               ENDIF
            ENDDO
         ENDIF

 140  CONTINUE
      ENDDO

      END SUBROUTINE
      
! -----------------------------------------------------

      ! added by Deniz
      ! imposes periodic BCs by strictly coupling the DOFs on the opposite sides of a cubic domain.
      ! in general, a node on one side of the domain, may not coincide with another node in unstructured meshes.
      ! in this case, the node is projected to the opposite face, and the triangualar face the node coincides with, is identified.
      ! then the node is coupled with the three nodes of this triangular face, with weights determined from interp. functions.
      ! this couples the node, to the projected point on the face, assuming the face stays linear/planar.
      SUBROUTINE IMPOSE_PERIODIC_BCs & 
                        (coupledDOFs,couplingWeights,nCoupledDOFs, & 
                        IBC,N_UPDATE,IN_ST,IROWPTR,IBINDX,VAL, & 
                        RESID,NNZ_LOC,IFLAG)
      ! EXAMPLE
      ! coupledDOFs  = (u0 u1 u2 u3  v0 v1 v2 v3 ...) here, u0 u1 u2 etc. are EQN numbers of DOFs
      ! couplingWeights=( w1 w2 w3    w1 w2 w3 ...) : u0 is coupled to u1, u2, and u3 with weights w1, w2, w3. (w1+w2+w3=1)
      !
      ! coupledDOFs  = (u0 u1 0 0  v0 v1 v2 v3 ...)
      ! couplingWeights=( 1 0 0    w1 w2 w3 ...) : u0 is directly coupled to u1.
      implicit none
      INCLUDE 'PARDIS.H'
      integer, intent(in):: coupledDOFs(nCoupledDOFs*4)
      real(8), intent(in):: couplingWeights(nCoupledDOFs*3)
      integer, intent(in):: nCoupledDOFs
      integer, intent(in):: N_UPDATE,IN_ST,NNZ_LOC,IFLAG
      integer, intent(in):: IBC(MDOFX),IROWPTR(N_UPDATE+1), &
                            IBINDX(NNZ_LOC)
      real(8), intent(in):: VAL(NNZ_LOC),RESID(N_UPDATE)
      
      integer :: NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDF,NGAUSS
      
      integer :: IROWNO, IS, IE, I, J, iDOF, jDOF, pos1D
      real(8) :: k_inf, fourWeights(4)
      
      k_inf = 1.D5
      
      do I=1,nCoupledDOFs
         do J=1,3
            iDOF = coupledDOFS((I-1)*4+1)
            jDOF = coupledDOFS((I-1)*4+1+J)
            ! iDOF is coupled to jDOF.

            if ((iDOF.GE.IN_ST+1).AND.  & 
                (iDOF.LE.IN_ST+N_UPDATE)) then
               pos1D = IROWPTR(iDOF-IN_ST)
               ! VAL IS COMPACT -- no zeroes.
               ! new couplings will overwrite some zeroes. consider re-shaping val beforehand.
            else
               ! think about this..
               ! CYCLE !(J)
               CONTINUE
            endif
         enddo
      enddo
            
      END SUBROUTINE

!*************************END*******************************************
!******PERFORMS DIVISION OF DOF TO EACH PROCESSOR
      SUBROUTINE NODE_PARTITION(N_UPDATE,IN_ST,IPROC,NPROCS,NX,NDF)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'

      ! nodes are equally distributed to processors
      ! if not divisible, IB number of processors with highest id numbers, are assigned one less node.
      ! example, 14 nodes / 4 processors: 4 4 3 3. IB = 2
      !          15 nodes / 4 processors: 4 4 4 3. IB = 1


      I1=NX/NPROCS+1
      I2=NX/NPROCS

      IA=NX-I2*NPROCS
      IB=NPROCS-IA      ! IB: number of 'empty slots' the last processor has

      I1=I1*NDF
      I2=I2*NDF

      IF((IPROC+1).LE.IA)THEN
         IN_ST=IPROC*I1
         N_UPDATE=I1

      ELSE
         IN_ST=IA*I1+(IPROC-IA)*I2
         N_UPDATE=I2
      ENDIF

      RETURN
      END SUBROUTINE
!***********************END*********************************************


      SUBROUTINE MAKEGRAD_TET4(XB_ST,XBTR_ST,XG_ST,SHPTAU_ST,DVOLTAU_ST, & 
                   NNODE,MCRD,NDOFEL,NGAUSS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XB_ST(6,NDOFEL,NGAUSS),DVOLTAU_ST(NGAUSS), & 
            SHPTAU_ST(MCRD,NNODE,NGAUSS), & 
           NSHR(3),XB_AV(NDOFEL),XBTR_ST(NDOFEL,NGAUSS), & 
           XG_ST(9,NDOFEL,NGAUSS),NSHR0(3)


      NSHR(1)=2
      NSHR(2)=3
      NSHR(3)=3

      NSHR0(1)=1
      NSHR0(2)=1
      NSHR0(3)=2

      IF(MCRD.EQ.3)THEN
         NSTR=6
      ENDIF
      IF(MCRD.EQ.2)THEN
         NSTR=3
      ENDIF
      XB_AV(1:NDOFEL)=0.D0

      DO IPT=1,NGAUSS
         DO I=1,NSTR
            DO J=1,NDOFEL
          XB_ST(I,J,IPT)=0.D0
          JTEMP1=(J-1)/MCRD
          JTEMP2=J-JTEMP1*MCRD
          IF(I.EQ.JTEMP2.AND.I.LE.MCRD) THEN
          XB_ST(I,J,IPT)=SHPTAU_ST(I,JTEMP1+1,IPT)
          ENDIF
          IF(I.GT.MCRD)THEN
             IS=I-MCRD
             ISHR=NSHR(IS)
             IS=NSHR0(IS)
             IF(JTEMP2.EQ.IS)XB_ST(I,J,IPT)=SHPTAU_ST(ISHR,JTEMP1+1,IPT)
             IF(JTEMP2.EQ.ISHR)XB_ST(I,J,IPT)=SHPTAU_ST(IS,JTEMP1+1,IPT)
          ENDIF
            ENDDO
         ENDDO


         DO J=1,NDOFEL
            XBTR_ST(J,IPT)=0.D0
            DO I=1,MCRD
          XBTR_ST(J,IPT)=XBTR_ST(J,IPT)+XB_ST(I,J,IPT)
            ENDDO
            DO I=1,MCRD
          XB_ST(I,J,IPT)=XB_ST(I,J,IPT)-1.D0/3.D0*XBTR_ST(J,IPT)
            ENDDO
            XB_AV(J)=XB_AV(J)+XBTR_ST(J,IPT)*DVOLTAU_ST(IPT)
         ENDDO
         DO J=1,NDOFEL
            JTEMP=(J-1)/MCRD
            JTEMP1=JTEMP+1
            JTEMP2=J-JTEMP*MCRD
            DO I=1,MCRD**2
          ITEMP=(I-1)/MCRD
          ITEMP1=ITEMP+1
          ITEMP2=I-ITEMP*MCRD
          XG_ST(I,J,IPT)=0.D0
         IF(JTEMP2.EQ.ITEMP1)XG_ST(I,J,IPT)=SHPTAU_ST(ITEMP2,JTEMP1,IPT)
          IF(ITEMP2.EQ.1)THEN
          XG_ST(I,J,IPT)=XG_ST(I,J,IPT)-1.D0/3.D0*XBTR_ST(J,IPT)
         ENDIF
            ENDDO
         ENDDO
      ENDDO
      TV=SUM(DVOLTAU_ST(1:NGAUSS))
      XB_AV(1:NDOFEL)=XB_AV(1:NDOFEL)/TV
      DO IPT=1,NGAUSS
         DO J=1,NDOFEL
            DO I=1,MCRD
          XB_ST(I,J,IPT)=XB_ST(I,J,IPT)+1.D0/3.D0*XB_AV(J)
          ITEMP=(I-1)*MCRD+1
          XG_ST(ITEMP,J,IPT)=XG_ST(ITEMP,J,IPT)+1.D0/3.D0*XB_AV(J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE


      SUBROUTINE CALC_DFGRD_TET4(IPT,U,MCRD,NNODE,SHP_ST,DFGRD,NGAUSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(MCRD,NNODE),DFGRD(3,3),SHP_ST(MCRD,NNODE,NGAUSS)

      DO I=1,MCRD
         DO J=1,MCRD
            DFGRD(I,J)=0.D0
            IF(I.EQ.J)DFGRD(I,J)=1.D0
            DO IN=1,NNODE
          DFGRD(I,J)=DFGRD(I,J)+SHP_ST(J,IN,IPT)*U(I,IN)
            ENDDO
         ENDDO
      ENDDO
      IF(MCRD.EQ.2)THEN
         DFGRD(1:3,3)=0.D0
         DFGRD(3,1:3)=0.D0
         DFGRD(3,3)=1.D0
      ENDIF
      RETURN
      END SUBROUTINE



      ! returns the gradient of TET shape function wrt spatial coordinates
      ! input: 
      !     XEL(3,NODES) : element nodal coordinates
      !     IPT          : 
      !  output:
      !     dSHP_dx          : shape function gradient d(N_i)/d(s_j)
      !     DVOL         : volume of the element
      !     XC           :
      SUBROUTINE SHAPE_TET4_new(XEL,IPT,dSHP_dx,DVOL,XC)
      IMPLICIT NONE

      REAL(8):: dSHP_ds(3,4),dSHP_dx(3,4),XJAC(3,3),XJACINV(3,3)
      REAL(8):: XEL(3,4),XC(3)
      REAL(8):: X(4),Y(4),Z(4),edge1(3),edge2(3),edge3(3)
      REAL(8):: DVOL
      INTEGER:: I,J, INODE, IPT

      DO I=1,3
      DO J=1,4
      dSHP_ds(I,J)=0.D0
      dSHP_dx(I,J)=0.D0
      ENDDO
      ENDDO
      
      
      ! dSHP_ds is the gradient of shape functions wrt natural coordinates:
      ! dSHP_ds_ij = d(N_j)/d(s_i)
      ! where N_j are shape functions, and s_i are natural coordinates
      dSHP_ds(1,1)=1.D0

      dSHP_ds(2,2)=1.D0

      dSHP_ds(3,3)=1.D0


      dSHP_ds(1,4)=-1.D0

      dSHP_ds(2,4)=-1.D0

      dSHP_ds(3,4)=-1.D0


   !		WRITE(*,*) 'XEL   IS:  ', XEL

      DO I=1,4
         X(I)=XEL(1,I)
         Y(I)=XEL(2,I)
         Z(I)=XEL(3,I)
      ENDDO

      DO INODE=1,3
         edge1(:)=XEL(:,INODE)-XEL(:,4)
         edge2(:)=XEL(:,INODE)-XEL(:,4)
         edge3(:)=XEL(:,INODE)-XEL(:,4)
      ENDDO

      ! jacobian corresponding to the natural to spatial mapping of the element
      ! XJAC_ij = d(x_i)/d(s_j)
      XJAC(:,1)=edge1(:)
      XJAC(:,2)=edge2(:)
      XJAC(:,3)=edge3(:)

      CALL MATINV3_UEL(XJAC,XJACINV,DVOL)
      IF(DVOL.LT.1.D-16)THEN
         WRITE(*,*)'NEGATIVE VOLUME'
         STOP
      ENDIF


      DVOL=1.D0/6.D0*DABS(DVOL)

      DO INODE=1,4
         dSHP_dx(:,INODE)=MATMUL(dSHP_ds(:,INODE),XJACINV)
      ENDDO

      RETURN
      END SUBROUTINE
      
      SUBROUTINE SHAPE_TET4(XEL,IPT,SHP,DVOL,XC)
      IMPLICIT NONE

      REAL(8):: SHP1(3,4),SHP(3,4),XJAC(3,3),XJACINV(3,3)
      REAL(8):: XEL(3,4),XC(3)
      REAL(8):: X(4),Y(4),Z(4),XX(3),YY(3),ZZ(3)
      REAL(8):: DVOL
      INTEGER:: I,J, INODE, IPT

      DO I=1,3
      DO J=1,4
      SHP1(I,J)=0.D0
      SHP(I,J)=0.D0
      ENDDO
      ENDDO

      SHP1(1,1)=1.D0

      SHP1(2,2)=1.D0

      SHP1(3,3)=1.D0


      SHP1(1,4)=-1.D0

      SHP1(2,4)=-1.D0

      SHP1(3,4)=-1.D0


   !		WRITE(*,*) 'XEL   IS:  ', XEL

      DO I=1,4
         X(I)=XEL(1,I)
         Y(I)=XEL(2,I)
         Z(I)=XEL(3,I)
      ENDDO

      DO I=1,3
         XX(I)=XEL(1,I)-XEL(1,4)
         YY(I)=XEL(2,I)-XEL(2,4)
         ZZ(I)=XEL(3,I)-XEL(3,4)
      ENDDO

      XJAC(1:3,1)=XX(1:3)
      XJAC(1:3,2)=YY(1:3)
      XJAC(1:3,3)=ZZ(1:3)



      CALL MATINV3_UEL(XJAC,XJACINV,DVOL)
      IF(DVOL.LT.1.D-16)THEN
         WRITE(*,*)'NEGATIVE VOLUME'
         STOP
      ENDIF


      DVOL=1.D0/6.D0*DABS(DVOL)

      DO INODE=1,4
         SHP(:,INODE)=MATMUL(XJACINV,SHP1(:,INODE))
      ENDDO



      RETURN
      END SUBROUTINE

      ! calculates the TET interpolation of a field Phi(s) at the natural coordinates s,
      ! given the nodal values
      SUBROUTINE interpolateTETs(sPoint,interpValue,nodalValuesElem, & 
                                 nQuantDOF)
      implicit none
      real(8), intent(in) :: sPoint(3)
      real(8), intent(out):: interpValue(nQuantDOF)
      real(8), intent(in) :: nodalValuesElem(nQuantDOF,4)
      integer, intent(in) :: nQuantDOF
      
      real(8) :: N(4)
      
      N(1) = sPoint(1)
      N(2) = sPoint(2)
      N(3) = sPoint(3)
      N(4) = 1-sPoint(1)-sPoint(2)-sPoint(3)
      
      interpValue = MATMUL(nodalValuesElem,N)
      
      END SUBROUTINE
      
      ! calculates the TET interpolation of a field Phi(s) at the spatial coordinates x,
      ! given the nodal values Phi_n and element nodal positions x_n
      SUBROUTINE interpolateTETx(sPoint,interpValue,nodalValuesElem, & 
                                 nQuantDOF)
      implicit none
      real(8), intent(in) :: sPoint(3)
      real(8), intent(out):: interpValue(nQuantDOF)
      real(8), intent(in) :: nodalValuesElem(nQuantDOF,4)
      integer, intent(in) :: nQuantDOF
      
      ! NOT YET CODED
      interpValue = 0.d0
      
      END SUBROUTINE

      SUBROUTINE CALC_DFG_TET4(UE,NNODE,NDOFEL,MCRD,DFGRD_ST,DVOL_ST, & 
       	XJAC_ST,SHP_ST,XJBAR,IBAR,NGAUSS)


      IMPLICIT NONE

      REAL(8)::  XJBAR,TV,DVOL,FAC,XJAC
      INTEGER::  MCRD, NNODE,NGAUSS, IBAR, NDOFEL,IPT
      REAL(8)::  UE(MCRD,NNODE), & 
         XJAC_ST(NGAUSS),DFGRD_ST(3,3,NGAUSS),DFGRD_INV(3,3), & 
         DVOL_ST(NGAUSS),SHP_ST(MCRD,NNODE,NGAUSS),DFGRD(3,3)


      XJBAR=0.0D0

      TV=0.0D0


      DO IPT=1,NGAUSS
         CALL CALC_DFGRD_TET4(IPT,UE,MCRD,NNODE,SHP_ST,DFGRD,NGAUSS)
         CALL MATINV3_UEL(DFGRD,DFGRD_INV,XJAC)
         XJAC_ST(IPT)=XJAC
         DVOL=DVOL_ST(IPT)
         XJBAR=XJBAR+XJAC*DVOL
         TV=TV+DVOL
         DFGRD_ST(1:3,1:3,IPT)=DFGRD(1:3,1:3)
      ENDDO
      XJBAR=XJBAR/TV
      IF(IBAR.EQ.1)THEN
         DO IPT=1,NGAUSS
            FAC=XJBAR**(1.D0/3.D0)*XJAC_ST(IPT)**(-1.D0/3.D0)
            DFGRD_ST(1:3,1:3,IPT)=DFGRD_ST(1:3,1:3,IPT)*FAC
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE

      SUBROUTINE MAKESPIN(W,V)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(3,3),V(3)
      W=0.D0
      W(1,2)=-V(3)
      W(1,3)=V(2)
      W(2,1)=V(3)
      W(2,3)=-V(1)
      W(3,1)=-V(2)
      W(3,2)=V(1)
      RETURN
      END SUBROUTINE
!------------------------------------------------
      SUBROUTINE MATINV3_UEL(A,AI,DET)
      IMPLICIT NONE
      REAL(8), intent(in) :: A(3,3)
      REAL(8), intent(out):: AI(3,3)
      REAL(8), intent(out):: DET

      DET=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2) & 
            *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1) & 
           *A(1,3)*A(2,2))

      AI(1,1) =  ( A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AI(1,2) = -( A(1,2)*A(3,3)-A(1,3)*A(3,2))/DET
      AI(1,3) =  ( A(1,2)*A(2,3)-A(1,3)*A(2,2))/DET
      AI(2,1) = -( A(2,1)*A(3,3)-A(2,3)*A(3,1))/DET
      AI(2,2) =  ( A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AI(2,3) = -( A(1,1)*A(2,3)-A(1,3)*A(2,1))/DET
      AI(3,1) =  ( A(2,1)*A(3,2)-A(2,2)*A(3,1))/DET
      AI(3,2) = -( A(1,1)*A(3,2)-A(1,2)*A(3,1))/DET
      AI(3,3) =  ( A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
      RETURN
      END SUBROUTINE
!---------------------------------------------------------------
      SUBROUTINE DOT3(DOTT,U,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3),V(3),W(3)
      W(1)=U(1)*V(1)
      W(2)=U(2)*V(2)
      W(3)=U(3)*V(3)
      DOTT=SUM(W)
      RETURN
      END SUBROUTINE
!----------------------------------------------
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
      
!-----------------------------------------------
      SUBROUTINE RUDCMP_UEL(F,RD,UD)
      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER*8 (I-N)
!      REAL *8 LI1,LI2,LI3,LAMDA1,LAMDA2,LAMDA3
!
      DIMENSION RD(3,3),F(3,3),UD(3,3)

!
!	WRITE(6,*)'ENTERING RU'
      O3=1.0D0/3.0D0
      ROOT3=DSQRT(3.D0)
      C11=F(1,1)*F(1,1)+F(2,1)*F(2,1)+F(3,1)*F(3,1)
      C12=F(1,1)*F(1,2)+F(2,1)*F(2,2)+F(3,1)*F(3,2)
      C13=F(1,1)*F(1,3)+F(2,1)*F(2,3)+F(3,1)*F(3,3)
      C23=F(1,2)*F(1,3)+F(2,2)*F(2,3)+F(3,2)*F(3,3)
      C22=F(1,2)*F(1,2)+F(2,2)*F(2,2)+F(3,2)*F(3,2)
      C33=F(1,3)*F(1,3)+F(2,3)*F(2,3)+F(3,3)*F(3,3)
      C1212=C12*C12
      C1313=C13*C13
      C2323=C23*C23
      C2313=C23*C13
      C1223=C12*C23
      C1213=C12*C13
      S11=C22*C33-C2323
      UI1=O3*(C11+C22+C33)
      UI2=S11+C11*C22+C33*C11-C1212-C1313
      UI3=C11*S11+C12*(C2313-C12*C33)+C13*(C1223-C22*C13)
      UI1S=UI1*UI1
      Q    =DSQRT(-DMIN1(O3*UI2-UI1S,0.D0))
      R    =0.5D0*(UI3-UI1*UI2)+UI1*UI1S
      XMOD =Q*Q*Q
      SCL1 =.5D0+DSIGN(.5D0,XMOD-1.D-30)
      SCL2 =.5D0+DSIGN(.5D0,XMOD-DABS(R))
      SCL0 =DMIN1(SCL1,SCL2)
      SCL1 =1.D0-SCL0
      SDETM=DACOS(R/(XMOD+SCL1))*O3
      Q  =SCL0*Q
      CT3=Q*DCOS(SDETM)
      ST3=Q*ROOT3*DSIN(SDETM)
      SDETM=SCL1*DSQRT(DMAX1(0.0D0,R))
      AA=2.D0*(CT3+SDETM)+UI1
      BB=-CT3+ST3-SDETM+UI1
      CC=-CT3-ST3-SDETM+UI1
      XLAMDA1=DSQRT(DMAX1(AA,0.D0))
      XLAMDA2=DSQRT(DMAX1(BB,0.D0))
      XLAMDA3=DSQRT(DMAX1(CC,0.D0))
      SDETM=XLAMDA1*XLAMDA2
      XLI1=XLAMDA1+XLAMDA2+XLAMDA3
      XLI2= SDETM+XLAMDA2*XLAMDA3+XLAMDA3*XLAMDA1
      XLI3= SDETM*XLAMDA3/XLI1
      S11= C11+XLI3
      S22= C22+XLI3
      S33= C33+XLI3
      S12= C2313-C12*S33
      S13= C1223-S22*C13
      S23=-C2323+S22*S33
      SDETM=1.D0/(XLI1*(S11*S23+C12*S12+C13*S13))
      C11=C11+XLI2
      C22=C22+XLI2
      C33=C33+XLI2
      SI11=SDETM*S23
      SI12=SDETM*S12
      SI13=SDETM*S13
      SI22=SDETM*( S11*S33-C1313)
      SI23=SDETM*(-S11*C23+C1213)
      SI33=SDETM*( S11*S22-C1212)
      S12=C12*SI12
      S13=C13*SI13
      S23=C23*SI23
      UI11=C11*SI11+S12+S13
      UI22=S12+C22*SI22+S23
      UI33=S13+S23+C33*SI33
      UI12=C11*SI12+C12*SI22+C13*SI23
      UI13=C11*SI13+C12*SI23+C13*SI33
      UI23=C12*SI13+C22*SI23+C23*SI33
      RD(1,1)=F(1,1)*UI11+F(1,2)*UI12+F(1,3)*UI13
      RD(1,2)=F(1,1)*UI12+F(1,2)*UI22+F(1,3)*UI23
      RD(1,3)=F(1,1)*UI13+F(1,2)*UI23+F(1,3)*UI33
      RD(2,1)=F(2,1)*UI11+F(2,2)*UI12+F(2,3)*UI13
      RD(2,2)=F(2,1)*UI12+F(2,2)*UI22+F(2,3)*UI23
      RD(2,3)=F(2,1)*UI13+F(2,2)*UI23+F(2,3)*UI33
      RD(3,1)=F(3,1)*UI11+F(3,2)*UI12+F(3,3)*UI13
      RD(3,2)=F(3,1)*UI12+F(3,2)*UI22+F(3,3)*UI23
      RD(3,3)=F(3,1)*UI13+F(3,2)*UI23+F(3,3)*UI33

      DO I=1,3
         DO J=1,3
            UD(I,J)=0.D0
            DO K=1,3
               UD(I,J)=UD(I,J)+RD(K,I)*F(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
!--------------------------------------------------
      SUBROUTINE DROTMAT(R,A_2D)
         IMPLICIT REAL*8(A-H,O-Z)
         DIMENSION R(3,3),A_2D(6,6),A_4D(3,3,3,3)
         DO I=1,3
            DO J=1,3
               DO K=1,3
             DO L=1,3
                A_4D(I,J,K,L)=R(K,I)*R(L,J)
             ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO I=1,3
            A_2D(I,1)=A_4D(I,I,1,1)
            A_2D(I,2)=A_4D(I,I,2,2)
            A_2D(I,3)=A_4D(I,I,3,3)
            A_2D(I,4)=A_4D(I,I,1,2)+A_4D(I,I,2,1)
            A_2D(I,5)=A_4D(I,I,1,3)+A_4D(I,I,3,1)
            A_2D(I,6)=A_4D(I,I,2,3)+A_4D(I,I,3,2)
         ENDDO
         
         DO I=1,3
            A_2D(4,I)=A_4D(1,2,I,I)
            A_2D(5,I)=A_4D(1,3,I,I)
            A_2D(6,I)=A_4D(2,3,I,I)
         ENDDO
         
         A_2D(4,4)=A_4D(1,2,1,2)+A_4D(1,2,2,1)
         A_2D(4,5)=A_4D(1,2,1,3)+A_4D(1,2,3,1)
         A_2D(4,6)=A_4D(1,2,2,3)+A_4D(1,2,3,2)
         
         A_2D(5,4)=A_4D(1,3,1,2)+A_4D(1,3,2,1)
         A_2D(5,5)=A_4D(1,3,1,3)+A_4D(1,3,3,1)
         A_2D(5,6)=A_4D(1,3,2,3)+A_4D(1,3,3,2)
         
         A_2D(6,4)=A_4D(2,3,1,2)+A_4D(2,3,2,1)
         A_2D(6,5)=A_4D(2,3,1,3)+A_4D(2,3,3,1)
         A_2D(6,6)=A_4D(2,3,2,3)+A_4D(2,3,3,2)

         A_2D(1:6,4:6)=A_2D(1:6,4:6)/2.D0

         A_2D(4:6,1:6)=2.D0*A_2D(4:6,1:6)
            
         RETURN
      END SUBROUTINE
      
      SUBROUTINE solveAXb(A,N,b,err)
      integer, intent(in):: N
      real(8), intent(in):: A(N,N)
      real(8), intent(inout):: b(N)
      integer, intent(out):: err
      ! locals
      real(8) :: A_local(N,N) ! DPOSV overwrites the matrix !!
      integer*4 :: err_4
      
      A_local = A
      CALL DPOSV('U',N,1,A_local,N,b,N,err_4)
      err = err_4
      END SUBROUTINE
      
      SUBROUTINE eigen_Lapack(valEigen,vecEigen,A,N,INFO)
      implicit none
      integer, intent(in) :: N
      real(8), intent(in) :: A(N,N)
      real(8),intent(out) :: valEigen(N),vecEigen(N,N)
      integer, intent(out):: INFO
      ! locals
      character, parameter :: UPLO = 'U'
      character, parameter :: Operation = 'V'
      real(8), allocatable :: AP(:)
      integer, allocatable :: IWORK(:)
      real(8), allocatable :: WORK(:)
      integer(4):: INFO_4
      integer :: LDA, LWORK, LIWORK
      INTEGER :: sizeUpper
      integer :: i,j,iPos

      sizeUpper = N*(N+1)/2
      LWORK=N*N+6*N+1
      LIWORK=5*N+3
      
      allocate(AP(sizeUpper))
      allocate(IWORK(LIWORK))
      allocate(WORK(LWORK))
      
      ! pack to upper diagonal column-wise into linear array
      AP = 0.d0
      iPos = 0
      do j=1,N
         do i=1,j
            iPos = iPos + 1
            AP(iPos) = A(i,j)
         enddo
      enddo
      
      ! see: http://www.netlib.org/lapack/lapack-3.1.1/html/dspevd.f.html
      CALL DSPEVD(Operation,UPLO,N,AP, &
                  valEigen,vecEigen,   &
                  N,WORK,LWORK,IWORK,LIWORK,INFO_4)
                  
      INFO = INFO_4
                  
      END SUBROUTINE

      
      ! solve symmetric, positive semi-definite systems
      SUBROUTINE solveAXb_LDL(A,N,b,iError)
      integer, intent(in):: N
      real(8), intent(in):: A(N,N)
      real(8), intent(inout):: b(N)
      integer, intent(out):: iError
      ! locals
      real(8) :: A_local(N,N)
      integer :: ipiv(N)
      integer :: lwork
      real(8) :: lworkOptimal
      real(8), allocatable :: work(:)
      integer(4) :: err_4
      
      A_local = A ! lapack overwrites A
      
      ! find optimal work buffer size
      lwork = -1
      call dsysv('L',N,1,A_local,N,ipiv,b,N,lworkOptimal,lwork,err_4)
      lwork = nint(lworkOptimal)
      if (err_4/=0) then
         write(*,*) 'dsysv query error',err_4
      endif
      allocate(work(lwork))
      
      ! solve symmetric system using LDL**T decomposition
      call dsysv('L',N,1,A_local,N,ipiv,b,N,work,lwork,err_4)
      
      if (err_4/=0) then
         write(*,*) 'dsysv error',err_4
      endif
      
      iError = err_4
      END SUBROUTINE
      
      ! LU solver for non-symmetric, positive semi-definite matrices
      SUBROUTINE solveAXb_LU(A,N,b,iError)
      integer, intent(in):: N
      real(8), intent(in):: A(N,N)
      real(8), intent(inout):: b(N)
      integer, intent(out):: iError
      ! locals
      real(8) :: A_local(N,N)
      integer :: ipiv(N)
      integer(4) :: err_4
      
      A_local = A ! lapack overwrites A
      
      ! compute LU factorization
      call dgetrf(N,N,A_local,N,ipiv,err_4)
      if (err_4/=0) then
         write(*,*) 'dgetrf query error',err_4
      endif
      
      !CALL printMatrix(A(:,:),6,6, & 
      !'KKT matrix:',0)

      ! solve non-symmetric system using LDL**T decomposition
      call dgetrs('N',N,1,A_local,N,ipiv,b,N,err_4)      
      if (err_4/=0) then
         write(*,*) 'dsysv error',err_4
      endif
      
      iError = err_4
      END SUBROUTINE
      
      
      SUBROUTINE getGaussValues(gaussValues,SVARS_all, & 
                                 posQuantity,nQuantDOF, & 
                                  NSTATE,NELX,NGAUSS)
      implicit none
      INCLUDE 'PARDIS.H'
      real(8), intent(in):: SVARS_all(NELX*NGAUSS*NSTATE)
      integer, intent(in):: NSTATE
      integer, intent(in):: NELX
      integer, intent(in):: NGAUSS
      integer, intent(in):: posQuantity, nQuantDOF
      real(8), intent(out):: gaussValues(nQuantDOF,NELX*NGAUSS)
      !locals
      integer iGauss,iQuantDOF
      do iGauss=1,NELX*NGAUSS
         gaussValues(:,iGauss) =  & 
            SVARS_all((iGauss-1)*NSTATE+posQuantity : &
                      (iGauss-1)*NSTATE+posQuantity+nQuantDOF-1)
      enddo
      END SUBROUTINE

      SUBROUTINE getGaussValuesProc(gaussValues,SVARS, & 
                                 posQuantity,nQuantDOF, & 
                                 NSTATE,NELX,NGAUSS, &
                                 NELST,NELNOX)
      implicit none
      INCLUDE 'PARDIS.H'
      real(8), intent(in):: SVARS(NELNOX*NGAUSS*NSTATE)
      integer, intent(in):: NSTATE
      integer, intent(in):: NELX,NELST,NELNOX
      integer, intent(in):: NGAUSS
      integer, intent(in):: posQuantity, nQuantDOF
      real(8), intent(out):: gaussValues(nQuantDOF,NELX*NGAUSS)
      !locals
      integer iGauss,iGaussProc,iQuantDOF
      gaussValues = 0.d0
      do iGaussProc=1,NELNOX*NGAUSS
         iGauss = (NELST-1)*NSTATE + iGaussProc
         gaussValues(:,iGauss) =  & 
            SVARS((iGaussProc-1)*NSTATE+posQuantity : &
                      (iGaussProc-1)*NSTATE+posQuantity+nQuantDOF-1)
      enddo
      END SUBROUTINE
      
      function getNodeIDinCube(i,j,k,Lcube) result(nodeID)
         integer, intent(in) :: i,j,k,Lcube
         integer:: nodeID
         nodeID = (k-1)*(Lcube+1)*(Lcube+1)+(j-1)*(Lcube+1)+i
      end function
      
      function getBlockIDinCube(i,j,k,Lcube) result(blockID)
         integer, intent(in) :: i,j,k,Lcube
         integer:: blockID
         blockID = (k-1)*(Lcube)*(Lcube)+(j-1)*(Lcube)+i
      end function
      
      function factorial(n) result(n_fact)
         integer, intent(in) :: n
         integer:: n_fact
         n_fact = 1
         do i=1,n
            n_fact=n_fact*i
         enddo
      end function factorial


      
      SUBROUTINE lineSearch(STEP,G0,G,F,F_ext, &
                            stressMacroscopic,ddsddeMacroscopic, & 
                            defGradMacroscopic,dpMacroscopic, &
                            volumeRVE,volume0RVE, &
                            VT,VTAU,DVTAU,incVTAU_init,G0XYZ, & 
                            NEQ,WLF_c2,STEP_max,G_history,G_n,method,info, & 
                            XUND,IJK,IMID,NPROPSG,PROMAT2D,DELTME, & 
                            NBC1,IBC1,KINC,NDF,N_UPDATE,IN_ST,IROWPTR, & 
                            IBINDX,VAL,RESID,NNZ_LOC, & 
                 VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                            NGAUSS,PNEWDT,IPROC)
      use fBarPatch
      implicit none
      include 'PARDIS.H'
      real(8), intent(inout):: STEP
      real(8), intent(in)   :: G0
      real(8), intent(inout):: G
      real(8), intent(inout):: F(NEQ),F_ext(NEQ)
      real(8), intent(out)  :: stressMacroscopic(6),ddsddeMacroscopic(6,6)
      real(8), intent(out)  :: defGradMacroscopic(3,3),dpMacroscopic(6)
      real(8), intent(out)  :: volumeRVE,volume0RVE
      real(8), intent(in)   :: VT(NEQ)
      real(8), intent(inout):: VTAU(NEQ)
      real(8), intent(in)   :: DVTAU(NEQ)
      real(8), intent(in)   :: incVTAU_init(NEQ)
      real(8), intent(in)   :: G0XYZ(MAXCRD)
      integer, intent(in)   :: NEQ
      real(8), intent(in)   :: WLF_c2, STEP_max
      real(8), intent(out)  :: G_history(2,100)
      integer, intent(out)  :: G_n
      integer, intent(in)   :: method
      integer, intent(out)  :: info
      integer, intent(in)   :: IPROC

      REAL(8), intent(in)   :: XUND(MAXCRD)
      REAL(8), intent(in)   :: PROMAT2D(MAXPROPS,MAXGRP)
      integer, intent(in)   :: N_UPDATE
      REAL(8), intent(inout):: RESID(N_UPDATE),VAL(NNZ_LOC)
      real(8), intent(inout)::VAL_rxn(ISIZE_rxn)
      integer, intent(in)   ::IBINDX_rxn(ISIZE_rxn), & 
                              ICOLPTR_rxn(N_UPDATE_rxn)
      integer, intent(in)   ::ISIZE_rxn,N_UPDATE_rxn
      REAL(8), intent(in)   :: DELTME,PNEWDT

      ! error triggers
      logical :: listLocalUmatErrorsOccurred(7)
      logical :: trigger_umatError
      logical :: trigger_BroyDidNotConv
      logical :: trigger_cutBack
      logical :: trigger_reformK
      COMMON/triggers/trigger_umatError,listLocalUmatErrorsOccurred, &
                      trigger_BroyDidNotConv,trigger_cutBack,trigger_reformK

      INTEGER::IMID(MAXDIM*MAXEL),IJK(MNELX),NBC1(MBOUND1)
      INTEGER::NPROPSG(MAXGRP),IBC1(MDOFX),KINC,NDF,IN_ST,NNZ_LOC
      INTEGER::NGAUSS,IFLAG,IBINDX(NNZ_LOC),IROWPTR(N_UPDATE+1)      

      logical :: calcStiffnessMatrix, success
      real(8) :: STEP_a, STEP_b, STEP_c, Ga, Gb, Gc
      real(8) :: qB, qC
      info = 0
      
      IFLAG = 1   ! during line search LoadVC should never re-calculate the stiffness matrix
      
      ! an initial unit step u0 is already taken, and 
      ! F is the residual, and G is the directional derivative, u0.F, at (u + u0)
      STEP_c = 0.0
      Gc = G0
      STEP_b = 0.0
      Gb = G0
      STEP_a = STEP
      Ga = G
      
      G_n = 2
      G_history = 0.0
      G_history(1,1) = G0
      G_history(2,1) = STEP_b
      G_history(1,2) = G
      G_history(2,2) = STEP_a
      
      if (IPROC.EQ.0) write(*,*) 'G0, G:', G0, G

      do while (abs(Ga).GT.abs(WLF_c2*G0) & ! continue if: too steep (Wolfe-2 not satisfied) 
           .AND.Ga*Gb.GT.0.D0             & ! and also on the same side of the minimum 
           .AND.Ga.LE.Gb                  & ! and also the curvature is retained 
           .AND.STEP_a.LT.STEP_max)         ! and also STEP has not exceeded the max value

         ! shift and expand [b a]
         STEP_c = STEP_b
         STEP_b = STEP_a
         STEP_a = STEP_a * 2
         Gc = Gb
         Gb = Ga
         
         ! calculate residual F, and directional derivative G, at u_a
         VTAU = incVTAU_init + STEP_a*DVTAU
         
         calcStiffnessMatrix = (IFLAG == 0)
         call calcFBarStabilizers(VT,VTAU, &! displacements over cycle, time increment info
            IJK,G0XYZ,                     &                     ! mesh,parallelism etc.
            calcStiffnessMatrix,success)

         ! deniz - added output variable F_ext
         CALL LOADVC(F,F_ext, &
                     stressMacroscopic,ddsddeMacroscopic, & 
                     defGradMacroscopic,dpMacroscopic, &
                     volumeRVE,volume0RVE, &
                     VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                     PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE, & 
                     IN_ST,IROWPTR,IBINDX,VAL,RESID,NNZ_LOC, & 
                     VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                     IFLAG,PNEWDT)
         if(trigger_umatError) then    ! UMAT error occured on one of the processors. R is corrupted.
            info = -1                  ! signal to fem_uel
            return
         endif
            
         Ga = DOT_PRODUCT(DVTAU(1:NEQ),F(1:NEQ))
         G_n = G_n + 1
         if(G_n.LE.100) then
            G_history(1,G_n) = Ga
            G_history(2,G_n) = STEP_a
         endif

         if (IPROC.EQ.0) write(*,*) 'LSCH:', STEP_a, Ga
         
      enddo
      
      STEP = STEP_a
      G = Ga
      ! now STEP, VTAU, F, G: has the most recent values
      ! we either have: 
      ! (*) Wolfe-2 satisfied :
      !     can return VTAU, STEP, F now
      if (abs(G).LT.abs(WLF_c2*G0)) then
         info = 1
         return
      elseif (Ga.GT.Gb) then 
            ! an unexpected situation.. isn't the potential positive definite?
            ! for some reason, G will never be zero along this direction.
            ! get the lowest G possible by a quadratic fit..
            
            qB = -((Gc-Gb)*(STEP_a-STEP_c)**2+ & 
                  (Ga-Gc)*(STEP_b-STEP_c)**2)/ & 
                  ((STEP_a-STEP_b)*(STEP_a-STEP_c)*(STEP_b-STEP_c))
            qC = -((Gb-Gc)*(STEP_a-STEP_c)+ & 
                  (Gc-Ga)*(STEP_b-STEP_c))/ & 
                  ((STEP_a-STEP_c)*(STEP_b-STEP_c)*(STEP_a-STEP_b))

            STEP = STEP_c - qB/(2*qC)
            
            VTAU = incVTAU_init + STEP*DVTAU
         
            calcStiffnessMatrix = (IFLAG == 0)
            call calcFBarStabilizers(VT,VTAU, &! displacements over cycle, time increment info
            IJK,G0XYZ,                     &                     ! mesh,parallelism etc.
            calcStiffnessMatrix,success)

            ! deniz - added output variable F_ext         
            CALL LOADVC(F,F_ext, &
                        stressMacroscopic,ddsddeMacroscopic, & 
                        defGradMacroscopic,dpMacroscopic, &
                        volumeRVE,volume0RVE, &
                        VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                        PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE, & 
                        IN_ST,IROWPTR,IBINDX,VAL,RESID,NNZ_LOC, & 
                  VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                        IFLAG,PNEWDT)
            if(trigger_umatError) then    ! UMAT error. R is corrupted.
               info = -1                  ! signal to fem_uel
               return
            endif
            G = DOT_PRODUCT(DVTAU(1:NEQ),F(1:NEQ))
            G_n = G_n + 1
            if(G_n.LE.100) then
               G_history(1,G_n) = G
               G_history(2,G_n) = STEP
            endif
            if (IPROC.EQ.0) write(*,*) 'LSCH invalid curv ! ! :', STEP, G
            info = 2
            return
         
      ! OR (*) Ga landed on the other side of the minimum, and on a non-acceptable region
      !        do a zoomed interval search in [b a]
      elseif (Ga*Gb.LE.0.D0) then
      
            ! # # # WHAT if G_a > G_b ? revert to G_a and return ??? # # # 
            ! # # # WHAT if G_a > G_b ? revert to G_a and return ??? # # # 
         
         do while (abs(G).GT.abs(WLF_c2*G0))
            ! assume total potential has quadratic shape 
            !     => G is linear along search direction. minimize G.
            STEP = STEP_b + (STEP_a-STEP_b)/(Gb - Ga)*Gb
            VTAU = incVTAU_init + STEP*DVTAU
            
            calcStiffnessMatrix = (IFLAG == 0)
            call calcFBarStabilizers(VT,VTAU, &! displacements over cycle, time increment info
            IJK,G0XYZ,                     &                     ! mesh,parallelism etc.
            calcStiffnessMatrix,success)

            ! deniz - added output variable F_ext
            CALL LOADVC(F,F_ext,&
                        stressMacroscopic,ddsddeMacroscopic, & 
                        defGradMacroscopic,dpMacroscopic, &
                        volumeRVE,volume0RVE, &
                        VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                        PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE, & 
                        IN_ST,IROWPTR,IBINDX,VAL,RESID,NNZ_LOC, & 
                        VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                        IFLAG,PNEWDT)
            G = DOT_PRODUCT(DVTAU(1:NEQ),F(1:NEQ))
            
            if (IPROC.EQ.0) write(*,*) 'zoom:', STEP_b, Gb, STEP_a, Ga
            if (IPROC.EQ.0) write(*,*) '-->',STEP, G
            
            G_n = G_n + 1
            if(G_n.LE.100) then
               G_history(1,G_n) = G
               G_history(2,G_n) = STEP
            endif
            
            if (G*Gb.GT.0.D0) then  ! G landed between b and the minimum
               Gb = G
               STEP_b = STEP
            else                    ! G landed between the minimum and a
               Ga = G
               STEP_a = STEP                  
            endif
         enddo
         
         ! # # # maybe put an iteration counter up here? # # # 
         ! # # # maybe put an iteration counter up here? # # # 
         
         ! we have reached the Wolfe-2 condition. accept
         info=3
         return
                     
      ! OR (*) we have exceeded STEP_max (we're on a non-acceptable region and before the minimum)
      else
         info = 4
         return
      endif
      END SUBROUTINE
      
      SUBROUTINE lineSearchSTRANG(STEP,G0,G,F,F_ext, &
                        stressMacroscopic,ddsddeMacroscopic, & 
                        defGradMacroscopic,dpMacroscopic, &
                        volumeRVE,volume0RVE, &
                        VT,VTAU,DVTAU,G0XYZ, & 
              incVTAU_init,NEQ,WLF_c2,STEP_max,G_history,G_n,method,info, & 
                        XUND,IJK,IMID,NPROPSG,PROMAT2D,DELTME, & 
                        NBC1,IBC1,KINC,NDF,N_UPDATE,IN_ST,IROWPTR, & 
                        IBINDX,VAL,RESID,NNZ_LOC, & 
                VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                        NGAUSS,PNEWDT,IPROC)
      use fBarPatch
      implicit none
      include 'PARDIS.H'
      real(8), intent(inout):: STEP
      real(8), intent(in)   :: G0
      real(8), intent(inout):: G
      real(8), intent(inout):: F(NEQ),F_ext(NEQ)
      real(8), intent(out)  :: stressMacroscopic(6),ddsddeMacroscopic(6,6)
      real(8), intent(out)  :: defGradMacroscopic(3,3),dpMacroscopic(6)
      real(8), intent(out)  :: volumeRVE,volume0RVE
      real(8), intent(in)   :: VT(NEQ)
      real(8), intent(inout):: VTAU(NEQ)
      real(8), intent(in)   :: DVTAU(NEQ)
      real(8), intent(in)   :: incVTAU_init(NEQ)
      real(8), intent(in)   :: G0XYZ(MAXCRD)
      integer, intent(in)   :: NEQ
      real(8), intent(in)   :: WLF_c2, STEP_max
      real(8), intent(out)  :: G_history(2,100)
      integer, intent(out)  :: G_n
      integer, intent(in)   :: method
      integer, intent(out)  :: info
      integer, intent(in)   :: IPROC

      REAL(8), intent(in)   :: XUND(MAXCRD)
      REAL(8), intent(in)   :: PROMAT2D(MAXPROPS,MAXGRP)
      integer, intent(in):: N_UPDATE
      REAL(8), intent(inout):: RESID(N_UPDATE),VAL(NNZ_LOC)
      real(8), intent(inout)::VAL_rxn(ISIZE_rxn)
      integer, intent(in)   ::IBINDX_rxn(ISIZE_rxn), & 
                              ICOLPTR_rxn(N_UPDATE_rxn)
      integer, intent(in)   ::ISIZE_rxn,N_UPDATE_rxn
      REAL(8), intent(in)   :: DELTME,PNEWDT
      
      ! error triggers
      logical :: listLocalUmatErrorsOccurred(7)
      logical :: trigger_umatError
      logical :: trigger_BroyDidNotConv
      logical :: trigger_cutBack
      logical :: trigger_reformK
      COMMON/triggers/trigger_umatError,listLocalUmatErrorsOccurred, &
                      trigger_BroyDidNotConv,trigger_cutBack,trigger_reformK

      INTEGER::IMID(MAXDIM*MAXEL),IJK(MNELX),NBC1(MBOUND1)
      INTEGER::NPROPSG(MAXGRP),IBC1(MDOFX),KINC,NDF,IN_ST,NNZ_LOC
      INTEGER::NGAUSS,IFLAG,IBINDX(NNZ_LOC),IROWPTR(N_UPDATE+1)     

      logical :: calcStiffnessMatrix, success      

      real(8) :: STEP_a, STEP_b, STEP_c, Ga, Gb, Gc
      info = 0
      
      IFLAG = 1 ! never recalculate K during line search
      
      ! an initial unit step u0 is already taken, and 
      ! F is the residual, and G is the directional derivative, u0.F, at (u + u0)
      STEP_b = 0.0
      Gb = G0
      STEP_a = STEP
      Ga = G
      
      G_n = 2
      G_history = 0.0
      G_history(1,1) = G0
      G_history(2,1) = STEP_b
      G_history(1,2) = G
      G_history(2,2) = STEP_a
      
      do while (Ga*Gb.GT.0.0) ! continue until A and B are on different sides of the minimum
      
         ! shift and expand [b a]
         STEP_c = STEP_b
         STEP_b = STEP_a
         STEP_a = STEP_a * 2
         Gc = Gb
         Gb = Ga
         
         ! calculate residual F, and directional derivative G, at u_a
         VTAU = incVTAU_init + STEP_a*DVTAU

         calcStiffnessMatrix = (IFLAG == 0)
         call calcFBarStabilizers(VT,VTAU,   &  ! displacements over cycle, time increment info
                                 IJK,G0XYZ,  &  ! mesh,parallelism etc.
                                 calcStiffnessMatrix,success)

         ! deniz - added output variable F_ext
         CALL LOADVC(F,F_ext, &
                     stressMacroscopic,ddsddeMacroscopic, & 
                     defGradMacroscopic,dpMacroscopic, &
                     volumeRVE,volume0RVE, &
                     VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                     PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE, & 
                     IN_ST,IROWPTR,IBINDX,VAL,RESID,NNZ_LOC, & 
               VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                     IFLAG,PNEWDT)
         if(trigger_umatError) then    ! UMAT error. R is corrupted.
            info = -1                  ! signal to fem_uel
            return
         endif
         Ga = DOT_PRODUCT(DVTAU(1:NEQ),F(1:NEQ))
         G_n = G_n + 1
         if(G_n.LE.100) then
            G_history(1,G_n) = Ga
            G_history(2,G_n) = STEP_a
         endif
         
      enddo
      
      STEP = STEP_a
      G = Ga
      ! now STEP, VTAU, F, G: has the most recent values
      ! Illinois algorithm to find zero:
      ! convergence requires both (1) SLOPE and (2) STEP WIDTH to diminish below some tolerance
      !  or simply, for A and B to fall on the same side of the minimum.. (why?)
      do while ((Ga*Gb.LT.0.0).AND. & 
                 (abs(G).GT.abs(WLF_c2*G0).OR.      & 
                  abs(STEP_b-STEP_a).GT.WLF_c2*0.5D0*(STEP_a+STEP_b)))
              
         ! assume total potential has quadratic shape 
         !     => G is linear along search direction. minimize G.
         STEP = STEP_b + (STEP_a-STEP_b)/(Gb - Ga)*Gb
         VTAU = incVTAU_init + STEP*DVTAU
         
         calcStiffnessMatrix = (IFLAG == 0)
         call calcFBarStabilizers(VT,VTAU,   &  ! displacements over cycle, time increment info
                                 IJK,G0XYZ,  &  ! mesh,parallelism etc.
                                 calcStiffnessMatrix,success)
         
         ! deniz - added output variable F_ext
         CALL LOADVC(F,F_ext, &
                     stressMacroscopic,ddsddeMacroscopic, & 
                     defGradMacroscopic,dpMacroscopic, &
                     volumeRVE,volume0RVE, &
                     VTAU,VT,XUND,IJK,IMID,NPROPSG, & 
                     PROMAT2D,DELTME,NBC1,IBC1,KINC,N_UPDATE, & 
                     IN_ST,IROWPTR,IBINDX,VAL,RESID,NNZ_LOC, & 
                     VAL_rxn,ICOLPTR_rxn,IBINDX_rxn,N_UPDATE_rxn,ISIZE_rxn, & 
                     IFLAG,PNEWDT)
         if(trigger_umatError) then    ! UMAT error. R is corrupted.
            info = -1                  ! signal to fem_uel
            return
         endif
         G = DOT_PRODUCT(DVTAU(1:NEQ),F(1:NEQ))

         G_n = G_n + 1
         if(G_n.LE.100) then
            G_history(1,G_n) = G
            G_history(2,G_n) = STEP
         endif
         
         if (G*Ga.GT.0.D0) then  ! G landed between the minimum and a
            Gb = 0.5D0*Gb
         else                    ! G landed between b and the minimum
            STEP_b = STEP_a
            Gb = Ga                
         endif
         
         STEP_a = STEP
         Ga = G
      enddo
         
      END SUBROUTINE
      
      ! determines whether to execute a process, that is only executed periodically.
      ! execution periods can be specified in terms of time periods, or increment periods
      ! if (current time %MODULO% period) = 0 --> execute (do_it=.TRUE.)
      ! set period_N to 0, to disable increment trigger
      ! set period_T to 0.0, to disable time period trigger
      ! if both active, period_N is used.
      ! also provide real(8) :: last_T. this stores the last time the execution is made
      ! do NOT modify last_T outside the subroutine               
      SUBROUTINE doAtPeriod(TIME,STEP,period_T,last_T,period_N,do_it)
      implicit none
      real(8), intent(in) :: TIME,period_T
      real(8), intent(inout):: last_T
      integer, intent(in) :: STEP, period_N
      logical, intent(out):: do_it
      integer :: last_T_N, this_T_N
      logical :: lastTime
      
      do_it = .FALSE.
      
      ! if both period_N and period_T is given, decide according to period_N
      if(period_N.NE.0)then
         if(MOD(STEP,period_N).EQ.0) do_it = .TRUE.
      elseif(period_T.NE.0.0) then
         if (TIME.EQ.last_T) then
            do_it = .TRUE.
         else
            last_T_N = last_T / period_T
            this_T_N = TIME / period_T
            if (this_T_N.GT.last_T_N) then
               do_it = .TRUE.
               last_T = TIME
            endif
         endif
      endif
      
      END SUBROUTINE
      

      SUBROUTINE printIntMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      integer, intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      
      if(LEN_TRIM(matName).GT.0)then
      
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE printVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      real(8), intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(vecName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE
      
      SUBROUTINE printIntVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      integer, intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(vecName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE

      SUBROUTINE printMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      real(8), intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(matName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE printMatrix_Sparse(MAT,IBINDX,IROWPTR, & 
                                    n_rows,n_cols,n_entries, &
                                    indexing,matName,printTo)
      implicit none
      real(8), intent(in) :: MAT(n_entries)
      integer, intent(in) :: IBINDX(n_entries)      
      integer, intent(in) :: IROWPTR(n_rows+1)
      integer, intent(in) :: n_rows, n_cols, n_entries
      integer, intent(in) :: indexing
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer :: i_col,i_col_prev,i_char,i_row,k,k_sta,k_end
      character(len=16) :: formatStr
      
      if(LEN_TRIM(matName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      
      formatStr = '(F10.0,A)'
      
      do i_row = 1,n_rows
      
         k_sta = IROWPTR(i_row)
         k_end = IROWPTR(i_row+1)-1
         
         if(indexing.EQ.2) then
            k_sta = k_sta + 1  ! C-style index, convert to F-index.
            k_end = k_end + 1  
         endif
         
         i_col_prev = 0
         do k = k_sta, k_end

            i_col = IBINDX(k)
            
            if(indexing.EQ.2) i_col = i_col+1   ! C-style index, convert to F-index.
            
            ! print zeroes
            do i_char = i_col_prev+1,i_col-1
               if(printTo.GT.0) then
                  write(printTo,formatStr,advance='no') 0.D0, ' '
               else
                  write(*,formatStr,advance='no') 0.D0, ' '
               endif
            enddo
            
            ! print coefficient
            if(printTo.GT.0) then
               write(printTo,formatStr,advance='no') MAT(k),' '
            else
               write(*,formatStr,advance='no') MAT(k),' '
            endif
            
            i_col_prev = i_col
         enddo


         ! print remaining zeroes
         do i_char = i_col_prev+1,n_cols
            if(printTo.GT.0) then
               write(printTo,formatStr,advance='no') 0.D0, ' '
            else
               write(*,formatStr,advance='no') 0.D0, ' '
            endif
         enddo
         
         ! line feed
         if(printTo.GT.0) then
            write(printTo,*) ' '
         else
            write(*,*) ' '
         endif

            
      enddo      
      
      END SUBROUTINE
      

      SUBROUTINE printMatrix_Sparse_Compact( &
                                    MAT,IBINDX,IROWPTR, & 
                                    n_rows,n_cols,iRowStart,idProc,n_entries, &
                                    indexing,matName,printTo)
      implicit none
      real(8), intent(in) :: MAT(n_entries)
      integer, intent(in) :: IBINDX(n_entries)      
      integer, intent(in) :: IROWPTR(n_rows+1)
      integer, intent(in) :: n_rows, n_cols, iRowStart, idProc, n_entries
      integer, intent(in) :: indexing
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer :: i_col,i_col_prev,i_char,i_row,k,k_sta,k_end
      character(len=16) :: formatStr
      
      if(idProc==0)then
         if(printTo.GT.0) then
             write(printTo,*) n_cols,n_cols,n_entries
         else
             write(*,*) n_cols,n_cols,n_entries
         endif
      endif
      
      do i_row = 1,n_rows
      
         k_sta = IROWPTR(i_row)
         k_end = IROWPTR(i_row+1)-1
         
         if(indexing.EQ.2) then
            k_sta = k_sta + 1  ! C-style index, convert to F-index.
            k_end = k_end + 1  
         endif
         
         i_col_prev = 0
         do k = k_sta, k_end

            i_col = IBINDX(k)
            
            if(indexing.EQ.2) i_col = i_col+1   ! C-style index, convert to F-index.
            
            ! print coefficient
            if(printTo.GT.0) then
               write(printTo,*) iRowStart+i_row-1,i_col-1,MAT(k)
            else
               write(*,*) iRowStart+i_row-1,i_col-1,MAT(k)
            endif
            
            i_col_prev = i_col
         enddo

            
      enddo      
      
      END SUBROUTINE
      
      SUBROUTINE assertInt(int1,int2,msg)
         integer, intent(in) :: int1,int2
         character(len=*), intent(in) :: msg
         if (int1.NE.int2) then
            write(*,*) msg
            STOP
         endif
      END SUBROUTINE
      SUBROUTINE assertBool(bool1,bool2,msg)
         logical, intent(in) :: bool1,bool2
         character(len=*), intent(in) :: msg
         if (bool1.NEQV.bool2) then
            write(*,*) msg
            STOP
         endif
      END SUBROUTINE
      SUBROUTINE assertReal(real1,real2,msg)
         real(8), intent(in) :: real1,real2
         character(len=*), intent(in) :: msg
         if (abs(real1-real2).GT.10*tiny(real2)) then
            write(*,*) msg
            STOP
         endif
      END SUBROUTINE
      
      
      SUBROUTINE printSLUinputs(vecName,NSTEP,M,N,NEQ,VAL,NNZ_LOC,IN_ST, &
                 IROWPTR_LOC,IBINDX_LOC,N_UPDATE_LOC,NRHS, &
                 LUstruct, SOLVEstruct, indexing)

      implicit none
      integer, intent(in):: NSTEP,M,N,NEQ,NNZ_LOC,IN_ST,N_UPDATE_LOC, &
                            NRHS, LUstruct,SOLVEstruct
      integer, intent(in):: IROWPTR_LOC(N_UPDATE_LOC+1)
      integer, intent(in):: IBINDX_LOC(NNZ_LOC)
      real(8), intent(in):: VAL(NNZ_LOC)
      integer, intent(in) :: indexing
      character*(*),intent(in):: vecName

      integer:: i,j
      
      write(542,'(A,I0,A)') '------ STEP ',NSTEP,' ------'
      write(542,*) vecName
      write(542,'(9(A,I0))') 'M:',M,' N:',N,' NEQ:',NEQ,' NNZ_LOC:',NNZ_LOC, &
                           ' IN_ST:',IN_ST,' N_UPD_LOC:',N_UPDATE_LOC, &
                         ' NRHS:',NRHS,' LUstruct:',LUstruct,' SLSTR:',&
                           SOLVEstruct
                           
      CALL printIntMatrix(IROWPTR_LOC(:),N_UPDATE_LOC+1,1, & 
                         'IROWPTR',542)
      CALL printIntMatrix(IBINDX_LOC(:),NNZ_LOC,1, &
                         'IBINDX',542)

      CALL printMatrix_Sparse(VAL,IBINDX_LOC,IROWPTR_LOC,N_UPDATE_LOC,M,&
                           NNZ_LOC,indexing,'VAL',542)

      END SUBROUTINE
      
      SUBROUTINE setSLUoptions(SLUoptions,solver_CalledBefore,solver_LUisUpToDate)
         USE SUPERLU_MOD
         implicit none
         
         integer(SUPERLU_PTR), intent(in) :: SLUoptions
         logical, intent(in) :: solver_CalledBefore, solver_LUisUpToDate
         
         if(.NOT.solver_CalledBefore) then
            CALL set_superlu_options(SLUoptions,FACT=DoFact)   ! SLU has never been called before.
                                                               ! LU from scratch
            call set_superlu_options(SLUoptions,PrintStat=YES)
         else
            if(solver_LUisUpToDate) then        ! K has not changed since last SLU
               CALL set_superlu_options(SLUoptions,FACT=FACTORED)
               call set_superlu_options(SLUoptions,PrintStat=NO)
            else                                ! K was re-constructed since last SLU. Sparsity structure is the same.
               CALL set_superlu_options(SLUoptions,FACT=SamePattern_SameRowPerm)
               call set_superlu_options(SLUoptions,PrintStat=NO)
            endif
         endif
      END SUBROUTINE
      
      SUBROUTINE justCalledSLU(solver_CalledBefore,solver_LUisUpToDate)
      implicit none
      
      logical, intent(inout) :: solver_CalledBefore, solver_LUisUpToDate
      !logical, intent(inout) :: trigger_BroyDidNotConv, trigger_umatError
      
      solver_CalledBefore = .TRUE.
      solver_LUisUpToDate = .TRUE.

      END SUBROUTINE


      ! sets IFLAG and solver_LUisUpToDate 
      ! according to the state of the solver deduced from solver flags and triggers.
      ! IN:
      !     solverMethod:     1=Newton-Raphson
      !                       2=modified Newton (recalculate K0 at each increment)
      !                       3=modified Newton - Keep (keep K throughout increments)
      !                       4=BFGS            (recalculate K0 at each increment)
      !                       5=BFGS - Keep     (keep K throughout increments)
      !                       6=Broyden
      ! OUT:
      !     IFLAG               signals to LOADVC whether to re-calculate the stiffness matrix
      !     solver_LUisUpToDate signals to SLU that the stiffness matrix has been re-assembled
      !
      SUBROUTINE setFLAGS(IFLAG,solver_LUisUpToDate,stage,NSTEP,&
                          solverMethod,solver_CalledBefore, &
                          solver_reformedK_thisInc, &
                          trigger_cutBack,trigger_reformK)
      use options, only: solmeth_newtRaph,solmeth_modNewt,solmeth_modNewtKeep,solmeth_BFGS,solmeth_BFGSkeep,solmeth_Broyden
      implicit none
      include 'solver_stages.inc'

      integer, intent(out)  :: IFLAG
      logical, intent(out)  :: solver_LUisUpToDate
      integer, intent(in)   :: stage
      integer, intent(in)   :: NSTEP
      integer, intent(in)   :: solverMethod
      logical, intent(in)   :: solver_CalledBefore
      logical, intent(out)  :: solver_reformedK_thisInc
      logical, intent(in)   :: trigger_cutBack
      logical, intent(inout):: trigger_reformK

      if(.NOT.solver_CalledBefore) then
         IFLAG = 0
         solver_LUisUpToDate = .FALSE.
         solver_reformedK_thisInc = .TRUE. ! mark that we have reformed K_inv during this increment at least once
         return
      endif
      
      if(trigger_reformK) then     ! forced K-recalculation is triggered when 
         IFLAG = 0                     ! 1) convergence cannot be reached even after time step cut backs
         solver_LUisUpToDate = .FALSE. ! 2) umat error occurs during a K-recalculation
         
         trigger_reformK = .FALSE.     ! reset the trigger. it has done its job.
         solver_reformedK_thisInc = .TRUE. ! mark that we have reformed K_inv during this increment at least once         
         return
      endif         

      ! initialize
      IFLAG = 1
      solver_LUisUpToDate = .TRUE.

      if(solverMethod == solmeth_newtRaph) then                 ! newton-Raphson: always recalculate & invert the stiff. matrix
         IFLAG = 0
      elseif(solverMethod == solmeth_modNewt) then              ! modified Newton
         if(stage == stage_initVTAU) then                      ! calculate K at the beginning of each time increment
            IFLAG = 0
         endif
      elseif(solverMethod == solmeth_modNewtKeep) then          ! modified Newton - Keep
         if((stage == stage_initVTAU).AND.(NSTEP == 1)) then   ! calculate K only at the beginning of the simulation
            IFLAG = 0
         endif
      elseif(solverMethod == solmeth_BFGS) then                 ! BFGS update
         if(stage == stage_initVTAU) then                      ! calculate K at the beginning of each time increment
            IFLAG = 0
         endif
      elseif(solverMethod == solmeth_BFGSkeep) then             ! BFGS update - Keep K
         if((stage == stage_initVTAU).AND.(NSTEP == 1)) then   ! calculate K only at the beginning of the simulation
            IFLAG = 0
         endif
      elseif(solverMethod == solmeth_Broyden) then             ! Broyden update - Keep K
         if((stage == stage_initVTAU).AND.(NSTEP == 1)) then   ! calculate K only at the beginning of the simulation
            IFLAG = 0
         endif
      endif

      if (IFLAG.EQ.0) then                 ! if K is to be recalculated, 
         solver_LUisUpToDate = .FALSE.     ! signal to SuperLU that a refactorization is due.
         solver_reformedK_thisInc = .TRUE. ! mark that we have reformed K_inv during this increment at least once
      endif
                                   
      END SUBROUTINE
      

      SUBROUTINE checkConvergence(solver_Converged,convergenceNorm,  &
                                  eps_F,tol_F,F,F0,F_ext,F_ext_max,  &
                                  eps_u,tol_u,DVTAU,DVTAU0,VT,NEQ,   & 
                                  check_u,                           &
                                  conv_F, conv_u, belowNoise)
      implicit none
      logical, intent(out)  :: solver_Converged
      integer, intent(in)   :: convergenceNorm
      real(8), intent(in)   :: eps_F, eps_u
      real(8), intent(inout):: tol_F, tol_u
      real(8), intent(inout):: F_ext_max
      real(8), intent(in)   :: F(NEQ), F0(NEQ), F_ext(NEQ)
      real(8), intent(in)   :: DVTAU(NEQ), DVTAU0(NEQ), VT(NEQ)
      integer, intent(in)   :: NEQ
      logical, intent(in)   :: check_u  ! check convergence on the displacement too
      real(8), intent(out)  :: conv_F, conv_u
      logical, intent(out)  :: belowNoise
      
      real(8) :: iter_F, iter_u
      
      real(8), parameter :: tol_F_MIN = 1d-6 ! MICRONEWTONS - assuming this is the 'noise' external force for a load-free structure
      belowNoise = .false.
      
      solver_Converged = .FALSE.
      conv_F = 0.D0
      conv_u = 0.D0
      
      if (convergenceNorm.EQ.2) then ! L_2 norm
         !tol_F = max(tol_F,eps_F*DSQRT(DOT_PRODUCT(F_ext(1:NEQ),F_ext(1:NEQ))))
         tol_F = eps_F*DSQRT(DOT_PRODUCT(F_ext(1:NEQ),F_ext(1:NEQ)))
         if(tol_F < tol_F_MIN) then
            tol_F = tol_F_MIN ! below noise
            belowNoise = .true.
         endif
         iter_F = DSQRT(DOT_PRODUCT(F(1:NEQ),F(1:NEQ)))
         conv_F = iter_F / tol_F
         if(check_u) then
            tol_u = max(tol_u,eps_u*DSQRT(DOT_PRODUCT(DVTAU0(1:NEQ),DVTAU0(1:NEQ))))
            iter_u = DSQRT(DOT_PRODUCT(DVTAU(1:NEQ),DVTAU(1:NEQ)))
            conv_u = iter_u / tol_u
         endif
      else  ! convergenceNorm = 1, L_max norm
         !tol_F = max(tol_F,eps_F*MAXVAL(abs(F_ext(1:NEQ))))
         tol_F = eps_F*MAXVAL(abs(F_ext(1:NEQ)))
         if(tol_F < tol_F_MIN) then
            tol_F = tol_F_MIN ! below noise
            belowNoise = .true.
         endif
         iter_F = MAXVAL(abs(F(1:NEQ)))
         conv_F = iter_F / tol_F
         if(check_u) then
            tol_u = max(tol_u,eps_u*MAXVAL(abs(DVTAU0(1:NEQ))))
            iter_u = MAXVAL(abs(DVTAU(1:NEQ)))
            conv_u = iter_u / tol_u
         endif
      endif
      
      if((iter_F.LE.tol_F).AND.  &               ! convergence on the residual force
        ((iter_u.LE.tol_u).OR.(.NOT.check_u))) & ! convergence on the displacement
        solver_Converged = .TRUE.
         
      if(solver_Converged) then
         F_ext_max = MAXVAL(abs(F_ext(1:NEQ)))
      endif
      
      END SUBROUTINE
      
      SUBROUTINE BFGS_addWV(W,V,DVTAU,F,F_prev,G,G0,STEP,NEQ,MAXREF,N_BFGS)
      implicit none
      
      real(8), intent(inout):: W(NEQ,MAXREF+1),V(NEQ,MAXREF+1)
      real(8), intent(in)   :: DVTAU(NEQ)
      real(8), intent(in)   :: F(NEQ),F_prev(NEQ),G,G0,STEP
      integer, intent(in)   :: NEQ,MAXREF
      integer, intent(inout):: N_BFGS
      
      N_BFGS = N_BFGS + 1

      W(1:NEQ,N_BFGS) = - DVTAU(1:NEQ) / (G - G0)
      V(1:NEQ,N_BFGS) = F(1:NEQ) &
                      - F_prev(1:NEQ)* &
                     (1.D0+STEP*DSQRT((G-G0)/(-STEP*G0)))
                   
      END SUBROUTINE

      SUBROUTINE BFGS_pre_F(W,V,F_BFGS,NEQ,MAXREF,N_BFGS)
      implicit none
      
      real(8), intent(in)   :: W(NEQ,MAXREF+1),V(NEQ,MAXREF+1)
      real(8), intent(inout):: F_BFGS(NEQ)
      integer, intent(in)   :: NEQ,MAXREF,N_BFGS
      integer :: i_BFGS

      do i_BFGS = 1,N_BFGS
         F_BFGS(1:NEQ) = F_BFGS(1:NEQ) + &
                     V(1:NEQ,N_BFGS-i_BFGS+1)* &
                     DOT_PRODUCT(W(1:NEQ,N_BFGS-i_BFGS+1),F_BFGS(1:NEQ))
      enddo

      END SUBROUTINE

      SUBROUTINE BFGS_post_U(W,V,DVTAU,NEQ,MAXREF,N_BFGS)
      implicit none
      
      real(8), intent(in)   :: W(NEQ,MAXREF+1),V(NEQ,MAXREF+1)
      real(8), intent(inout):: DVTAU(NEQ)
      integer, intent(in)   :: NEQ,MAXREF,N_BFGS
      integer :: i_BFGS
      
      ! post-multiply the direction (after back-substitution)
      do i_BFGS = 1,N_BFGS
         DVTAU(1:NEQ) = DVTAU(1:NEQ) + &
         W(1:NEQ,i_BFGS)*DOT_PRODUCT(V(1:NEQ,i_BFGS),DVTAU(1:NEQ))
      enddo
                   
      END SUBROUTINE

      
      
      SUBROUTINE printSLUerror(INFO,BERR,NRHS,IPROC)
      implicit none
      integer, intent(in) :: INFO,NRHS,IPROC
      real(8), intent(in) :: BERR(NRHS)
      integer :: i
      if (IPROC.NE.0) return
      if (INFO.EQ.0) then
         !write (*,*) 'Backward error: ', &
         !             (BERR(i), i = 1, NRHS)
      else
         write(*,*) 'F_PDGSSVX error !'
         write(*,*) 'INFO from f_pdgssvx = ', INFO
      endif
      END SUBROUTINE

      SUBROUTINE task_Start_superLU(solver_LUisUpToDate,wtime)
      use task_Timing
      implicit none
      logical, intent(in) :: solver_LUisUpToDate
      real(8), intent(in) :: wtime
      real(8) :: task_time
      if(solver_LUisUpToDate) then
         CALL task_Start(task_LU_backsubs,wtime)
      else
         CALL task_Start(task_LU_refactorize,wtime)
      endif
      END SUBROUTINE task_Start_superLU
      SUBROUTINE task_Record_superLU(solver_LUisUpToDate,wtime,lastRefactTime)
      use task_Timing
      implicit none
      logical, intent(in) :: solver_LUisUpToDate
      real(8), intent(in) :: wtime
      real(8), intent(out) :: lastRefactTime
      real(8) :: task_time
      if(solver_LUisUpToDate) then
         CALL task_Record(task_LU_backsubs,wtime,task_time,-1)
      else
         CALL task_Record(task_LU_refactorize,wtime,task_time,-1)
         lastRefactTime = task_time
      endif
      END SUBROUTINE task_Record_superLU


		! overwrites any boundary conditions read in by READBC
		! and defines displacement BCs for the nodes on the domain boundary
		! for macroscopic-strain-controlled deformation
      SUBROUTINE setBC_StrainControlled(NBC1,VBC1,VBCS1,IDF,DIR_U,NBC2,VBC2,VBCS2,NBC3,VBC3,NDF,IBC,IBCT)
      
      use grains
      
      implicit none
      include 'PARDIS.H'

      REAL(8):: VBC1(MBOUND1),VBC2(MBOUND2),VBC3(MBOUND3)
      REAL(8):: VBCS1(MBOUND1),VBCS2(MBOUND2)
      REAL(8):: PRESSB(MAXEL)
      INTEGER:: NBC3(MBOUND3)
      INTEGER:: NBC1(MBOUND1),NBC2(MBOUND2)
      INTEGER:: IDF(MBOUND1),IBC(MDOFX),IBCT(MDOFX)
      INTEGER:: DIR_U(MBOUND1) 
      INTEGER:: ILF(MAXEL),IFACE_TET(12),IFACE_BRICK(24)
      INTEGER:: IBELEM(MAXEL), NDF
      INTEGER:: NX,NELX,NEQ,NPLANE,NODE,NDFF,NGAUSS
      COMMON/ELEMNT/NX,NELX,NEQ,NPLANE,NODE,NDFF,NGAUSS
      INTEGER:: NBOUD1,NBOUD2,NBOUD3,NEQ_uBC
      COMMON/BOUNDA/NBOUD1,NBOUD2,NBOUD3,NEQ_uBC ! NEQ_uBC: total # of DOFs with displ. BCs
      COMMON/PRESS/PRESSB,IBELEM,ILF,IFACE_TET,IFACE_BRICK
      
      integer :: iNodeIdx, nodeID

      !reset displ. BCs
      NBOUD1=0
      NBC1(:) = 0
      VBC1(:) =0.D0
      VBCS1(:) = 0.0
      NEQ_uBC = 0
      IBC(:) = 0
      ! IDF is reset below
      
      !reset point forces
      NBOUD2=0 
      NBC2(:)=0
      VBC2(:)=0.0D0
      VBCS2(:)=0.0D0
      
      ! reset pressure BCs
      NBOUD3 = 0
      IBCT(:)=0
      IBELEM(:) = 0
      NBC3(:) = 0
      VBC3(:) = 0.D0
      ILF(:)=0
      PRESSB(:)=0.D0

      ! set displacement BC information for strain-controlled deformation
      do iNodeIdx = 1,nDomainBoundaryNodes
         nodeID = domainBoundaryNodes(iNodeIdx)
         
         if (IBC(3*(nodeID-1)+1).EQ.1) then
            write(*,*) 'duplicate entries in domainBoundaryNodes'
            write(*,*) 'nodeID', nodeID
            STOP
         endif
         
         NBOUD1 = NBOUD1 + 1
         
         IF(4*NBOUD1.GT.MBOUND1)THEN  ! for NBC1()
            WRITE(*,*)'INSUFFICIENT MEMORY-INCREASE MBOUND1'
            STOP
         ENDIF

         IDF(NBOUD1) = 1
         DIR_U(NBOUD1) = 1
         NBC1(4*(NBOUD1-1)+1) = nodeID
         NBC1(4*(NBOUD1-1)+2) = 1	! u imposed
         NBC1(4*(NBOUD1-1)+3) = 1	! v imposed
         NBC1(4*(NBOUD1-1)+4) = 1	! w imposed
         
         ! mark the EQNs/DOFs as constrained
         IBC(3*(nodeID-1)+1) = 1 ! u imposed
         IBC(3*(nodeID-1)+2) = 1 ! v imposed
         IBC(3*(nodeID-1)+3) = 1 ! w imposed
         NEQ_uBC = NEQ_uBC + 1

      enddo

      END SUBROUTINE
      
      SUBROUTINE setDELTME(DELTME,PNEWDT,FINTME,XMXDT,XMINDT,TOTME,   &
                           adaptive_time_stepping,suppress_specialPoints,&
                           NITER,NITER_REF,   &
                           prevTimeCutBackForSpecial,DELTME_lastBeforeSpecTime, &
                           EXPORT_next_step)

      implicit none
      
      include 'PARDIS.H'
      real(8), intent(inout) :: DELTME
      real(8), intent(inout) :: PNEWDT
      real(8), intent(in) :: FINTME,XMXDT,XMINDT,TOTME
      logical, intent(in) :: adaptive_time_stepping,suppress_specialPoints ! some options
      logical, intent(inout) :: prevTimeCutBackForSpecial
      real(8), intent(inout) :: DELTME_lastBeforeSpecTime      
      logical, intent(out):: EXPORT_next_step
      integer, intent(in) :: NITER,NITER_REF

      integer :: loadType      ! 1: Creep, 2: CSR, 
                               ! 3: stress controlled cyclic
                               ! 4: displacement controlled cyclic
                               ! 5: load-history provided
                               ! 6: displ.-history provided
                               ! 7: stress-history provided
                               ! 8: strain-history provided
      integer :: dummy1,dummy2
      integer*4 :: padding
      real(8) :: P_cyclic_max  ! in MPa
      real(8) :: P_cyclic_min  ! in MPa
      real(8) :: P_cyclic_period ! in seconds
      integer :: P_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: P_dwell_ramptime
      integer :: P_history_repeat
      real(8) :: P_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: P_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: P_history_N    ! number of load-time points
      
      real(8) :: U_cyclic_max  !
      real(8) :: U_cyclic_min  !
      real(8) :: U_cyclic_period ! in seconds
      integer :: U_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: U_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: U_history_repeat    ! 1: repeat provided u-time history, 0: keep u at the last provided value
      real(8) :: U_history_values(MAX_HISTORY_N) ! provide u-time points. 
      real(8) :: U_history_times(MAX_HISTORY_N)  ! program will interpolate at integration time steps
      integer :: U_history_N    ! number of U-time points
      
      integer :: load_axis,iStressLateral,iStressMain
      real(8) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6) 
      logical :: specifiedLoadFrame
      logical :: strainControlledSimulation
      real(8) :: biaxial_StressRatioAngle      
      real(8) :: biaxial_StrainRate
      real(8) :: deformMultiplierTensionCompression
      real(8) :: macroStrain_t(6),macroStrain_tau(6)
      
      integer :: strain_history_repeat    ! 1: repeat provided strain-time history, 0: keep u at the last provided value
      real(8) :: strain_history_values(6,MAX_HISTORY_N) ! provide strain-time points. 
      real(8) :: strain_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: strain_history_N    ! number of strain-time points
      
      real(8) :: ksiCorner_init(3,8)  ! initial positions of corner nodes
      real(8) :: uCorner_history_values(3,8,MAX_HISTORY_N) ! provide displacement-time histories of corners of the cubic domain
      real(8) :: uCorner_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: uCorner_history_N           ! number of displacement-time points
      integer :: uCorner_remove_gradient     ! remove gradient-modes (trapezoids) from the deformation
      
      integer :: mark_for_export(MAX_HISTORY_N) ! mark time step for data-export. 0:default behavior, 1: force simple export, 2: force export all.
      integer :: export_special                 ! mark the special time points of cyclic loadings for data export

      common/load_cond/loadType,P_history_N,U_history_N,strain_history_N,uCorner_history_N, & 
                       P_history_repeat,U_history_repeat,P_dwell_load,U_dwell_load, &
                       P_dwell_ramptime,U_dwell_ramptime,P_cyclic_max, & 
                       P_cyclic_min,P_cyclic_period,P_history_values, & 
                       P_history_times,U_cyclic_max,U_cyclic_min, & 
                       U_cyclic_period,U_history_values,U_history_times, &
                       macroStrain_t, macroStrain_tau, load_axis,strainControlledSimulation, &
                       biaxial_StressRatioAngle,biaxial_StrainRate,deformMultiplierTensionCompression, &
                       iStressLateral,iStressMain,specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv, &
                       strain_history_values,strain_history_times,strain_history_repeat, &
                       uCorner_history_values,uCorner_history_times,uCorner_remove_gradient, &
                       ksiCorner_init,mark_for_export,export_special
                      
      REAL(8):: SAMP,SMIN,UAMP,UMIN
      REAL(8):: TPERIOD,TRAMP,T_NODWELL,T_DWELL
      COMMON/LOAD_DWELL/TPERIOD,TRAMP,T_NODWELL,T_DWELL,SAMP,SMIN,UAMP,UMIN

      real(8) :: PI
      parameter(PI=3.14159265359D0)
      
      !locals
      real(8) :: PNEWDT1,FINTME_try,XTOL
      real(8) :: tramp_cyc,t_cyc,P_cyclic_amp,P_cyclic_mean,U_cyclic_amp,U_cyclic_mean
      integer :: I
      
      EXPORT_next_step = .FALSE. ! reset the flag for exporting at the next time step
      
      XTOL = 1.0D-8
      
      ! adaptively increase time step, if larger than suggested step by UMAT
      if(adaptive_time_stepping) then
         PNEWDT1 = (1.D0*NITER_REF/(1.D0*NITER))**0.5D0
         ! --- special case --- !
         ! in dwell loading, increase adaptive time steps faster, if N cycles > 50
         if (loadType==3) then
            t_cyc = DMOD(FINTME,P_cyclic_period)
            if (P_dwell_load==1) then
               if(t_cyc > P_dwell_ramptime .and. &                  
                  t_cyc < P_cyclic_period-P_dwell_ramptime) then ! if on stress hold, boost.
                  PNEWDT1 = (4.0D0*NITER_REF/(1.D0*NITER))**0.5D0
               endif
            endif
         endif            
         PNEWDT = MAX(PNEWDT,PNEWDT1)
      endif
      DELTME=PNEWDT*DELTME
      
      if (prevTimeCutBackForSpecial) then ! adopt the last time increment that we got prior to special-time cut-back
         DELTME = MAX(DELTME,DELTME_lastBeforeSpecTime)
      endif
      prevTimeCutBackForSpecial = .false.
      
      ! but do not exceed the max allowed time step
      DELTME=MIN(DELTME,XMXDT)
      
      ! save this time-increment for overriding the adaptive time stepping next time,
      ! in case the coming time increment is cut-back in order to land on a special time point
      DELTME_lastBeforeSpecTime = DELTME
      
      if (loadType==1) then ! if monotonic creep, 
                              ! set the time step to the end of ramping time
         CALL T_CREEP_CSTRAINRATE(FINTME,DELTME,XMINDT,XMXDT)
         
      elseif(loadType==2) then ! if monotonic CSR, no special time-integration points..
         CONTINUE
      elseif(loadType==3) then ! traction-controlled cyclic loading
         if(P_dwell_load.EQ.1) then   ! dwell. 3 special points: end/beginning of ramping
         
            t_cyc = DMOD(FINTME,P_cyclic_period)   ! note: P_cyclic_period is the WHOLE load cycle period. i.e. in dwell, = t_hold+2*ramp
            
            if(t_cyc+XTOL.LT.P_dwell_ramptime) then   ! on ramp-up
               if(t_cyc+DELTME+XMINDT.GT.P_dwell_ramptime) then
                  DELTME = P_dwell_ramptime-t_cyc+XTOL
                  prevTimeCutBackForSpecial = .true.
                  
                  if(export_special.EQ.1 .and. &
                .not.suppress_specialPoints) EXPORT_next_step = .TRUE.
               endif
            elseif(t_cyc+XTOL.LT.P_cyclic_period-P_dwell_ramptime) then   ! on P_max
               if(t_cyc+DELTME+XMINDT.GT.P_cyclic_period-P_dwell_ramptime) then
                  
                  DELTME = P_cyclic_period-P_dwell_ramptime-t_cyc+XTOL
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1) EXPORT_next_step = .TRUE.
               endif
            elseif(t_cyc+XTOL.LT.P_cyclic_period) then ! on ramp-down
               if(t_cyc+DELTME+XMINDT.GT.P_cyclic_period) then
                  DELTME = P_cyclic_period-t_cyc
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1 .and. &
                .not.suppress_specialPoints) EXPORT_next_step = .TRUE.
               endif
            endif
         else              ! sinusoidal. 3 special points: max/min/zero

            P_cyclic_amp = 0.5D0*(P_cyclic_max-P_cyclic_min)
            P_cyclic_mean = 0.5D0*(P_cyclic_max+P_cyclic_min)
            tramp_cyc = P_cyclic_mean / (P_cyclic_amp*2*PI/P_cyclic_period)
            
            t_cyc = DMOD(FINTME-tramp_cyc, P_cyclic_period/4.d0)
            
            if(t_cyc+XMINDT > P_cyclic_period/4.d0) then  ! cannot cutback time to here. 
                                                          ! it is too close in front of the current solution.
                                                          ! (this should not be possible)
               if (export_special==1) write(*,*) 'unexpected landing spot -- cannot plot'
            else
               if(t_cyc+DELTME+XMINDT >= P_cyclic_period/4.d0) then
                  DELTME = P_cyclic_period/4.d0 - t_cyc + XTOL
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1) EXPORT_next_step = .TRUE.
               endif
            endif
            
         endif
      elseif(loadType==4) then ! displacement-controlled cyclic loading
         if(U_dwell_load.EQ.1) then   ! dwell. 3 special points: end/beginning of ramping
         
            t_cyc = DMOD(FINTME,U_cyclic_period)
            
            if(t_cyc+XTOL.LT.U_dwell_ramptime) then   ! on ramp-up
               if(t_cyc+DELTME+XMINDT.GT.U_dwell_ramptime) then
                  DELTME = U_dwell_ramptime-t_cyc+XTOL
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1 .and. &
                .not.suppress_specialPoints) EXPORT_next_step = .TRUE.
               endif
            elseif(t_cyc+XTOL.LT.U_cyclic_period-U_dwell_ramptime) then  ! on U_max
               if(t_cyc+DELTME+XMINDT.GT.U_cyclic_period-U_dwell_ramptime) then
                  DELTME = U_cyclic_period-U_dwell_ramptime-t_cyc+XTOL
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1) EXPORT_next_step = .TRUE.
               endif
            elseif(t_cyc+XTOL.LT.U_cyclic_period) then ! on ramp-down
               if(t_cyc+DELTME+XMINDT.GT.U_cyclic_period) then
                  DELTME = U_cyclic_period-t_cyc
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1 .and. &
                .not.suppress_specialPoints) EXPORT_next_step = .TRUE.
               endif
            endif
         else              ! sinusoidal. 3 special points: max/min/zero

            U_cyclic_amp = 0.5D0*(U_cyclic_max-U_cyclic_min)
            U_cyclic_mean = 0.5D0*(U_cyclic_max+U_cyclic_min)
            tramp_cyc = U_cyclic_mean / (U_cyclic_amp*2*PI/U_cyclic_period)
            
            t_cyc = DMOD(FINTME-tramp_cyc, U_cyclic_period / 4.d0)
            

            if(t_cyc+XMINDT > U_cyclic_period/4.d0) then  ! cannot cutback time to here. 
                                                          ! it is too close in front of the current solution.
                                                          ! (this should not be possible)
               if (export_special==1) write(*,*) 'unexpected landing spot -- cannot plot'
            else
               if(t_cyc+DELTME+XMINDT >= U_cyclic_period/4.d0) then
                  DELTME =  U_cyclic_period/4.d0 - t_cyc + XTOL
                  prevTimeCutBackForSpecial = .true.
                  if(export_special.EQ.1) EXPORT_next_step = .TRUE.
               endif
            endif
         endif         
      elseif(loadType==5) then   ! explicit traction-time history
         t_cyc = FINTME
         if(P_history_repeat.EQ.1) &
            t_cyc = DMOD(t_cyc,P_history_times(P_history_N))
         
         I=1
         do while(P_history_times(I).LE.t_cyc+XTOL)
            I=I+1
            if(I.GT.P_history_N) exit
         enddo
         
         if(.NOT.t_cyc.GT.P_history_times(I)) then ! if all time points have passed, skip.         
            if(t_cyc+DELTME+XMINDT.GT.P_history_times(I)) then
               DELTME = P_history_times(I) - t_cyc + XTOL
               prevTimeCutBackForSpecial = .true.
               if(mark_for_export(I).EQ.1) EXPORT_next_step = .TRUE.
            endif
         endif
      elseif(loadType==6) then   ! explicit displacement-time history
         t_cyc = FINTME
         if(U_history_repeat.EQ.1) &
            t_cyc = DMOD(t_cyc,U_history_times(U_history_N))
         
         I=1
         do while(U_history_times(I).LE.t_cyc+XTOL)
            I=I+1
            if(I.GT.U_history_N) exit
         enddo
         
         if(.NOT.t_cyc.GT.U_history_times(I)) then ! all time points have passed         
            if(t_cyc+DELTME+XMINDT.GT.U_history_times(I)) then
               DELTME = U_history_times(I) - t_cyc + XTOL
               prevTimeCutBackForSpecial = .true.
               if(mark_for_export(I).EQ.1) EXPORT_next_step = .TRUE.
            endif
         endif
      elseif(loadType==7) then
         ! not yet supported
         write(*,*) 'STRESS-TIME HISTORY NOT YET SUPPORTED'
         STOP
         
      elseif(loadType==8) then
         ! explicit macroscopic strain-time history
         t_cyc = FINTME
         if(strain_history_repeat.EQ.1) &
            t_cyc = DMOD(t_cyc,strain_history_times(strain_history_N))
         
         I=1
         do while(strain_history_times(I).LE.t_cyc+XTOL)
            I=I+1
            if(I.GT.strain_history_N) exit
         enddo
         
         if(.NOT.t_cyc.GT.strain_history_times(I)) then ! skip if all time points have passed
            if(t_cyc+DELTME+XMINDT.GT.strain_history_times(I)) then
               DELTME = strain_history_times(I) - t_cyc + XTOL
               prevTimeCutBackForSpecial = .true.
               if(mark_for_export(I).EQ.1) EXPORT_next_step = .TRUE.
            endif
         endif
      elseif(loadType==9 .or. &
             loadType==10) then  ! if mixed macroscopic stress-strain controlled, 
                                 ! deformation is effective-strain controlled and monotonic
                                 ! => no special time-integration points.
         CONTINUE
         
      endif
      
      ! must be larger than min allowed dt
      DELTME=MAX(DELTME,XMINDT)

      ! stop at the requested simulation time and also force export
      if (FINTME+DELTME.GT.TOTME) then
         EXPORT_next_step = .TRUE.
         DELTME = TOTME - FINTME + tiny(DELTME)
      endif

      END SUBROUTINE
      
      SUBROUTINE trap_FPE()
      implicit none
      logical :: trap_error
      COMMON/FPE_TRAP/trap_error

      trap_error = .TRUE.
      write(*,*) 'error trapped'
      
      END SUBROUTINE
      
      ! Calculates the equivalent finite strain (scalar) from a deformation gradient tensor:
      ! E_GL = 1/2*[F^T*F - I]
      ! E_equiv = Sqrt(2/3*||E_GL||)
      SUBROUTINE calcE_Equiv(F,E_Equiv)
      implicit none
      real(8), intent(in) :: F(9)
      real(8), intent(out):: E_Equiv
      
      real(8) :: F_2D(3,3), E_GreenLagrange(3,3)
      
      F_2D(1,1:3) = F(1:3)
      F_2D(2,1:3) = F(4:6)
      F_2D(3,1:3) = F(7:9)
      
      E_GreenLagrange = MATMUL(TRANSPOSE(F_2D),F_2D)
      E_GreenLagrange(1,1) = E_GreenLagrange(1,1) - 1.D0
      E_GreenLagrange(2,2) = E_GreenLagrange(2,2) - 1.D0
      E_GreenLagrange(3,3) = E_GreenLagrange(3,3) - 1.D0
      E_GreenLagrange = 0.5D0 * E_GreenLagrange
      
      E_Equiv = DSqrt(2.0D0/3.0D0*(E_GreenLagrange(1,1)**2 + &
                                   E_GreenLagrange(1,2)**2 + &
                                   E_GreenLagrange(1,3)**2 + &
                                   E_GreenLagrange(2,1)**2 + &
                                   E_GreenLagrange(2,2)**2 + &
                                   E_GreenLagrange(2,3)**2 + &
                                   E_GreenLagrange(3,1)**2 + &
                                   E_GreenLagrange(3,2)**2 + &
                                   E_GreenLagrange(3,3)**2 ))
      END SUBROUTINE
      

      SUBROUTINE calcElemGMATandVol(XBAR_ST,XG_ST,DVOLTAU, &
          COORDS_E,VTAU_E,DVTAU_E,DFGRD0,DFGRD1, &
          DVOL0,DVOLT,DET_FTAU,DET_FT,SHPTAU_ST,success)
     
      implicit none
      !include 'ivep.h'
      integer, parameter :: IPT=1 ! # of gauss pts
      integer, parameter :: NNODE=4 ! # nodes per element
      integer, parameter :: NDF=3
      integer, parameter :: MAXDIM=3
      
      ! arguments
      real(8), intent(in) :: COORDS_E(NDF,NNODE)
      real(8), intent(in) :: VTAU_E(NDF,NNODE), DVTAU_E(NDF,NNODE)
      real(8), intent(out) :: XBAR_ST(2*MAXDIM,MAXDIM*NNODE,IPT)
      real(8), intent(out) :: XG_ST(MAXDIM**2,MAXDIM*NNODE,IPT)
      real(8), intent(out) :: SHPTAU_ST(3,NNODE,IPT)
      real(8), intent(out) :: DVOLTAU,DVOLT,DVOL0
      real(8), intent(out) :: DFGRD0(3,3),DFGRD1(3,3),DET_FTAU,DET_FT
      logical, intent(out) :: success
      ! locals
      real(8) :: SHPTAU(3,NNODE),SHPT(3,NNODE)
      real(8) :: XE_TAU(MAXDIM,NNODE), XE_T(MAXDIM,NNODE)
      real(8) :: VT_E(NDF,NNODE)
      real(8) :: DFGRD1_ST(3,3,IPT),DFGRD0_ST(3,3,IPT)
      real(8) :: DVOL0_ST(1),DVOLT_ST(1),DVOLTAU_ST(1)
      real(8) :: SHP0(3,NNODE),SHP0_ST(3,NNODE,IPT), &
       XB_DEV(6,12,1),XB_VOL(6,12,1),XG_DEV(9,12,1), &
       XG_VOL(9,12,1)   

      success = .true.
      VT_E = VTAU_E - DVTAU_E
      
      XE_TAU = COORDS_E + VTAU_E
      XE_T = COORDS_E + VT_E

      CALL SHAPE_TET4_Kourosh(COORDS_E,SHP0,DVOL0)
      CALL SHAPE_TET4_Kourosh(XE_TAU,SHPTAU,DVOLTAU)
      CALL SHAPE_TET4_Kourosh(XE_T,SHPT,DVOLT)

      IF(DVOLTAU < 0.d0 .or. DVOL0 < 0.d0)THEN
         success = .false.
         return
      ENDIF
      
      DVOL0_ST(1)=DVOL0
      DVOLTAU_ST(1)=DVOLTAU
      DVOLT_ST(1)=DVOLT
      SHP0_ST(1:NDF,1:NNODE,1)=SHP0(1:NDF,1:NNODE)
      SHPTAU_ST(1:NDF,1:NNODE,1)=SHPTAU(1:NDF,1:NNODE)
 
      CALL CALC_DFG_TET4_Kourosh(VT_E,NNODE,NDF, &  ! check whether DVOL0, or DVOLT
                         DFGRD0_ST,DVOL0_ST,SHP0_ST,1,1)

      CALL CALC_DFG_TET4_Kourosh(VTAU_E,NNODE,NDF, & ! check whether DVOL0, or DVOLTAU
                         DFGRD1_ST,DVOL0_ST,SHP0_ST,1,1)
     
      CALL MAKEGRAD_TET4_Kourosh(XBAR_ST,XG_ST,SHPTAU_ST, &
                         XB_DEV,XB_VOL,XG_DEV,XG_VOL)

      DFGRD1(1:3,1:3)=DFGRD1_ST(1:3,1:3,1)
      DFGRD0(1:3,1:3)=DFGRD0_ST(1:3,1:3,1)

      CALL DET33(DFGRD1,DET_FTAU)
      CALL DET33(DFGRD0,DET_FT)
 
      RETURN
      END

!-------------------------

      SUBROUTINE DET33(A,DET)
      IMPLICIT NONE
      REAL(8):: DET
      REAL(8):: A(3,3)
      
      DET=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2) &
           *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1) &
           *A(1,3)*A(2,2))

      RETURN
      END
       
!------------------------------------------------------
      SUBROUTINE CALC_DFG_TET4_Kourosh( &
         UE,NNODE,MCRD,DFGRD_ST,        &
      	DVOL_ST,SHP_ST,IBAR,NGAUSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UE(MCRD,NNODE), &
         XJAC_ST(NNODE),DFGRD_ST(3,3,1),DFGRD_INV(3,3), &
         DVOL_ST(1),SHP_ST(3,4,1),DFGRD(3,3)

      XJBAR=0.D0
      TV=0.D0
      DO IPT=1,NGAUSS
         
         DO I=1,MCRD
         DO J=1,MCRD
            DFGRD(I,J)=0.D0
            IF(I.EQ.J)DFGRD(I,J)=1.D0
            DO IN=1,NNODE
            DFGRD(I,J)=DFGRD(I,J)+SHP_ST(J,IN,IPT)*UE(I,IN)
            ENDDO
         ENDDO
         ENDDO
         
         CALL MATINV3_UEL(DFGRD,DFGRD_INV,XJAC)
         XJAC_ST(IPT)=XJAC
         DVOL=DVOL_ST(IPT)
         XJBAR=XJBAR+XJAC*DVOL
         TV=TV+DVOL
         DFGRD_ST(1:3,1:3,IPT)=DFGRD(1:3,1:3)
      ENDDO
      XJBAR=XJBAR/TV
      IF(IBAR.EQ.1)THEN
         DO IPT=1,NGAUSS
            FAC=XJBAR**(1.D0/3.D0)*XJAC_ST(IPT)**(-1.D0/3.D0)
            DFGRD_ST(1:3,1:3,IPT)=DFGRD_ST(1:3,1:3,IPT)*FAC
         ENDDO
      ENDIF
      RETURN
      END   

      SUBROUTINE MAKEGRAD_TET4_Kourosh(XB_ST,XG_ST,SHPTAU_ST, &
                               XB_DEV,XB_VOL,XG_DEV,XG_VOL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XB_ST(6,12,1),SHPTAU_ST(3,4,1),XG_ST(9,12,1), &
                 H(6,9),XB_DEV(6,12,1),XB_VOL(6,12,1),XBTR_ST(12,1), &
                 XG_DEV(9,12,1),XG_VOL(9,12,1)

      XG_ST=0.D0
      XB_ST=0.D0

      XG_ST(1:3,1,1)=SHPTAU_ST(1:3,1,1)
      XG_ST(4:6,2,1)=SHPTAU_ST(1:3,1,1)
      XG_ST(7:9,3,1)=SHPTAU_ST(1:3,1,1)

      XG_ST(1:3,4,1)=SHPTAU_ST(1:3,2,1)
      XG_ST(4:6,5,1)=SHPTAU_ST(1:3,2,1)
      XG_ST(7:9,6,1)=SHPTAU_ST(1:3,2,1)

      XG_ST(1:3,7,1)=SHPTAU_ST(1:3,3,1)
      XG_ST(4:6,8,1)=SHPTAU_ST(1:3,3,1)
      XG_ST(7:9,9,1)=SHPTAU_ST(1:3,3,1)

      XG_ST(1:3,10,1)=SHPTAU_ST(1:3,4,1)
      XG_ST(4:6,11,1)=SHPTAU_ST(1:3,4,1)
      XG_ST(7:9,12,1)=SHPTAU_ST(1:3,4,1)

      H=0.D0
      H(1,1)=1.D0
      H(2,5)=1.D0
      H(3,9)=1.D0
      H(4,2)=1.D0
      H(4,4)=1.D0
      H(5,3)=1.D0
      H(5,7)=1.D0
      H(6,6)=1.D0
      H(6,8)=1.D0
      
      DO I=1,6
         DO J=1,12
            DO K=1,9
               XB_ST(I,J,1)=XB_ST(I,J,1)+H(I,K)*XG_ST(K,J,1)
            ENDDO
         ENDDO
      ENDDO
      
      	   
      XB_DEV=XB_ST  
      XB_VOL=0.0D0    

      DO J=1,12
         XBTR_ST(J,1)=0.D0
         DO I=1,3
            XBTR_ST(J,1)=XBTR_ST(J,1)+XB_ST(I,J,1)
         ENDDO
         DO I=1,3
            XB_DEV(I,J,1)=XB_DEV(I,J,1)-1.D0/3.D0*XBTR_ST(J,1)           !JIAHAO
            XB_VOL(I,J,1)=1.D0/3.D0*XBTR_ST(J,1)
         ENDDO
      ENDDO

      DO J=1,12
         JTEMP=(J-1)/3
         JTEMP1=JTEMP+1
         JTEMP2=J-JTEMP*3
         DO I=1,9
            ITEMP=(I-1)/3
            ITEMP1=ITEMP+1
            ITEMP2=I-ITEMP*3

            XG_DEV(I,J,1)=0.D0
            XG_VOL(I,J,1)=0.D0
            IF(JTEMP2.EQ.ITEMP1) THEN
               XG_DEV(I,J,1)=SHPTAU_ST(ITEMP2,JTEMP1,1)
            ENDIF

            IF(ITEMP2.EQ.ITEMP1)THEN
               XG_DEV(I,J,1)=XG_DEV(I,J,1)-1.D0/3.D0*XBTR_ST(J,1)
               XG_VOL(I,J,1)=XG_VOL(I,J,1)+1.D0/3.D0*XBTR_ST(J,1)
            ENDIF
         ENDDO
      ENDDO	      

      RETURN
      END
      
      SUBROUTINE SHAPE_TET4_Kourosh(XEL,SHP,DVOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHP1(3,4),SHP(3,4),XJAC(3,3),XJACINV(3,3),XEL(3,4), &
       XX(3),YY(3),ZZ(3)

      SHP1=0.D0
      SHP=0.D0

      SHP1(1,1)=1.D0
      SHP1(2,2)=1.D0
      SHP1(3,3)=1.D0

      SHP1(1,4)=-1.D0
      SHP1(2,4)=-1.D0
      SHP1(3,4)=-1.D0

      DO I=1,3
         XX(I)=XEL(1,I)-XEL(1,4)
         YY(I)=XEL(2,I)-XEL(2,4)
         ZZ(I)=XEL(3,I)-XEL(3,4)
      ENDDO

      XJAC(1:3,1)=XX(1:3)
      XJAC(1:3,2)=YY(1:3)
      XJAC(1:3,3)=ZZ(1:3)

      CALL MATINV3_UEL(XJAC,XJACINV,DVOL)
      IF(DVOL.LT.1D-16)THEN
         RETURN
      ENDIF

      DVOL=1.D0/6.D0*DABS(DVOL)
      
      DO INODE=1,4
         SHP(:,INODE)=MATMUL(XJACINV,SHP1(:,INODE))
      ENDDO

      RETURN
      END
      
      subroutine convertDDSDDE_EngToMat(DDSDDE)
      implicit none
      real(8), intent(inout) :: DDSDDE(6,6)
      ! strain definition -- mathematical vs. engineering strain
      DDSDDE(:,4) = DDSDDE(:,4) * 2.d0
      DDSDDE(:,5) = DDSDDE(:,5) * 2.d0
      DDSDDE(:,6) = DDSDDE(:,6) * 2.d0
      
      !DDSDDE(4,:) = DDSDDE(4,:) * 2.d0
      !DDSDDE(5,:) = DDSDDE(5,:) * 2.d0
      !DDSDDE(6,:) = DDSDDE(6,:) * 2.d0
      end subroutine

      subroutine updateStrainIterator( &
                     incStrainMacroscopic_inout,stressMacroscopic_in,ddsddeMacroscopic_in, &
                     stressResidual,stressReaction,errorVM,lambda, &
                     incrVonMises_imposed,biaxial_StressRatioAngle, &
                     iStressMain,iStressLateral, &
                     specifiedLoadFrame,R_loadVoigt,R_loadVoigt_Inv)
      use options, only: DEBUG_GENERAL, DEBUG_UMAT
      implicit none
      ! arguments
      real(8), intent(in)   :: stressMacroscopic_in(6), ddsddeMacroscopic_in(6,6)
      real(8), intent(out)  :: stressResidual(6), stressReaction, errorVM
      real(8), intent(inout):: incStrainMacroscopic_inout(6)
      real(8), intent(inout):: lambda
      real(8), intent(in) :: incrVonMises_imposed
      real(8), intent(in) :: biaxial_StressRatioAngle
      integer, intent(in) :: iStressMain, iStressLateral
      logical, intent(in) :: specifiedLoadFrame
      real(8), intent(in) :: R_loadVoigt(6,6),R_loadVoigt_Inv(6,6)
      
      ! locals
      real(8) :: stressMacroscopic(6),ddsddeMacroscopic(6,6)
      real(8) :: incStrainVM
      real(8) :: incStrainMacropscopic_prev(6),incStrainMacroscopic(6)
      integer, parameter :: nStrainConstraints=1,NTENS=6,nRHS=1
      real(8) :: maxRes(nStrainConstraints), maxCon(nStrainConstraints), maxUpd(nStrainConstraints)
      real(8) :: grad_strainConstraint(nStrainConstraints,NTENS)
      real(8) :: R_strainConstraint(nStrainConstraints,nRHS)
      
      ddsddeMacroscopic = ddsddeMacroscopic_in
      stressMacroscopic = stressMacroscopic_in
      incStrainMacroscopic = incStrainMacroscopic_inout
      
      call convertDDSDDE_EngToMat(ddsddeMacroscopic)

      ! (if loading is in arbitrary direction)
      ! switch to the coordinate system of the desired load frame
      !----------------------------------------------------------
      if (specifiedLoadFrame) then
         ! transform the strain tensor to the load coordinate frame
         ! in order to impose the constraints:
         stressMacroscopic = MATMUL(R_loadVoigt_Inv,stressMacroscopic)
         ! transform the tangent stiffness too
         ddsddeMacroscopic = MATMUL(MATMUL(R_loadVoigt_Inv,ddsddeMacroscopic),R_loadVoigt)
         ! transform the strain increment too
         incStrainMacroscopic = MATMUL(R_loadVoigt_Inv,incStrainMacroscopic)
      endif
      !----------------------------------------------------------
      
      stressResidual = - stressMacroscopic
      
      ! modify the residual and tangent matrix to impose constraint on biaxial stress ratio
      call ImposeBiaxialConstraints(stressResidual,ddsddeMacroscopic,stressReaction, &
           biaxial_StressRatioAngle,iStressMain,iStressLateral)
                                    
      ! calculate the residual and gradient of the VM-constraint
      call getVonMisesStrainVoigt_MathStrain(incStrainVM,incStrainMacroscopic)
      call getVonMisesStrainConstraintVoigt_MathStrain( &
               R_strainConstraint,grad_strainConstraint,incrVonMises_imposed,incStrainMacroscopic)

      ! apply the newton update for the constrained problem
      call incrementNRwConstraints( &
               incStrainMacroscopic,lambda, &
               stressResidual,ddsddeMacroscopic, &
               grad_strainConstraint,R_strainConstraint, &
               NTENS,nStrainConstraints,nRHS, &
               maxRes,maxCon,maxUpd)
               
      ! (if loading is in arbitrary direction)
      ! switch back to the global coordinate frame
      ! before returning output variables
      !----------------------------------------------------------
      if (specifiedLoadFrame) then
         ! transform the strain increment
         incStrainMacroscopic = MATMUL(R_loadVoigt,incStrainMacroscopic)
      endif
      !----------------------------------------------------------
               
      ! save updated strain increment to output variable
      incStrainMacroscopic_inout = incStrainMacroscopic

      ! error in von mises strain
      errorVM = abs(incStrainVM - incrVonMises_imposed)
                        
      end subroutine
      
      subroutine ImposeBiaxialConstraints(res,DDSDDE,R_rxn, &
                 biaxial_StressRatioAngle,iStressMain,iStressLateral)
      implicit none
      real(8), intent(inout) :: res(6), DDSDDE(6,6)
      real(8), intent(out)   :: R_rxn
      real(8), intent(in)    :: biaxial_StressRatioAngle
      integer, intent(in)    :: iStressMain,iStressLateral
      ! locals
      real(8) :: DDSDDE_row_Main(6),DDSDDE_row_Lateral(6)
      real(8) :: DDSDDE_StressConstraint(6)
      real(8) :: res_Main,res_Lateral,resStressConstraint
      
      res_Main = res(iStressMain)
      res_Lateral = res(iStressLateral)
      resStressConstraint = 0.5d0*(cos(biaxial_StressRatioAngle)*res_Lateral - &
                                   sin(biaxial_StressRatioAngle)*res_Main )
      res(iStressMain) = resStressConstraint
      res(iStressLateral) = resStressConstraint
      
      DDSDDE_row_Main(:) = DDSDDE(iStressMain,:)
      DDSDDE_row_Lateral(:) = DDSDDE(iStressLateral,:)
      DDSDDE_StressConstraint = 0.5d0*(cos(biaxial_StressRatioAngle)*DDSDDE_row_Lateral - &
                                       sin(biaxial_StressRatioAngle)*DDSDDE_row_Main )
      DDSDDE(iStressMain,:) = DDSDDE_StressConstraint
      DDSDDE(iStressLateral,:) = DDSDDE_StressConstraint
      
      R_rxn = sqrt(res_Main**2+res_Lateral**2)
      
      end subroutine
      

      ! applies a Newton-Raphson increment for a non-linear problem with
      ! linear constraints
      ! 
      ! Author: Deniz Ozturk, 2017
      SUBROUTINE incrementNRwConstraints(x,lambda,res,H,A,resCon,n,m,nP,maxRes,maxCon,maxUpd)      
         
!      use LibraryMath, only:     &
!         solveMatrixEquation
         
      implicit none
      ! arguments
      integer, intent(in) :: n,m,nP
      real(8), intent(inout) :: x(n,nP),lambda(m,nP)
      real(8), intent(in) :: res(n,nP) ! the residual of the problem. sign convention: res = -f_int + f_ext, i.e. external forces go into the residual with a positive sign.
      real(8), intent(in) :: H(n,n)    ! the Hessian (or the stiffness matrix) for the unconstrained problem
      real(8), intent(in) :: A(m,n)    ! the linear constraint-coefficient matrix for linear problems, or the gradient of the constraint function for nonlinear problems.
      real(8), intent(in) :: resCon(m,nP)  ! ... and the value of non-linear constraint function at the, c(x). m constraints are given as: c(x) = 0
      !real(8), intent(in) :: b(m,nP)   ! ... and the constraint vector. m constraints are given as: A.x - b = 0
      real(8), intent(out):: maxRes(nP),maxCon(nP),maxUpd(nP)  ! (OUT) max residuals of the constrained problem and max newton-updates used in convergence criterion
      ! locals
      real(8) :: KKT(n+m,n+m),rhs(n+m,nP)
      integer :: iP
      integer :: ierror
      character(len=3) :: strRepeat
      
      ! construct the KKT matrix
      KKT(1:n,1:n) = H
      KKT(1:n,n+1:n+m) = transpose(A)
      KKT(n+1:n+m,1:n) = A
      KKT(n+1:n+m,n+1:n+m) = 0.d0
      
      ! construct the RHS for the constrained increment
      do iP = 1, nP 
         rhs(1:n,iP) = res(1:n,iP) - MATMUL(transpose(A),lambda(:,iP))
         rhs(n+1:n+m,iP) = -(resCon(:,iP))
      enddo
      
      ! find the maximum residual for the constrained problem
      do iP = 1, nP
         maxRes(iP) = maxval(abs(res(1:n,iP)))
         maxCon(iP) = maxval(abs(resCon(:,iP)))
      enddo
      
      !write(900,*) 'KKT RHS:'
      !write(900,*) rhs(:,1)
      
      do iP = 1, nP       
         CALL solveAXb_LU(KKT,n+m,rhs(:,iP),ierror)
         if (ierror /= 0) then
            write(*,*) 'lapack error in KKT solver'
            stop
         endif
!         rhs(:,iP) = solveMatrixEquation(KKT,rhs(:,iP))     
      enddo
      
      ! find the maximum update of the DOFs and constraint forces
      do iP = 1, nP
         maxUpd(iP) = maxval(abs(rhs(:,iP)))
      enddo
      
      !write(900,*) 'KKT update:'
      !write(900,*) rhs(:,1)
      
      ! apply the newton update
      do iP = 1, nP 
         x(1:n,iP) = x(1:n,iP) + rhs(1:n,iP)
         lambda(1:m,iP) = lambda(1:m,iP) + rhs(n+1:n+m,iP)
      enddo
      
      END SUBROUTINE
      

      subroutine transformAffineDisplacementField(v,macroStrain_old,macroStrain_new,G0XYZ,domainRange0,nNodes)
      use options, only: useLogStrainRate
      implicit none
      integer, intent(in)   :: nNodes
      real(8), intent(inout):: v(3*nNodes)
      real(8), intent(in)   :: macroStrain_old(6),macroStrain_new(6)
      real(8), intent(in)   :: domainRange0(3,2)
      real(8), intent(in)   :: G0XYZ(3*nNodes)
      
      ! locals
      real(8) :: macroF_old(3,3),macroF_new(3,3),dF(3,3)
      real(8) :: domainCenter(3)
      real(8) :: nodePosRel(3)
      integer :: iNode,offsetDOFnode
      
      domainCenter(:) = 0.5D0*(domainRange0(:,1) + domainRange0(:,2))
      
      
      ! calculate macroF_tau associated with macroStrain_tau
      if (useLogStrainRate) then
         call getFfromE_log(macroF_new,macroStrain_new)
         call getFfromE_log(macroF_old,macroStrain_old)
      else
         call getFfromE_GL(macroF_new,macroStrain_new)
         call getFfromE_GL(macroF_old,macroStrain_old)
      endif
      
      dF = macroF_new - macroF_old

      do iNode=1,nNodes
         offsetDOFnode = 3*(iNode-1)
         nodePosRel(1:3) = G0XYZ(offsetDOFnode+1:offsetDOFnode+3) - domainCenter(1:3)
         
         ! apply affine transformation association with increment of deformation gradient
         v(offsetDOFnode+1:offsetDOFnode+3) = &
               v(offsetDOFnode+1:offsetDOFnode+3) + matmul(dF,nodePosRel)
      enddo
         
      end subroutine
                  
      ! calculates F from Green-Lagrange strain, assuming F=RU and R=I
      ! E_voigtMath(6) is Green-Lagrange in Voigt form with mathematical shear components
      SUBROUTINE getFfromE_GL(F,E_voigtMath)
      
      implicit none
      real(8), intent(out) :: F(3,3)
      real(8), intent(in) :: E_voigtMath(6)
      
      ! locals
      real(8), parameter :: Identity(3,3) = &
         reshape((/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0 /), (/3,3/))
      real(8) :: E(3,3),C(3,3),stretch(3)
      real(8) :: U(3,3)
      real(8) :: valEigen_C(3),vecEigen_C(3,3)
      integer :: i,j,k
      integer :: iError
         
      ! use F=R.U, where R=I. express U in terms of the given E.
      
      Call sigmat(E_voigtMath,E)
      
      C = 2*E+Identity
      
      CALL eigen_Lapack(valEigen_C,vecEigen_C,C,3,iError)
      
      IF (iError /= 0) THEN
         write(*,*) 'Lapack EigenProblem Error:',iError
         write(*,*) 'Matrix (C):',C
         stop
      endif
      
      stretch(1) = sqrt(valEigen_C(1))
      stretch(2) = sqrt(valEigen_C(2))
      stretch(3) = sqrt(valEigen_C(3))
      
      U=0.d0
      do k=1,3 ! eigenvalues
         do i=1,3
            do j=1,3
               U(i,j) = U(i,j) + stretch(k)*vecEigen_C(i,k)*vecEigen_C(j,k)
            enddo
         enddo
      enddo
      
      F = U ! F=R*U, where R=I.
      
      END SUBROUTINE

      ! calculates F from logarithmic strain, assuming F=RU and R=I
      ! E_voigtMath(6) is Log-Strain in Voigt form with mathematical shear components
      SUBROUTINE getFfromE_log(F,E_voigtMath)
      
      implicit none
      real(8), intent(out) :: F(3,3)
      real(8), intent(in) :: E_voigtMath(6)
      
      ! locals
      real(8) :: E(3,3),C(3,3),stretch(3)
      real(8) :: U(3,3)
      real(8) :: valEigen_E(3),vecEigen_E(3,3)
      integer :: i,j,k
      integer :: iError
         
      ! use F=R.U, where R=I. express U in terms of the given E.
      
      Call sigmat(E_voigtMath,E)
      
      CALL eigen_Lapack(valEigen_E,vecEigen_E,E,3,iError)
      
      IF (iError /= 0) THEN
         write(*,*) 'Lapack EigenProblem Error:',iError
         write(*,*) 'Matrix (LogStrain):',E
         stop
      endif
      
      stretch(1) = exp(valEigen_E(1))
      stretch(2) = exp(valEigen_E(2))
      stretch(3) = exp(valEigen_E(3))
      
      U=0.d0
      do k=1,3 ! eigenvalues
         do i=1,3
            do j=1,3
               U(i,j) = U(i,j) + stretch(k)*vecEigen_E(i,k)*vecEigen_E(j,k)
            enddo
         enddo
      enddo
      
      F = U ! F=R*U, where R=I.
      
      END SUBROUTINE
      

      ! returns the euler angles Phi1 and Theta,
      ! corresponding to a transformation that 
      ! brings the Z axis along a given vector n
      !
      ! Author: Deniz Ozturk, 2016
      !----------------------------
      SUBROUTINE getPhi1andTheta(nx,ny,nz,phi1,theta)
      implicit none
      real(8), intent(in) :: nx,ny,nz
      real(8), intent(out) :: phi1,theta
      real(8), parameter :: PI=3.14159265359D0
      real(8) :: r

      if (nx==0.d0 .and. ny==0.d0 .and. nz==0.d0) then
         phi1 = 0.d0
         theta = 0.d0
         return
      endif

      r = dsqrt(nx*nx+ny*ny+nz*nz)
      
      theta = acos(nz/r)
      
      if (nx==0.d0 .and. ny==0.d0) then   ! singularity here
         phi1=0.d0
      else
         phi1 = acos(-ny/DSQRT(nx*nx+ny*ny))
      endif
      if(nx.LT.0.D0) phi1 = 2.D0*PI-phi1

      END SUBROUTINE