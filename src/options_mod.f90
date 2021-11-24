      ! Options Module for the CPFE code
      !
      ! Author: Deniz Ozturk (2015)
      !-----------------------------------------------------------------
      ! input file format:
      ! [variable name] [value]
      ! to leave comments in the option file, at any line:
      ! // comments...
      !
      ! notes: 
      !   supported input types: Integer, Floating Pt Number, String, Boolean (True/False)
      !   for boolean (true/false) variables, enter 1, 0 as values.
      !   there are default values for each option/variable (set these in options_SetDefaults()),
      !   options take these default values unless their values are set in the options/input file.
      !   use double slash to insert comments anywhere into your options/input file.
      !   these parameters can be modified to meet requirements:
      !        maxNumOptions     : max. # of options read from the options file
      !        maxStringLength   : max. length of option names and string input-values
      !
      MODULE options

      integer, parameter :: maxNumOptions = 200
      integer, parameter :: maxStringLength = 32
      
      private :: options_SetDefaults
      private :: parseString
      private :: toLowerCase

      ! program options (will be imported from options.inp)
      character(len=maxStringLength) :: option_names(maxNumOptions)
      integer :: options_Int(maxNumOptions)
      real(8) :: options_Float(maxNumOptions)
      logical :: options_Bool(maxNumOptions)
      character(len=maxStringLength) :: options_String(maxNumOptions)
      integer :: nOptions
      
      ! -------------------------------------------
      ! Option Values
      ! -------------------------------------------
      ! StiffnessMatrixFormulation:
      integer, parameter :: SM_delSdelLnU = 1, &   
                            SM_Jiahao = 2, &
                            SM_delPdelF = 3
                            
      ! solverMethod     
      integer, parameter :: solmeth_newtRaph=1, &
                            solmeth_modNewt=2, &
                            solmeth_modNewtKeep=3, &
                            solmeth_BFGS=4, &
                            solmeth_BFGSkeep=5, &
                            solmeth_Broyden=6
      ! specified the nonlinear solver method: 1:newton-Raphson 2:modif Newt 3:modif Newt (Keep K) 4: BFGS 5: BFGS (Keep K) 6: Broyden?
      ! -------------------------------------------
      ! PROGRAM options
      ! -------------------------------------------
      logical :: DEBUG_general     ! export debugging information during calculation of nodal values
      logical :: DEBUG_UMAT        ! export UMAT debugging information
      logical :: DEBUG_grains      ! export debug info for grain identification
      logical :: DEBUG_errInd      ! export debug info for error indicators
      logical :: DEBUG_nodalValues ! export debugging information during calculation of nodal values
      logical :: DEBUG_timer       ! export timing information
      logical :: DEBUG_RESID       ! export residual for debugging
      logical :: DEBUG_linesearch  ! export line search reports for debugging
      real(8) :: postProcT, dataExportT
      integer :: postProcN, dataExportN
      logical :: do_errAnalysis    ! whether to do (on the run) error analysis/adaptivity
      logical :: exportSolution    ! export solution at the end of the calculation to continuing later on
      logical :: continueSolution  ! read-in and continue solution
      logical :: verboseUmat       ! umat exports detailed error messages for each element
      logical :: conventionDeka    ! use deka's convention for Euler angles (euler angles describe an active rotation of the crystal frame onto the global frame)
      logical :: EXPORT_watmus     ! export inputs fir WATMUS-acceleration
      logical :: EXPORT_stiffnessMatrix ! export the initial stiffness matrix and quit program
      logical :: EXPORT_delta_gamma! export delta_gamma.out
      logical :: EXPORT_backstress ! export back-stress on the 3 basal slip systems
      logical :: EXPORT_temperature! export local temperature
      logical :: EXPORT_wp         ! export dissipated energies: dissipated energy density, back-stress dissipation on basal plane and the 3 prismatic planes
      logical :: EXPORT_damage     ! export damage.out, damage related variables
      logical :: EXPORT_fp         ! export fp.out
      logical :: EXPORT_creep      ! export local GL plastic strain (creep.out)
      logical :: EXPORT_creepRate  ! export rate of local GL plastic strain (creepRate.out)
      logical :: EXPORT_strainRate ! export rate of local total strain
      logical :: EXPORT_nye        ! export Nye tensor
      logical :: EXPORT_nodePos    ! export nodal positions
      logical :: EXPORT_cauchy     ! export nodal positions
      logical :: EXPORT_grainAVG   ! export grain-averaged quantities
      logical :: EXPORT_grainRot   ! export grain rotation matrices (for plotting textures)
      logical :: EXPORT_hardening  ! export max basal/prismatic hardening parameters (CRSS)
      logical :: EXPORT_specialOnly! data export is done only at the special time points of a load form
                     !  (extrema/zeroes of sine, begin/end/zeros of dwell and times marked as 'export' in explicit load histories)
                     !  this option overrides dataExportT/dataExportN
      logical :: suppress_specialPoints   ! if set, the special points of a load waveform are limited to a single, maximum-load point
      real(8) :: EXPORT_allAtLast  ! export all variables at the end of the simulation
                     ! this can be used to limit disk usage. set EXPORT_ variables of only the essential variables to true
      logical :: UMAT_lsch         ! enables line search in umat solver
      logical :: UMAT_halfStep     ! dampens NR iteration step by half for the first iteration to avoid overshooting the flow rule (for CPFE umat)
      logical :: UMAT_stabilize    ! if second level iteration takes long to converge, dampens the delta_g update steps. also suppresses h0 values if g is close to g_saturation.
      real(8) :: TAU_LIMIT         ! threshold value of tau/g over which the flow rule will be assumed linear, to accelerate the convergence of iterations.
      integer :: improve_U0        ! improve initial trial disp. field. 0: no improv. 1: Jiahao's 2. Deniz's.
      real(8) :: STEP_v            ! velocity field factor for initial disp. field guess
      real(8) :: WLF_c2            ! Wolfe-2 condition coefficient
      real(8) :: STEP_max          ! max allowed step length in line search
      integer :: lschMethod        ! specifies the line search algorithm to use. 1: deniz 2: matthies, strang, 1979
      integer :: solverMethod      ! specified the nonlinear solver method: 1:newton-Raphson 2:modif Newt 3:modif Newt (Keep K) 4: BFGS 5: BFGS (Keep K) 6: Broyden?
      logical :: useTangentStiffnessStrainUpdates
      logical :: useLogStrainRate
      real(8) :: endSimulationAtWp ! ends simulation when plastic work exceeds specified value
      integer :: StiffnessMatrixFormulation  ! tangent matrix 1=ddS/ddE, 2=Not-Supported, 3=delPdelF (as in the FBar paper, default)
      logical :: SymmetrizeStiffnessMatrix   ! done at the element-level
      logical :: adaptive_time_stepping
      integer :: convergenceNorm   ! 1: L_max, 2: L_2 norm used in calculating tolerance for convergence
      real(8) :: eps_F, eps_u      ! convergence tolerance/precision on force and displacement      

      ! -------------------------------------------
      ! options for MATERIAL BEHAVIOR
      ! -------------------------------------------
      character(len=maxStringLength) :: materialName  ! material properties file to use. default: Ti6424 
      logical :: ELASTIC_ONLY ! disable plastic behavior (works by setting a_dot=0)
      integer :: methodTensCompCriterion ! method used for determining the tension-compression state.
      logical :: smoothTCtransition ! whether to use smooth T/C transition or discrete T/C state
      integer :: tensOnly ! for stress-dependent material properties, only use params for tension(+1), compression(-1), automatic(0)
      integer :: softOnly ! for flow-direction dependent material properties (BCC), only use values for 'soft'(+1) direction, hard (-1), automatic (0)
      logical :: HARDENING_DISABLED ! disable hardening
      logical :: BACKSTRESS_DISABLED ! disable back stress (e.g. single crystals)
      logical :: HALLPETCH_DISABLED ! disable hall petch effect
      logical :: PRISMATIC_ONLY ! disable all slip systems but prismatic-a's
      logical :: BASAL_ONLY     ! disable all slip systems but basal-a's
      logical :: GND_DISABLED ! disable GND hardening
      logical :: USE_POWER_LAW! use power law, instead of the thermally activated flow
      logical :: thermalSoftening_enabled ! decrease the intrinsic slip resistance in the power law flow rule, for thermal softening
      logical :: adiabaticHeating_enabled ! plastic dissipation locally heats up the material (valid at very high strain rates only -- check the heat diffusion time scale)
      integer :: backStressDamageMode  ! introduce backstress-induced cyclic softening at load reversals. 0: disabled. 1: reset to 0 if load reversed. 2: decay faster if load reversed.
      real(8) :: backStressDecayTime   ! characteristic decay time for back-stress parameter c
      logical :: damage_enabled  ! introduce damage (softening) to slip sytstems g, after exhausting plastic deformation lfie
      logical :: damage_gsat_enabled  ! introduce damage (softening) to slip sytstems g_sat, as voids form
      logical :: damage_nonlocal  ! regularize damage variable
      real(8) :: wt_BCC
      integer :: modelYieldPointPhenomenon  ! model used for introducing yield point phenomenon (0=off, 1=Ghosh 1998, 2=Shade 2017)
      
      ! applying pressure/creep loads
      logical :: pressureDirectionZ

      ! PROGRAM and POST-PROCESSING options.
      ! (values are imported from options.inp)
      ! -----------------------------------------
      ! error analysis / nodal value calculation
      ! -----------------------------------------
      integer :: posQuantity,nQuantDOF ! position of the quantity of interest, in the state var. vector
      character(len=maxStringLength) :: errInd_method ! identifies the method used in calculating the error indicator: jump, ZZ, norm_of_gradient
      integer :: p_deg
      integer :: recoveryMethod        ! recovery method to use to calc nodal Fp. 1: nodalAVG, 2: SPR. 3: LeastSQ
      logical :: normalize_SPR         ! normalize patch position to improve the scaling of the SPR matrix
      logical :: regularizeBoundaries  ! whether to do averaging of the SPR-extrapolated values at boundary nodes
      integer :: localNonlocal         ! 1: use local SPR-patches (formed by elements directly connected to nodes), 2: use distant elements to form patches

      contains
            
      SUBROUTINE readOptions(strOptionsFile)
      
      implicit none
      
      character(len=*), intent(in) :: strOptionsFile
      ! locals
      integer :: error
      character(len=maxStringLength) :: strOptionName
      character(len=maxStringLength) :: strOptionValue
      character(len=300) :: strLine
      real(8) :: floBuffer
      integer :: wPos(10,2), wCount, I
      
      CALL options_SetDefaults()

      open(104,file=trim(strOptionsFile))
      nOptions=0
      do while(.TRUE.)
         read(104,'(A)',iostat=error) strLine
         if (error.NE.0) exit             ! EOF / Error
         if (strLine(1:2).EQ.'//') cycle  ! comment. skip.
         
         ! extract the option name and value
         read(strLine,*,iostat=error) strOptionName, floBuffer
         if (error.NE.0) floBuffer=0      ! option value is not a number
         ! extract option value as string

         CALL parseString(strLine,' ',wPos(:,:),10,wCount)
         if(wCount.NE.2) cycle            ! incomplete input. skip.
         
         nOptions=nOptions+1
         option_names(nOptions) = toLowerCase(strOptionName)
         options_Float(nOptions) = floBuffer
         options_Int(nOptions) = floBuffer
         options_Bool(nOptions) = (options_Int(nOptions).NE.0)
         options_String(nOptions) = strLine(wPos(2,1):wPos(2,2))
      enddo
      close(104)

      if(nOptions.GT.maxNumOptions) then
         write(*,*) 'INCREASE maxNumOptions to ',nOptions
         STOP
      endif
      
      do I=1,nOptions
         if (trim(option_names(I)).EQ.toLowerCase('postProcT')) then
            postProcT = options_Float(I)
            postProcN = 0
         elseif (trim(option_names(I)).EQ.toLowerCase('postProcN')) then
            postProcN = options_Int(I)
            postProcT = 0.0
         elseif (trim(option_names(I)).EQ.toLowerCase('dataExportT')) then
            dataExportT = options_Float(I)
            dataExportN = 0
         elseif (trim(option_names(I)).EQ.toLowerCase('dataExportN')) then
            dataExportN = options_Int(I)            
            dataExportT = 0.0
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_general')) then
            DEBUG_general = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_UMAT')) then
            DEBUG_UMAT = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_timer')) then
            DEBUG_timer = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_RESID')) then
            DEBUG_RESID = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_linesearch')) then
            DEBUG_linesearch = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_grains')) then
            DEBUG_grains = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_errInd')) then
            DEBUG_errInd = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('do_errAnalysis')) then
            do_errAnalysis = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('exportSolution')) then
            exportSolution = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('continueSolution')) then
            continueSolution = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('verboseUmat')) then
            verboseUmat = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('conventionDeka')) then
            conventionDeka = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_watmus')) then
            EXPORT_watmus = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_stiffnessMatrix')) then
            EXPORT_stiffnessMatrix = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_delta_gamma')) then
            EXPORT_delta_gamma = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_backstress')) then
            EXPORT_backstress = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_wp')) then
            EXPORT_wp = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_temperature')) then
            EXPORT_temperature = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_damage')) then
            EXPORT_damage = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_fp')) then
            EXPORT_fp = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_creep')) then
            EXPORT_creep = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_creepRate')) then
            EXPORT_creepRate = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_strainRate')) then
            EXPORT_strainRate = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_nye')) then
            EXPORT_nye = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_nodePos')) then
            EXPORT_nodePos = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_cauchy')) then
            EXPORT_cauchy = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_grainAVG')) then
            EXPORT_grainAVG = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_grainRot')) then
            EXPORT_grainRot = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_hardening')) then
            EXPORT_hardening = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('suppress_specialPoints')) then
            suppress_specialPoints = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_specialOnly')) then
            EXPORT_specialOnly = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('EXPORT_allAtLast')) then
            EXPORT_allAtLast = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('DEBUG_nodalValues')) then
            DEBUG_nodalValues = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('posQuantity')) then
            posQuantity = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('nQuantDOF')) then
            nQuantDOF = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('errInd_method')) then
            errInd_method = options_String(I)
            ! 'norm_of_gradient'
            ! 'jump'
            ! 'ZZ'
         elseif (trim(option_names(I)).EQ.toLowerCase('recoveryMethod')) then
            ! 1: nodal AVG
            ! 2: SPR
            ! 3: least SQ
            recoveryMethod = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('regularizeBoundaries')) then
            regularizeBoundaries = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('localNonlocal')) then
            localNonlocal = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('p_deg')) then
            p_deg = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('normalize_SPR')) then
            normalize_SPR = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('materialName')) then
            materialName = options_String(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('elastic_only')) then
            ELASTIC_ONLY = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('tensOnly')) then
            tensOnly = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('methodTensCompCriterion')) then
            methodTensCompCriterion = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('smoothTCtransition')) then
            smoothTCtransition = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('softOnly')) then
            softOnly = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('prismatic_only')) then
            PRISMATIC_ONLY = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('basal_only')) then
            BASAL_ONLY = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('hardening_disabled')) then
            HARDENING_DISABLED = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('backstress_disabled')) then
            BACKSTRESS_DISABLED = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('hallpetch_disabled')) then
            HALLPETCH_DISABLED = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('GND_disabled')) then
            GND_DISABLED = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('use_power_law')) then
            USE_POWER_LAW = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('thermalSoftening_enabled')) then
            thermalSoftening_enabled = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('adiabaticHeating_enabled')) then
            adiabaticHeating_enabled = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('backStressDamageMode')) then
            backStressDamageMode = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('backStressDecayTime')) then
            backStressDecayTime = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('damage_enabled')) then
            damage_enabled = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('damage_gsat_enabled')) then
            damage_gsat_enabled = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('damage_nonlocal')) then
            damage_nonlocal = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('wt_BCC')) then
            wt_BCC = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('pressureDirectionZ')) then
            pressureDirectionZ = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('UMAT_lsch')) then
            UMAT_lsch = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('UMAT_halfStep')) then
            UMAT_halfStep = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('UMAT_stabilize')) then
            UMAT_stabilize = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('TAU_LIMIT')) then
            TAU_LIMIT = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('improve_U0')) then
            improve_U0 = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('STEP_v')) then
            STEP_v = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('WLF_c2')) then
            WLF_c2 = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('lschMethod')) then
            lschMethod = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('STEP_max')) then
            STEP_max = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('solverMethod')) then
            solverMethod = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('useTangentStiffnessStrainUpdates')) then
            useTangentStiffnessStrainUpdates = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('useLogStrainRate')) then
            useLogStrainRate = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('endSimulationAtWp')) then
            endSimulationAtWp = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('StiffnessMatrixFormulation')) then
            StiffnessMatrixFormulation = options_Int(I)
            if (StiffnessMatrixFormulation == SM_delPdelF) SymmetrizeStiffnessMatrix = .false.
            if (StiffnessMatrixFormulation == SM_delSdelLnU) SymmetrizeStiffnessMatrix = .true.
         elseif (trim(option_names(I)).EQ.toLowerCase('adaptive_time_stepping')) then
            adaptive_time_stepping = options_Bool(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('convergenceNorm')) then
            convergenceNorm = options_Int(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('eps_F')) then
            eps_F = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('eps_u')) then
            eps_u = options_Float(I)
         elseif (trim(option_names(I)).EQ.toLowerCase('modelYieldPointPhenomenon')) then   
            modelYieldPointPhenomenon = options_Int(I)
         endif
      enddo
      
      END SUBROUTINE readOptions
      
      SUBROUTINE options_SetDefaults()
         implicit none
               ! initialize to default values
         postProcT = 0.0
         postProcN = 2
         dataExportT = 0.0
         dataExportN = 2
         DEBUG_general = .FALSE.
         DEBUG_UMAT = .FALSE.
         DEBUG_timer = .FALSE.
         DEBUG_grains = .FALSE.
         DEBUG_errInd = .FALSE.
         DEBUG_RESID = .FALSE.
         DEBUG_linesearch = .FALSE.
         DEBUG_nodalValues=.FALSE.
         do_errAnalysis = .FALSE.
         exportSolution = .FALSE.
         continueSolution = .FALSE.
         verboseUmat      = .FALSE.
         conventionDeka   = .FALSE.
         EXPORT_watmus = .FALSE.
         EXPORT_stiffnessMatrix = .false.
         EXPORT_delta_gamma = .FALSE.
         EXPORT_backstress = .FALSE.
         EXPORT_wp = .FALSE.
         EXPORT_temperature = .FALSE.
         EXPORT_damage = .FALSE.
         EXPORT_fp = .FALSE.
         EXPORT_creep = .FALSE.
         EXPORT_creepRate = .FALSE.
         EXPORT_strainRate = .FALSE.
         EXPORT_nye = .FALSE.
         EXPORT_nodePos = .FALSE.
         EXPORT_cauchy = .TRUE.
         EXPORT_grainAVG = .TRUE.
         EXPORT_grainRot = .FALSE.
         EXPORT_hardening = .FALSE.
         EXPORT_specialOnly = .FALSE.
         suppress_specialPoints = .FALSE.
         EXPORT_allAtLast = 1.D0
         posQuantity = 211
         nQuantDOF = 6
         p_deg=1
         recoveryMethod = 3   ! least Square Fit
         normalize_SPR = .FALSE.
         regularizeBoundaries = .FALSE.
         localNonlocal = 1
         errInd_method = 'none'
         materialName = 'Ti6242'
         ELASTIC_ONLY = .FALSE.
         tensOnly = 0
         methodTensCompCriterion = 4   ! use Lode angle (clockwise from principal stress axis S3)
         smoothTCtransition = .false.
         softOnly = 0
         HARDENING_DISABLED = .FALSE. 
         BACKSTRESS_DISABLED = .FALSE.
         HALLPETCH_DISABLED = .FALSE.
         PRISMATIC_ONLY = .FALSE.
         BASAL_ONLY = .FALSE.
         GND_DISABLED = .FALSE.
         USE_POWER_LAW = .FALSE.
         thermalSoftening_enabled = .FALSE.
         adiabaticHeating_enabled = .FALSE.
         backStressDamageMode = 0
         backStressDecayTime = 10.d0
         damage_enabled = .FALSE.
         damage_gsat_enabled = .FALSE.
         damage_nonlocal = .FALSE.
         wt_BCC = 0.12d0   ! % volume fraction of BCC phase in trans-B grains
         pressureDirectionZ = .FALSE.
         UMAT_lsch = .FALSE.
         UMAT_halfStep = .TRUE.
         UMAT_stabilize = .TRUE.
         TAU_LIMIT = 1.5d0
         improve_U0 = 4
         STEP_v = 1.0
         WLF_c2 = 0.9
         STEP_max = 10.0
         lschMethod = 0
         solverMethod = solmeth_BFGSkeep         ! BFGS, keep K
         useTangentStiffnessStrainUpdates = .false.
         useLogStrainRate = .true.
         endSimulationAtWp = -1.d0
         StiffnessMatrixFormulation = SM_delPdelF ! default with FBar method
         SymmetrizeStiffnessMatrix = .false.
         adaptive_time_stepping = .TRUE.
         convergenceNorm = 1
         eps_F = 0.007    ! convergence tolerance/precision on force and displacement
         eps_u = 0.028   
         modelYieldPointPhenomenon = 0 ! off
         
      END SUBROUTINE options_SetDefaults
      
      SUBROUTINE summarizeOptions()
      
         if (EXPORT_stiffnessMatrix) write(*,*) 'will export stiffness matrix and quit'
      
         if (postProcT.NE.0.D0) write(*,*)  & 
               'post processing time period:',postProcT
         if (postProcN.NE.0) write(*,*)  & 
               'post processing increment period:',postProcN
         if (dataExportT.NE.0.D0) write(*,*)  & 
               'data export time period:',dataExportT
         if (dataExportN.NE.0) write(*,*)  & 
               'data export increment period:',dataExportN
         
         write(*,*) 'using material:',materialName
         if(ELASTIC_ONLY)  & 
            write(*,*) 'plasticity disabled'
         if(GND_DISABLED) & 
            write(*,*) 'GND hardening disabled'
         if(.not.USE_POWER_LAW) & 
            write(*,*) 'Thermally Activated Flow rule instead of power law'
         if(thermalSoftening_enabled) & 
            write(*,*) 'thermal softening is enabled !'
         if(adiabaticHeating_enabled) & 
            write(*,*) 'adiabatic heating is enabled !'
         if(backStressDamageMode /= 0) & 
            write(*,*) 'load reversals induce back-stress loss !'
         if(recoveryMethod.EQ.1) then
            write(*,*) 'Recovery method: nodal averaging'
         elseif(recoveryMethod.EQ.2) then
            write(*,*) 'Recovery method: SPR'
         elseif(recoveryMethod.EQ.3) then
            write(*,*) 'Recovery method: least square fit'
         endif
         
         if(do_errAnalysis) then
            write(*,*) 'error analysis ON'
         endif

         if(improve_U0.EQ.0) then
            write(*,*) & 
       'no improvement for initial displacement field'
         elseif(improve_U0.EQ.1) then
            write(*,*) & 
       'improve initial displacement fields using: Jiahao''s method'
         elseif(improve_U0.EQ.3) then
            write(*,*) & 
       'improve initial displacement fields using: (quadratic) Jiahao''s method'
         elseif(improve_U0.EQ.2) then
            write(*,*) & 
       'improve initial displacement fields using: Deniz''s method'
         elseif(improve_U0.EQ.4) then
            write(*,*) & 
       'improve initial displacement fields using most recent tangent matrix'
         else
            write(*,*) & 
       'improve initial guess - method not recognised:',improve_U0
         endif
         
         if (StiffnessMatrixFormulation /= SM_delPdelF) &
            write(*,*) 'WARNING ! old tangent matrix formulation is used:',StiffnessMatrixFormulation            
         if (SymmetrizeStiffnessMatrix) &
            write(*,*) 'WARNING ! symmetrizing stiffness matrix'

         if(lschMethod.EQ.0) write(*,*) & 
       'line search disabled.'
         if(lschMethod.EQ.1) write(*,*) & 
       'line search enabled. will use method: Deniz'
         if(lschMethod.EQ.2) write(*,*) & 
       'line search enabled. will use method: Matthies, Strang, 1979'
       
         if (conventionDeka) then
            write(*,*) 'WARNING ! Deka''s convention of euler angles are used'
         else
            write(*,*) 'WARNING ! ''Active'' convention of euler angles are used'
         endif
         
         if (backStressDamageMode /= 0) then
            write(*,*) 'back Stress cyclic damage enabled:',backStressDamageMode
            write(*,*) 'B.S. decay time:',backStressDecayTime
         endif

         write(*,*) '# options imported:',nOptions
         
      END SUBROUTINE
            
      function toLowerCase(strIn) result(strOut)
      ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
      ! Original author: Clive Page

      implicit none

      character(len=*), intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j

      do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
               strOut(i:i) = strIn(i:i)
          end if
      end do

      end function toLowerCase

      ! Author: Deniz Ozturk (2014)
      ! parses a string into words separated by character 'del'
      SUBROUTINE parseString(str, del, wordPos, arrSize, wordCount)
      implicit none
      character(len=*), intent(in) :: str
      character, intent(in) :: del
      integer, intent(in) :: arrSize
      integer, intent(out) :: wordPos(arrSize,2)
      integer, intent(out) :: wordCount
      
      integer :: i,wordCountMax
      character :: prevChar, thisChar, del2
      
      ! space and tab alternatively
      del2=del
      if(del.EQ.' ') del2=ACHAR(9)
      if(IACHAR(del).EQ.9) del2=' '
      
      wordCountMax = wordCount
      wordCount = 0
      wordPos(:,:) = 0
      prevChar = del
      do i=1,len(str)
         thisChar = str(i:i)
         if((prevChar.NE.del.AND.prevChar.NE.del2)   &
       .AND.(thisChar.EQ.del.OR.thisChar.EQ.del2)) then
            ! end of a word
            wordPos(wordCount,2) = i - 1
         elseif((prevChar.EQ.del.OR.prevChar.EQ.del2) &
       .AND.(thisChar.NE.del.AND.thisChar.NE.del2))  then
            ! start of a word
            wordCount = wordCount + 1
            wordPos(wordCount,1) = i
         endif
         prevChar = thisChar
      enddo
      if(prevChar.NE.del.AND.prevChar.NE.del2) then
         !end of string ends the word
         wordPos(wordCount,2) = len(str)
      endif
      END SUBROUTINE

      END MODULE options
      