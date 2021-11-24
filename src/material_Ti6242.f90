      module material_Ti6242
      
      use crystal_Titanium    ! slip system definitions
      
      implicit none
      
      integer, parameter :: mat_nstate = 377
      
      ! unit conversion
      real(8), parameter :: nm_to_um = 1.d-3 ! nanometer to micrometer unit conversion for the burgers vector
      real(8), parameter :: Jm_to_MPa = 1.d-6 ! J/m^3 to MPa unit conversion
      
      ! MATERIAL PARAMETERS
      !-------------------------
      ! variables below store parameters for tension/compression and pri-A/trans-B:
      ! 1st index: slip system, or slip family (..._fam)
      ! 2nd index: 1: tension. 2: compression
      ! 3rd index: Phase: 1: HCP in primary Alpha. 2: HCP in transformed Beta
      
      ! HCP - hardening related params
      ! - per family
      real(8) :: g0_bar_famHCP(nslipfam,2,2) ! internal slip resistance
      real(8) :: xh0_famHCP(nslipfam,2,2)    ! max (initial) hardening rates
      real(8) :: xs_bar_famHCP(nslipfam,2,2) ! max (reference) limit slip resistances
      real(8) :: xn_famHCP(nslipfam,2,2)     ! sensitivity of limit strength to shear rate
      real(8) :: xr_famHCP(nslipfam,2,2)     ! sensitivity of hardening rate to slip resistance
      ! - per slip system
      real(8) :: g0_bar_HCP(nslip_max,2,2)    ! internal slip resistances
      real(8) :: xh0_HCP(nslip_max,2,2)
      real(8) :: xs_bar_HCP(nslip_max,2,2)
      real(8) :: xn_HCP(nslip_max,2,2)        
      real(8) :: xr_HCP(nslip_max,2,2)        

      ! BCC - hardening related parameters
      ! - per family
      real(8) :: g0_bar_famBCC(nslipfam,2)     ! internal slip resistance:
                                               ! 1st index: slip families: {110},{112}soft,{112}hard,{123}soft,{123}hard
                                               ! 2nd index: 1: b=[-111] systems, 2: all others
      real(8) :: xh_0_famBCC(nslipfam)         ! BCC hardening params per slip family
      real(8) :: xh_inf_famBCC(nslipfam)       ! BCC hardening params per slip family
      real(8) :: xg_0_famBCC(nslipfam)         ! BCC hardening params per slip family
      real(8) :: xg_inf_famBCC(nslipfam)       ! BCC hardening params per slip family
      ! - per slip system.
      ! 1st index: slip system
      ! 2nd index: 1: soft 2: hard direction
      real(8) :: g0_bar_BCC(nslip_max,2)    ! internal slip resistance (for hard/soft directions)
      real(8) :: xh_0_BCC(nslip_max,2)      ! BCC hardening params per slip family
      real(8) :: xh_inf_BCC(nslip_max,2)    ! BCC hardening params per slip family
      real(8) :: xg_0_BCC(nslip_max,2)      ! BCC hardening params per slip family
      real(8) :: xg_inf_BCC(nslip_max,2)    ! BCC hardening params per slip family
      
      ! flow related parameters
      real(8) :: burgers_famHCP(nslipfam)     ! HCP - burgers vector lengths
      real(8) :: burgers_HCP(nslip_max)       
      real(8) :: burgers_famBCC(nslipfam)     ! BCC - burgers vector lengths
      real(8) :: burgers_BCC(nslip_max,2)
      real(8) :: coeffGND                    ! calibration coefficient for GND hardening
      ! HCP - per family
      real(8) :: a_dot_famHCP(nslipfam)      ! HCP reference slip rate (power law formulation)
      real(8) :: a_dot_therm_famHCP(nslipfam)! HCP reference slip rate (thermally activated flow)
      real(8) :: rate_exp_famHCP(nslipfam)   ! HCP flow sensitivity    (power law formulation)
      real(8) :: E_act_famHCP(nslipfam)      ! HCP slip activation energy (thermally activated flow)
      ! HCP - per slip system
      real(8) :: a_dot_HCP(nslip_max)         ! reference slip rate (power law formulation)
      real(8) :: a_dot_therm_HCP(nslip_max)   ! reference slip rate (thermally activated flow)
      real(8) :: rate_exp_HCP(nslip_max)      ! flow sensitivity (power law formulation)
      real(8) :: E_act_HCP(nslip_max)         ! Slip activation energy at each slip system
      real(8) :: softening_T_char_HCP(nslip_max)  ! characteristic temperature for drag-dominated flow (full thermal activation)
      real(8) :: softening_g_therm_HCP(nslip_max) ! thermally activated (local) stress barrier to flow
      real(8) :: g_drop_c_a_300_400               ! drop in the internal slip resistance of the c+a systems rom 300K to 400K (MPa)
      
      ! BCC - per family
      real(8) :: a_dot_famBCC(nslipfam)             ! BCC reference slip rate  (power law formulation)
      real(8) :: a_dot_therm_famBCC(nslipfam)       ! BCC reference slip rate  (thermally activated flow)
      real(8) :: rate_exp_famBCC(nslipfam)          ! BCC flow sensitivity     (power law formulation)
      real(8) :: E_act_famBCC(nslipfam)             ! BCC slip activation energy (thermally activated flow)
      real(8) :: softening_T_char_BCC(nslip_max)    ! characteristic temperature for drag-dominated flow (full thermal activation)
      real(8) :: softening_g_therm_BCC(nslip_max)   ! thermally activated (local) stress barrier to flow
      ! BCC - per slip system
      ! 2nd index: 1: soft direction (-), 2: hard direction (+)
      real(8) :: a_dot_BCC(nslip_max,2)       ! reference slip rate (soft/hard)(power law formulation)
      real(8) :: a_dot_therm_BCC(nslip_max,2) ! reference slip rate (soft/hard)(thermally activated flow)
      real(8) :: rate_exp_BCC(nslip_max,2)    ! flow sensitivity (soft/hard)(power law formulation)
      real(8) :: E_act_BCC(nslip_max,2)       ! Slip activation energy at each slip system (soft/hard)
      
      ! Hall-Petch parameters
      ! if, D_hp < 0, grain size is used for that slip system
      real(8) :: K_hp_HCP(nslipfam,2,2)   ! Hall Petch coefficient HCP tens/Comp, pri-A/trans-B
      real(8) :: D_hp_HCP(nslipfam,2)     ! charactrs. slip length HCP pri-A/trans-B
      real(8) :: K_hp_BCC(2)              ! Hall Petch coefficient BCC
      real(8) :: D_hp_BCC(2)              ! charactrs. slip length BCC
      
      integer, parameter :: STATE_TENSION=1, STATE_COMPRESSION=2
      ! characteristic stress for a smooth transition of the T/C state
      real(8) :: characteristicStressTCtransition
      
      ! Material Internal Variables - to be passed down to flow rule
      ! (generic, for HCP/BCC)
      real(8) :: a_dot(nslip_max)
      real(8) :: a_dot_star(nslip_max)
      real(8) :: a_dot_therm(nslip_max)
      real(8) :: rate_exp(nslip_max)
      real(8) :: E_act(nslip_max)
      real(8) :: burgers(nslip_max)
      ! hardening variables
      real(8) :: g0_BCC(nslip_max,2)      ! BCC, initial slip resistance (g0_bar + g_hp)
      real(8) :: g_BCC(nslip_max,2)       ! BCC, current slip resistance (g0_bar + g_hp + g_hardening)
      real(8) :: g_bar_HCP(nslip_max,2)   ! HCP, initial slip resistance (g0_bar + g_hp)
      real(8) :: g0_HCP(nslip_max,2)      ! HCP, initial slip resistance (g0_bar + g_hp)
      real(8) :: g_HCP(nslip_max,2)       ! HCP, current slip resistance (g0_bar + g_hp + g_hardening)
      real(8) :: dd_g(nslip_max)          ! change in work-hardening over this time step (tension/compression independent, per slip system)
      real(8) :: delta_g_BCC(nslip_max)   ! work-hardening part of slip resistance (slip systems)   (read from statev)
      real(8) :: delta_g_HCP(nslip_max)   ! work-hardening part of slip resistance (slip systems)   (read from statev)
      real(8) :: chi_t(nslip_max)         ! back stress at time t
      real(8) :: chi_tau(nslip_max)       ! back stress at time t
      real(8) :: lambda_beta_t(nslip_max) ! lattice incompatibilities on crystal planes
      
!      integer :: iCrystal
!      integer :: iPhase
      
      real(8) :: xc  ! back stress evolution: kinematic hardening coefficient
      real(8) :: xd  ! back stress evolution: dynamic recovery coefficient
      
      real(8) :: hardeningSelf         ! self-hardening coefficient q(i,i)
      real(8) :: hardeningLatent       ! latent-hardening coefficient q(i,j)
      
      real(8) :: p_therm, q_therm      ! flow exponents p and q for the thermally activated flow rule
      
      real(8) :: heatCapacityVolumetric! volumetric heat capacity (J/Kelvins/m^3)
      
      real(8) :: c11_HCP,c12_HCP,c13_HCP,c33_HCP,c55_HCP ! HCP elasticity coeffs
      real(8) :: c1_BCC,c2_BCC,c3_BCC                    ! BCC elasticity coeffs

      real(8) :: c11_k_T_HCP,c12_k_T_HCP,c13_k_T_HCP, &  ! rate of change of HCP elastic coefficients with temperature
                 c33_k_T_HCP,c55_k_T_HCP 
      real(8) :: c1_k_T_BCC,c2_k_T_BCC,c3_k_T_BCC        ! rate of change of BCC elastic coefficients with temperature

      real(8) :: thrm_C(3,2)   ! thermal expansion coefficients along crystal axes 1,2,3 (HCP/BCC)
      
      real(8) :: matrixB_HCP(63,9) ! the B-matrix, used in calculating GND densities from the Nye tensor
      real(8) :: screwGND(nslip_max)  ! density of screw-type GNDs in each slip sys
      real(8) :: edgeNGND(nslip_max)  ! density of edge-N-type GNDs in each slip sys
      real(8) :: edgeTGND(nslip_max)  ! density of edge-T-type GNDs in each slip sys
      real(8) :: forestGND(nslip_max) ! forest GND density for each slip system
      real(8) :: parallelGND(nslip_max)! parallel GND density for each slip system
      real(8) :: stressCut(nslip_max) ! cutting stress for each slip system
      real(8) :: stressPass(nslip_max)! cutting stress for each slip system
      real(8) :: K_star,lp,a_dot_init ! yield point phenomenon properties (model 1)
      real(8) :: ep0,gYPP_basal,gYPP_prism,gYPP_ca ! yield point phenomenon properties (model 2)
      
      contains
      
      SUBROUTINE initialize_Parameters(props,nprops,matName)
            
      implicit none
      
      real(8), intent(in) :: props(nprops)
      integer, intent(in) :: nprops
      character(*), intent(in) :: matName
      
      integer :: iTC
      
      ! initialize crystal slip systems
      call init_crys_HCP()
      call init_crys_BCC()
      
      xc=props(137)  ! back-stress accummulation coefficient
      xd=props(138)  ! back-stress dynamic recovery  coefficient

      hardeningSelf = props(20)   ! self and latent hardening coefficients for the phenomenological hardening rule
      hardeningLatent = props(21)
      
      p_therm = props(22)  ! p and q - rate exponents for the thermally activated flow rule
      q_therm = props(23)
      
      if (hardeningSelf == 0.d0 .OR. hardeningLatent == 0.d0) then
         write(*,*) 'WARNING !!! - set hardeningSelf, hardeningLatent in props.dat'
         stop
      endif
            
      heatCapacityVolumetric = props(24)     ! volumetric heat capacity (megaJ/Kelvins/m^3) -- pure Titanium: 2.43 (megaJ/K/m^3=MPa/K) --> wp(MPa=MJ/m^3) / VHC (MJ/K/m^3) = dT (K)
      
      ! read in elasticity coefficients
      ! -------------------------------
      ! HCP
      c11_HCP=props(1)
      c12_HCP=props(2)
      c13_HCP=props(3)
      c33_HCP=props(4)
      c55_HCP=props(5)
      ! BCC
      c1_BCC=props(297)
      c2_BCC=props(298)
      c3_BCC=props(299)
      
      ! read in thermal softening rate of elastic constants
      ! HCP
      c11_k_T_HCP=props(139)
      c12_k_T_HCP=props(140)
      c13_k_T_HCP=props(141)
      c33_k_T_HCP=props(142)
      c55_k_T_HCP=props(143)
      ! BCC
      c1_k_T_BCC=props(365)
      c2_k_T_BCC=props(366)
      c3_k_T_BCC=props(367)
      
      
      ! read in thermal expansion coefficients         
      ! -------------------------------
      ! HCP - transverse isotropy
      thrm_C(1,HCP) = props(6)    
      thrm_C(2,HCP) = props(6)
      thrm_C(3,HCP) = props(7)         
      ! BCC - isotropic elasticity
      thrm_C(1:3,BCC) = props(300)
      
      coeffGND = props(17)
      
      ! read in flow rule parameters
      ! -------------------------------      
      ! HCP
      burgers_famHCP(1:8)    = props(9:16)
      a_dot_famHCP(1:8)      = props(25:32)  
      rate_exp_famHCP(1:8)   = props(33:40)
      a_dot_therm_famHCP(1:8)= props(41:48)
      E_act_famHCP(1:8)      = props(49:56)
      ! distribute to slip systems
      call slipFamily_to_SystemsHCP(burgers_famHCP,    burgers_HCP,    nslips_famHCP)      
      call slipFamily_to_SystemsHCP(a_dot_famHCP,      a_dot_HCP,      nslips_famHCP)
      call slipFamily_to_SystemsHCP(a_dot_therm_famHCP,a_dot_therm_HCP,nslips_famHCP)
      call slipFamily_to_SystemsHCP(rate_exp_famHCP,   rate_exp_HCP,   nslips_famHCP)
      call slipFamily_to_SystemsHCP(E_act_famHCP,      E_act_HCP,      nslips_famHCP)
      
      ! thermal softening law:
      ! delta_g(dT) = -[1-exp(dT/T_c)]*Dg_therm where dT = T - 300K
      ! HCP - T_c
      softening_T_char_HCP(1:3)   = props(41)   ! 3 x basal a
      softening_T_char_HCP(4:6)   = props(42)   ! 3 x prism a
      softening_T_char_HCP(7:12)  = props(43)   ! 6 x pyramidal a
      softening_T_char_HCP(13:30) = props(44)   !18 x 1st & 2nd pyramidal c+a
      ! HCP - Dg_therm
      softening_g_therm_HCP(1:3)   = props(45)  ! 3 x basal a
      softening_g_therm_HCP(4:6)   = props(46)  ! 3 x prism a
      softening_g_therm_HCP(7:12)  = props(47)  ! 6 x pyramidal a
      softening_g_therm_HCP(13:30) = props(48)  !18 x 1st & 2nd pyramidal c+a
      ! drop in the internal slip resistance of the c+a systems rom 300K to 400K (MPa)
      g_drop_c_a_300_400 = props(8)
      
      ! BCC
      burgers_famBCC(1:5)    = props(301:305)
      a_dot_famBCC(1:5)      = props(313:317)
      rate_exp_famBCC(1:5)   = props(318:322)
      a_dot_therm_famBCC(1:5)= props(323:327)
      E_act_famBCC(1:5)      = props(328:332)
      ! distribute to slip systems
      call slipFamily_to_SystemsBCC(burgers_famBCC,    burgers_BCC,    nslips_famBCC)
      call slipFamily_to_SystemsBCC(a_dot_famBCC,      a_dot_BCC,      nslips_famBCC)
      call slipFamily_to_SystemsBCC(a_dot_therm_famBCC,a_dot_therm_BCC,nslips_famBCC)
      call slipFamily_to_SystemsBCC(rate_exp_famBCC,   rate_exp_BCC,   nslips_famBCC)
      call slipFamily_to_SystemsBCC(E_act_famBCC,      E_act_BCC,      nslips_famBCC)
      
      ! thermal softening law:
      ! delta_g(dT) = -[1-exp(dT/T_c)]*Dg_therm where dT = T - 300K
      ! BCC - T_c
      softening_T_char_BCC(1:48)  = props(41)   ! use the same for 3 x basal a
      ! BCC - Dg_therm
      softening_g_therm_BCC(1:48)   = props(41) ! use the same for 3 x basal a
      
      ! read in hardening parameters for HCP in primary-alpha grains
      ! ------------------------------------------------------------
      g0_bar_famHCP(1:8,1,pri_A) = props(57:64)   ! internal slip resistances (excluding HallPetch) for HCP/primary-a in tension    
      g0_bar_famHCP(1:8,2,pri_A) = props(65:72)   ! internal slip resistances (excluding HallPetch) for HCP/primary-a in compression    
      xn_famHCP(1:8,1,pri_A) = props(73:80)  ! limit CRSS strain-rate response exponent // HCP in primary-a // tension 
      xn_famHCP(1:8,2,pri_A) = props(81:88)  ! limit CRSS strain-rate response exponent // HCP in primary-a // compression 
      xr_famHCP(1:8,1,pri_A) = props(89:96)  ! CRSS saturation exponent // HCP in primary-a // tension 
      xr_famHCP(1:8,2,pri_A) = props(97:104) ! CRSS saturation exponent // HCP in primary-a // compression 
      xh0_famHCP(1:8,1,pri_A) = props(105:112) ! h0: initial hardening rate // HCP in primary-a // tension
      xh0_famHCP(1:8,2,pri_A) = props(113:120) ! h0: initial hardening rate // HCP in primary-a // compression
      xs_bar_famHCP(1:8,1,pri_A)=props(121:128) ! HCP in primary-a // tension
      xs_bar_famHCP(1:8,2,pri_A)=props(129:136) ! HCP in primary-a // compression
      ! distribute to slip systems
      do iTC = 1,2   ! tension/compression
         call slipFamily_to_SystemsHCP(g0_bar_famHCP(:,iTC,pri_A),g0_bar_HCP(:,iTC,pri_A),nslips_famHCP)
         call slipFamily_to_SystemsHCP(xn_famHCP(:,iTC,pri_A),    xn_HCP(:,iTC,pri_A),    nslips_famHCP)
         call slipFamily_to_SystemsHCP(xr_famHCP(:,iTC,pri_A),    xr_HCP(:,iTC,pri_A),    nslips_famHCP)
         call slipFamily_to_SystemsHCP(xh0_famHCP(:,iTC,pri_A),   xh0_HCP(:,iTC,pri_A),   nslips_famHCP)
         call slipFamily_to_SystemsHCP(xs_bar_famHCP(:,iTC,pri_A),xs_bar_HCP(:,iTC,pri_A),nslips_famHCP)
      enddo
      
      ! read in hardening parameters for HCP in transormed-Beta grains
      ! --------------------------------------------------------------
      g0_bar_famHCP(1:8,1,trans_B) = props(193:200)   ! internal slip resistances (excluding HallPetch) // HCP in trans-B // tension    
      g0_bar_famHCP(1:8,2,trans_B) = props(201:208)   ! internal slip resistances (excluding HallPetch) // HCP in trans-B // compression    
      xn_famHCP(1:8,1,trans_B) = props(209:216)  ! n: limit CRSS strain-rate response exponent // HCP in trans-B // tension 
      xn_famHCP(1:8,2,trans_B) = props(217:224)  ! n: limit CRSS strain-rate response exponent // HCP in trans-B // compression 
      xr_famHCP(1:8,1,trans_B) = props(225:232)  ! r: CRSS saturation exponent // HCP in trans-Beta // tension 
      xr_famHCP(1:8,2,trans_B) = props(233:240)  ! r: CRSS saturation exponent // HCP in trans-Beta // compression 
      xh0_famHCP(1:8,1,trans_B) = props(241:248)    ! h0: initial hardening rate // HCP in trans-Beta // tension
      xh0_famHCP(1:8,2,trans_B) = props(249:256)    ! h0: initial hardening rate // HCP in trans-Beta // compression
      xs_bar_famHCP(1:8,1,trans_B) = props(257:264) ! HCP in trans-B // tension
      xs_bar_famHCP(1:8,2,trans_B) = props(265:272) ! HCP in trans-B // compression
      
      ! distribute to slip systems
      do iTC = 1,2   ! tension/compression
         call slipFamily_to_SystemsHCP(g0_bar_famHCP(:,iTC,trans_B),g0_bar_HCP(:,iTC,trans_B),nslips_famHCP)
         call slipFamily_to_SystemsHCP(xn_famHCP(:,iTC,trans_B),    xn_HCP(:,iTC,trans_B),    nslips_famHCP)
         call slipFamily_to_SystemsHCP(xr_famHCP(:,iTC,trans_B),    xr_HCP(:,iTC,trans_B),    nslips_famHCP)
         call slipFamily_to_SystemsHCP(xh0_famHCP(:,iTC,trans_B),   xh0_HCP(:,iTC,trans_B),   nslips_famHCP)
         call slipFamily_to_SystemsHCP(xs_bar_famHCP(:,iTC,trans_B),xs_bar_HCP(:,iTC,trans_B),nslips_famHCP)
      enddo
      
      ! BCC slip families (used in indexing of g0_fam,xh_0,xh_inf,xg_0,xg_inf)
      ! [1] = all n={101}, all b=<111> systems
      ! [2] = all n={112}, b=(soft) directions (-)
      ! [3] = all n={112}, b=(hard) directions (+)
      ! [4] = all n={123}, b=(soft) directions (-)
      ! [5] = all n={123}, b=(hard) directions (+)
      
      ! read in hardening parameters for BCC
      ! ------------------------------------------------------------
      g0_bar_famBCC(1:5,1)= props(333:337)   ! internal slip resistances (excluding HallPetch) // BCC in trans-B, b1[-111]
      g0_bar_famBCC(1:5,2)= props(338:342)   ! internal slip resistances (excluding HallPetch) // BCC in trans-B, all others
      xh_0_famBCC(1:5)  = props(345:349)    
      xh_inf_famBCC(1:5)= props(111:115)    
      xg_0_famBCC(1:5)  = props(355:359)    
      xg_inf_famBCC(1:5)= props(360:364)    
      ! distribute to slip systems
      call slipFamily_to_SystemsBCC(g0_bar_famBCC(:,2),g0_bar_BCC,nslips_famBCC)    ! distribute to all slip systems
      call slipFamily_to_SystemsBCC_b1(g0_bar_famBCC(:,1),g0_bar_BCC,nslips_famBCC) ! distribute only to b1=[-111] slip systems
      call slipFamily_to_SystemsBCC(xh_0_famBCC,  xh_0_BCC,  nslips_famBCC)
      call slipFamily_to_SystemsBCC(xh_inf_famBCC,xh_inf_BCC,nslips_famBCC)
      call slipFamily_to_SystemsBCC(xg_0_famBCC,  xg_0_BCC,  nslips_famBCC)
      call slipFamily_to_SystemsBCC(xg_inf_famBCC,xg_inf_BCC,nslips_famBCC)
      
      ! Hall-Petch parameters
      ! if, D_hp < 0, grain size is used for that slip system
      K_hp_HCP(1:8,1,pri_A)  = props(145:152)    ! Hall Petch coefficient HCP, tension, pri-A
      K_hp_HCP(1:8,2,pri_A)  = props(153:160)    ! Hall Petch coefficient HCP, compression, pri-A
      D_hp_HCP(1:8,pri_A)    = props(161:168)    ! charactrs. slip length HCP, pri-A
      K_hp_HCP(1:8,1,trans_B)= props(273:280)    ! Hall Petch coefficient HCP, tension, trans-B
      K_hp_HCP(1:8,2,trans_B)= props(281:288)    ! Hall Petch coefficient HCP, compression, trans-B
      D_hp_HCP(1:8,trans_B)  = props(289:296)    ! charactrs. slip length HCP, trans-B
      K_hp_BCC(1)            = props(369)        ! Hall Petch coefficient BCC, b1[-111]
      K_hp_BCC(2)            = props(370)        ! Hall Petch coefficient BCC, all others
      D_hp_BCC(1)            = props(371)        ! charactrs. slip length BCC, b1[-111]
      D_hp_BCC(2)            = props(372)        ! charactrs. slip length BCC, all others

      ! characteristic stress for a smooth transition of the T/C state
      characteristicStressTCtransition = props(144)
      
      ! parameters in yield point phenomenon
      ! model = 1
      K_star = props(373)
      lp = props(374)
      a_dot_init = props(375)
      ! model = 2
      ep0        = props(373)
      gYPP_basal = props(374)
      gYPP_prism = props(375)
      gYPP_ca    = props(376)
      
      ! initialize the B-matrix, used in calculating GND densities from the Nye tensor.
      CALL B_MATRIX_HCP(matrixB_HCP)
      
      

      END SUBROUTINE
      
      SUBROUTINE initialize_StateVariables(stateVars,nsvars,elemID,grainSizes,grainTexture,grainPhases,nGrains,grainIDelem,NELX)
      
      implicit none
      
      real(8), intent(out):: stateVars(nsvars)
      real(8), intent(in) :: grainSizes(nGrains)
      real(8), intent(in) :: grainTexture(3,nGrains)
      integer, intent(in) :: grainPhases(nGrains)
      integer, intent(in) :: elemID,nGrains,NELX,nsvars
      integer, intent(in) :: grainIDelem(NELX)

      
      integer :: grainID, iFam, phase, ps
      real(8) :: euler(3), ckrone(9), grainSize
      real(8) :: g_hp_fam(nslipfam),g_hp(nslip_max)
      real(8) :: g_hp_famBCC(nslipfam,2)
      
      ckrone = 0.d0
      ckrone(1) = 1.d0 ! ROW-MAJOR
      ckrone(5) = 1.d0 ! ROW-MAJOR
      ckrone(9) = 1.d0 ! ROW-MAJOR
   
      grainID = grainIDelem(elemID)
      grainSize = grainSizes(grainID)
      euler(1:3) = grainTexture(1:3,grainID)
      phase = grainPhases(grainID)
      
      stateVars(1 :30)  = 0.d0  ! delta_g, work hardening   HCP
      stateVars(31:60)  = 0.d0  ! tot_gamma                 HCP
      stateVars(61:66)  = 0.d0  ! 2nd PK   (Voigt)          HCP
      stateVars(67:75)  = ckrone! Fp       (Row Major)      HCP
      stateVars(76:123) = 0.d0  ! delta_g, work hardening   BCC
      stateVars(124:171)= 0.d0  ! tot_gamma                 BCC      
      stateVars(172:177)= 0.d0  ! 2nd PK   (Voigt)          BCC      
      stateVars(178:186)= ckrone! Fp       (Row Major)      BCC         
      stateVars(187:192)= 0.d0  ! Logarithmic Strain (Voigt)-generic, output only-
      stateVars(193)    = 0.d0  ! Local Temperature (ambient+adiabatic) -generic, output only-
      stateVars(211:216)= 0.d0  ! Cauchy Stress      (Voigt)-generic, output only-
      stateVars(221:250)= 0.d0  ! Chi, back stress          HCP
      stateVars(251:298)= 0.d0  ! Chi, back stress          BCC
      stateVars(299)       = 0.d0  ! w_p, plastic work         -generic-
      stateVars(300)       = 0.d0  ! volume                    -generic, used by FEM-UEL only-
      stateVars(301:303)= euler ! Crys.Orient.(Euler angles)-generic, output only-
      stateVars(304)       = phase ! phase (1:priA 2:traB)     -generic, input  only-
         
      ! determine hall-petch contributions to strength
      if (phase==pri_A) then
         ! HCP in pri-A, tension
         do iFam=1,nslipfam
            CALL calc_HallPetchStrength(K_hp_HCP(iFam,1,pri_A),D_hp_HCP(iFam,pri_A), &
                                        grainSize,g_hp_fam(iFam))
         enddo
         stateVars(305:312) = g_hp_fam(1:8)  ! HCP in primary-a, Tension
         ! HCP in pri-A, compression
         do iFam=1,nslipfam
            CALL calc_HallPetchStrength(K_hp_HCP(iFam,2,pri_A),D_hp_HCP(iFam,pri_A), &
                                        grainSize,g_hp_fam(iFam))
         enddo
         stateVars(313:320) = g_hp_fam(1:8)  ! HCP in primary-a, Compression
         
      elseif(phase==trans_B) then
         ! HCP in trans_B, tension
         do iFam=1,nslipfam
            CALL calc_HallPetchStrength(K_hp_HCP(iFam,1,trans_B),D_hp_HCP(iFam,trans_B), &
                                        grainSize,g_hp_fam(iFam))
         enddo
         stateVars(305:312) = g_hp_fam(1:8)  ! HCP in trans-B, Tension
         ! HCP in trans_B, compression
         do iFam=1,nslipfam
            CALL calc_HallPetchStrength(K_hp_HCP(iFam,2,trans_B),D_hp_HCP(iFam,trans_B), &
                                        grainSize,g_hp_fam(iFam))
         enddo
         stateVars(313:320) = g_hp_fam(1:8)  ! HCP in trans-B, Compression

         ! BCC in trans_B, only b=[-111] slip systems
         CALL calc_HallPetchStrength(K_hp_BCC(1),D_hp_BCC(1), &
                                     grainSize,g_hp_famBCC(1,1))
         g_hp_famBCC(2:5,1)=g_hp_famBCC(1,1)
         stateVars(337:341) = g_hp_famBCC(1:5,1)   ! hall-petch contributions for systems with b=[-111]
         
         ! BCC in trans_B, all other slip systems
         CALL calc_HallPetchStrength(K_hp_BCC(2),D_hp_BCC(2), &
                                     grainSize,g_hp_famBCC(1,2))
         g_hp_famBCC(2:5,2)=g_hp_famBCC(1,2)
         stateVars(342:346) = g_hp_famBCC(1:5,2)   ! hall-petch contributions for all others
            
      endif
      
      stateVars(321:326) = 0.d0     ! Green-Lagrange Creep Strain
      stateVars(327:332) = 0.d0     ! dp, plastic rate of deformation in current configuration (push fwd of Dp=symm(Lp))
      stateVars(333:336) = 0.d0     ! Free
      
      stateVars(347:351) = 0.d0  ! Free
      
      stateVars(352)       = 0.d0   ! normal stress on basal   HCP, output only
      stateVars(353)       = 0.d0   ! shear stress on basal    HCP, output only
      stateVars(354:356)= 0.d0   ! c-axis dir               HCP, output only
      stateVars(357:365)= 0.d0   ! Nye tensor               generic, input only
      stateVars(366)       = grainSize
      stateVars(367)       = 0.d0   ! NOT USED
      stateVars(368)       = 0.d0   ! norm of delta_gamma     generic, output only
      stateVars(369:377)= 0.d0   ! R_cryst, cur.lattice ori. (Row Major) -generic, output only
      
      END SUBROUTINE
      
      SUBROUTINE calc_HallPetchStrength(K_hp,D_hp,grainSize,g_hp)
      implicit none
      real(8), intent(in) :: K_hp, D_hp, grainSize
      real(8), intent(out):: g_hp
      if (D_hp < 0.d0) then
         g_hp = K_hp / dsqrt(grainSize)
      else
         g_hp = K_hp / dsqrt(D_hp)
      endif
      END SUBROUTINE
      
      SUBROUTINE calc_SlipSysStrengthsBCC(g0_BCC,g_BCC,g_hp_famBCC,delta_g_BCC_l, &
                                          temperature,thermalSoftening)
      implicit none
      real(8), intent(in) :: g_hp_famBCC(nslipfam,2)
      real(8), intent(in) :: delta_g_BCC_l(nslip_max)
      real(8), optional, intent(in) :: temperature
      logical, optional, intent(in) :: thermalSoftening
      real(8), intent(out):: g0_BCC(nslip_max,2)     ! for each slip sys, soft/hard directions
      real(8), intent(out):: g_BCC(nslip_max,2)
      ! locals
      real(8) :: g0_hp_BCC(nslip_max,2)      
      real(8) :: g0_softening_BCC(nslip_max)
      integer :: isys
      
      ! calculate thermal softening at each slip system
      g0_softening_BCC(:) = 0.d0
      if(present(temperature).and.present(thermalSoftening)) then
         if(thermalSoftening) then
            do isys = 1,nslip_BCC
               g0_softening_BCC(isys) = (1.d0 - dexp(-(temperature-300.d0)/softening_T_char_BCC(isys)))*softening_g_therm_BCC(isys)
            enddo
         endif
      endif
      
      ! distribute the Hall Petch contributions to slip systems
      call slipFamily_to_SystemsBCC(g_hp_famBCC(:,2),g0_hp_BCC,nslips_famBCC)    ! distribute to all slip systems
      call slipFamily_to_SystemsBCC_b1(g_hp_famBCC(:,1),g0_hp_BCC,nslips_famBCC) ! distribute only to b1=[-111] slip systems
      
      ! now add hall petch contributions to internal resistances of each slip system
      g0_BCC = g0_bar_BCC + g0_hp_BCC
      ! add work-hardening to get the current slip system resistance
      g_BCC(:,1) = g0_BCC(:,1) + delta_g_BCC_l(:) ! accumulated work-hardening is 
      g_BCC(:,2) = g0_BCC(:,2) + delta_g_BCC_l(:) ! independent of (soft/hard) slip direction
      
      END SUBROUTINE

      SUBROUTINE calc_SlipSysStrengthsHCP(g0_HCP,g_HCP,phase,g_hp_famHCP,delta_g_HCP_l, &
                                          temperature,thermalSoftening)

      implicit none
      real(8), intent(in) :: g_hp_famHCP(nslipfam,2) ! for each slip fam, tension/compression
      real(8), intent(in) :: delta_g_HCP_l(nslip_max)
      integer, intent(in) :: phase
      real(8), optional, intent(in) :: temperature
      logical, optional, intent(in) :: thermalSoftening
      real(8), intent(out):: g0_HCP(nslip_max,2)     ! for each slip sys, tension/compression
      real(8), intent(out):: g_HCP(nslip_max,2)
      ! locals
      real(8) :: g0_hp_HCP(nslip_max,2)      ! tension/compression   
      real(8) :: g0_softening_HCP(nslip_max)
      integer :: isys
      
      ! calculate thermal softening at each slip system
      g0_softening_HCP(:) = 0.d0
      if(present(temperature).and.present(thermalSoftening)) then
         if(thermalSoftening) then
            do isys = 1,12
               g0_softening_HCP(isys) = (1.d0 - dexp(-(temperature-300.d0)/softening_T_char_HCP(isys))) &
                                        *softening_g_therm_HCP(isys)
            enddo
            do isys = 13,30   ! <c+a> systems
               if(temperature < 400.d0) then
                  g0_softening_HCP(isys) = g_drop_c_a_300_400*(temperature-300.d0)/100.d0
               else
                  g0_softening_HCP(isys) = g_drop_c_a_300_400 + &
                              (1.d0 - dexp(-(temperature-400.d0)/softening_T_char_HCP(isys))) &
                              *softening_g_therm_HCP(isys)
               endif
            enddo
         endif
      endif
      
      ! distribute the Hall Petch contributions to slip systems
      call slipFamily_to_SystemsHCP(g_hp_famHCP(:,1),g0_hp_HCP(:,1),nslips_famHCP) ! tension values
      call slipFamily_to_SystemsHCP(g_hp_famHCP(:,2),g0_hp_HCP(:,2),nslips_famHCP) ! compression values
      
      ! now add hall petch contributions to internal resistances of each slip system
      g0_HCP(:,:) = g0_bar_HCP(:,:,phase) + g0_hp_HCP(:,:)
      ! add work-hardening to get the current slip system resistance
      g_HCP(:,1) = g0_HCP(:,1) + delta_g_HCP_l(:) - g0_softening_HCP(:) ! accumulated work-hardening and thermal softening are  
      g_HCP(:,2) = g0_HCP(:,2) + delta_g_HCP_l(:) - g0_softening_HCP(:) ! independent of tension/compression
      
      END SUBROUTINE
      
      
      SUBROUTINE calc_ElasticConstantsTemp_HCP(c11_T_HCP,c12_T_HCP,c13_T_HCP, &
                                               c33_T_HCP,c55_T_HCP, &
                                               temperature,thermalSoftening)

      implicit none
      real(8), intent(out) :: c11_T_HCP,c12_T_HCP,c13_T_HCP, &
                              c33_T_HCP,c55_T_HCP                 ! HCP elasticity coeffs at ambient temperature
      real(8), optional, intent(in) :: temperature
      logical, optional, intent(in) :: thermalSoftening
      
      if ( thermalSoftening ) then
         c11_T_HCP = c11_HCP - (temperature-300.d0)*c11_k_T_HCP
         c12_T_HCP = c12_HCP - (temperature-300.d0)*c12_k_T_HCP
         c13_T_HCP = c13_HCP - (temperature-300.d0)*c13_k_T_HCP
         c33_T_HCP = c33_HCP - (temperature-300.d0)*c33_k_T_HCP
         c55_T_HCP = c55_HCP - (temperature-300.d0)*c55_k_T_HCP
      else
         c11_T_HCP = c11_HCP
         c12_T_HCP = c12_HCP
         c13_T_HCP = c13_HCP
         c33_T_HCP = c33_HCP
         c55_T_HCP = c55_HCP
      endif
      END SUBROUTINE
      
      SUBROUTINE calc_ElasticConstantsTemp_BCC(c1_T_BCC,c2_T_BCC,c3_T_BCC, &
                                               temperature,thermalSoftening)

      implicit none
      real(8), intent(out) :: c1_T_BCC,c2_T_BCC,c3_T_BCC    ! BCC elasticity coeffs at ambient temperature
      real(8), optional, intent(in) :: temperature
      logical, optional, intent(in) :: thermalSoftening
      
      if ( thermalSoftening ) then
         c1_T_BCC = c1_BCC - (temperature-300.d0)*c1_k_T_BCC
         c2_T_BCC = c2_BCC - (temperature-300.d0)*c2_k_T_BCC
         c3_T_BCC = c3_BCC - (temperature-300.d0)*c3_k_T_BCC
      else
         c1_T_BCC = c1_BCC
         c2_T_BCC = c2_BCC
         c3_T_BCC = c3_BCC
      endif
      END SUBROUTINE
      
      
      
      SUBROUTINE calc_AdibaticHeating(w_p_t, dT_wp_t)
      implicit none
      real(8), intent(in) :: w_p_t
      real(8), intent(out):: dT_wp_t
      
      dT_wp_t = 0.9d0 * w_p_t / heatCapacityVolumetric
      
      END SUBROUTINE

      ! taken from from Kourosh's code
      ! usage: rho_GND = [B].Nye_vector,
      !   where Nye_vector is row-major ordered components 
      !   of the Nye's tensor in crystal coordinate frame
      SUBROUTINE B_MATRIX_HCP(matB)
      implicit none
      real(8), intent(out) :: matB(63,9)

      matB (  1 ,  1 ) =    1.128880E-01
      matB (  1 ,  2 ) =   -8.166213E-02
      matB (  1 ,  3 ) =    0.000000E+00
      matB (  1 ,  4 ) =   -8.166213E-02
      matB (  1 ,  5 ) =    2.071833E-01
      matB (  1 ,  6 ) =    0.000000E+00
      matB (  1 ,  7 ) =    0.000000E+00
      matB (  1 ,  8 ) =    0.000000E+00
      matB (  1 ,  9 ) =    5.776589E-02
      matB (  2 ,  1 ) =    1.128880E-01
      matB (  2 ,  2 ) =    8.166213E-02
      matB (  2 ,  3 ) =    0.000000E+00
      matB (  2 ,  4 ) =    8.166213E-02
      matB (  2 ,  5 ) =    2.071833E-01
      matB (  2 ,  6 ) =    0.000000E+00
      matB (  2 ,  7 ) =    0.000000E+00
      matB (  2 ,  8 ) =    0.000000E+00
      matB (  2 ,  9 ) =    5.776589E-02
      matB (  3 ,  1 ) =    2.543310E-01
      matB (  3 ,  2 ) =    0.000000E+00
      matB (  3 ,  3 ) =    0.000000E+00
      matB (  3 ,  4 ) =    0.000000E+00
      matB (  3 ,  5 ) =    6.574036E-02
      matB (  3 ,  6 ) =    0.000000E+00
      matB (  3 ,  7 ) =    0.000000E+00
      matB (  3 ,  8 ) =    0.000000E+00
      matB (  3 ,  9 ) =    5.776589E-02
      matB (  4 ,  1 ) =    7.333343E-02
      matB (  4 ,  2 ) =    2.306295E-02
      matB (  4 ,  3 ) =    4.359713E-02
      matB (  4 ,  4 ) =    2.306295E-02
      matB (  4 ,  5 ) =    9.996424E-02
      matB (  4 ,  6 ) =    7.551245E-02
      matB (  4 ,  7 ) =    2.376822E-02
      matB (  4 ,  8 ) =    4.116776E-02
      matB (  4 ,  9 ) =    1.377837E-01
      matB (  5 ,  1 ) =    7.333343E-02
      matB (  5 ,  2 ) =   -2.306295E-02
      matB (  5 ,  3 ) =   -4.359713E-02
      matB (  5 ,  4 ) =   -2.306295E-02
      matB (  5 ,  5 ) =    9.996424E-02
      matB (  5 ,  6 ) =    7.551245E-02
      matB (  5 ,  7 ) =   -2.376822E-02
      matB (  5 ,  8 ) =    4.116776E-02
      matB (  5 ,  9 ) =    1.377837E-01
      matB (  6 ,  1 ) =    1.132796E-01
      matB (  6 ,  2 ) =    0.000000E+00
      matB (  6 ,  3 ) =   -8.719426E-02
      matB (  6 ,  4 ) =    0.000000E+00
      matB (  6 ,  5 ) =    6.001803E-02
      matB (  6 ,  6 ) =    0.000000E+00
      matB (  6 ,  7 ) =   -4.753643E-02
      matB (  6 ,  8 ) =    0.000000E+00
      matB (  6 ,  9 ) =    1.377837E-01
      matB (  7 ,  1 ) =    7.333343E-02
      matB (  7 ,  2 ) =    2.306295E-02
      matB (  7 ,  3 ) =   -4.359713E-02
      matB (  7 ,  4 ) =    2.306295E-02
      matB (  7 ,  5 ) =    9.996424E-02
      matB (  7 ,  6 ) =   -7.551245E-02
      matB (  7 ,  7 ) =   -2.376822E-02
      matB (  7 ,  8 ) =   -4.116776E-02
      matB (  7 ,  9 ) =    1.377837E-01
      matB (  8 ,  1 ) =    7.333343E-02
      matB (  8 ,  2 ) =   -2.306295E-02
      matB (  8 ,  3 ) =    4.359713E-02
      matB (  8 ,  4 ) =   -2.306295E-02
      matB (  8 ,  5 ) =    9.996424E-02
      matB (  8 ,  6 ) =   -7.551245E-02
      matB (  8 ,  7 ) =    2.376822E-02
      matB (  8 ,  8 ) =   -4.116776E-02
      matB (  8 ,  9 ) =    1.377837E-01
      matB (  9 ,  1 ) =    1.132796E-01
      matB (  9 ,  2 ) =    0.000000E+00
      matB (  9 ,  3 ) =    8.719426E-02
      matB (  9 ,  4 ) =    0.000000E+00
      matB (  9 ,  5 ) =    6.001803E-02
      matB (  9 ,  6 ) =    0.000000E+00
      matB (  9 ,  7 ) =    4.753643E-02
      matB (  9 ,  8 ) =    0.000000E+00
      matB (  9 ,  9 ) =    1.377837E-01
      matB ( 10 ,  1 ) =    0.000000E+00
      matB ( 10 ,  2 ) =   -1.452474E-01
      matB ( 10 ,  3 ) =    8.376144E-02
      matB ( 10 ,  4 ) =   -2.046507E-02
      matB ( 10 ,  5 ) =    0.000000E+00
      matB ( 10 ,  6 ) =    0.000000E+00
      matB ( 10 ,  7 ) =    8.704741E-03
      matB ( 10 ,  8 ) =    0.000000E+00
      matB ( 10 ,  9 ) =    0.000000E+00
      matB ( 11 ,  1 ) =    7.175562E-02
      matB ( 11 ,  2 ) =   -2.096306E-02
      matB ( 11 ,  3 ) =    4.188072E-02
      matB ( 11 ,  4 ) =    1.038193E-01
      matB ( 11 ,  5 ) =   -7.175562E-02
      matB ( 11 ,  6 ) =    7.253953E-02
      matB ( 11 ,  7 ) =    4.352370E-03
      matB ( 11 ,  8 ) =    7.538527E-03
      matB ( 11 ,  9 ) =    0.000000E+00
      matB ( 12 ,  1 ) =   -7.175562E-02
      matB ( 12 ,  2 ) =   -2.096306E-02
      matB ( 12 ,  3 ) =   -4.188072E-02
      matB ( 12 ,  4 ) =    1.038193E-01
      matB ( 12 ,  5 ) =    7.175562E-02
      matB ( 12 ,  6 ) =    7.253953E-02
      matB ( 12 ,  7 ) =   -4.352370E-03
      matB ( 12 ,  8 ) =    7.538527E-03
      matB ( 12 ,  9 ) =    0.000000E+00
      matB ( 13 ,  1 ) =    0.000000E+00
      matB ( 13 ,  2 ) =   -1.452474E-01
      matB ( 13 ,  3 ) =   -8.376144E-02
      matB ( 13 ,  4 ) =   -2.046507E-02
      matB ( 13 ,  5 ) =    0.000000E+00
      matB ( 13 ,  6 ) =    0.000000E+00
      matB ( 13 ,  7 ) =   -8.704741E-03
      matB ( 13 ,  8 ) =    0.000000E+00
      matB ( 13 ,  9 ) =    0.000000E+00
      matB ( 14 ,  1 ) =    7.175562E-02
      matB ( 14 ,  2 ) =   -2.096306E-02
      matB ( 14 ,  3 ) =   -4.188072E-02
      matB ( 14 ,  4 ) =    1.038193E-01
      matB ( 14 ,  5 ) =   -7.175562E-02
      matB ( 14 ,  6 ) =   -7.253953E-02
      matB ( 14 ,  7 ) =   -4.352370E-03
      matB ( 14 ,  8 ) =   -7.538527E-03
      matB ( 14 ,  9 ) =    0.000000E+00
      matB ( 15 ,  1 ) =   -7.175562E-02
      matB ( 15 ,  2 ) =   -2.096306E-02
      matB ( 15 ,  3 ) =    4.188072E-02
      matB ( 15 ,  4 ) =    1.038193E-01
      matB ( 15 ,  5 ) =    7.175562E-02
      matB ( 15 ,  6 ) =   -7.253953E-02
      matB ( 15 ,  7 ) =    4.352370E-03
      matB ( 15 ,  8 ) =   -7.538527E-03
      matB ( 15 ,  9 ) =    0.000000E+00
      matB ( 16 ,  1 ) =   -3.224887E-03
      matB ( 16 ,  2 ) =   -3.859451E-02
      matB ( 16 ,  3 ) =    2.225672E-02
      matB ( 16 ,  4 ) =   -5.437889E-03
      matB ( 16 ,  5 ) =   -7.949124E-02
      matB ( 16 ,  6 ) =    2.497762E-02
      matB ( 16 ,  7 ) =    2.312986E-03
      matB ( 16 ,  8 ) =   -6.102021E-02
      matB ( 16 ,  9 ) =    4.509505E-02
      matB ( 17 ,  1 ) =   -4.135806E-02
      matB ( 17 ,  2 ) =    2.745409E-02
      matB ( 17 ,  3 ) =   -1.050289E-02
      matB ( 17 ,  4 ) =    6.061071E-02
      matB ( 17 ,  5 ) =   -4.135806E-02
      matB ( 17 ,  6 ) =    3.176370E-02
      matB ( 17 ,  7 ) =    5.400154E-02
      matB ( 17 ,  8 ) =   -2.850700E-02
      matB ( 17 ,  9 ) =    4.509505E-02
      matB ( 18 ,  1 ) =   -7.949124E-02
      matB ( 18 ,  2 ) =   -3.859451E-02
      matB ( 18 ,  3 ) =   -3.275961E-02
      matB ( 18 ,  4 ) =   -5.437889E-03
      matB ( 18 ,  5 ) =   -3.224887E-03
      matB ( 18 ,  6 ) =    6.786080E-03
      matB ( 18 ,  7 ) =    5.168856E-02
      matB ( 18 ,  8 ) =    3.251321E-02
      matB ( 18 ,  9 ) =    4.509505E-02
      matB ( 19 ,  1 ) =   -3.224887E-03
      matB ( 19 ,  2 ) =   -3.859451E-02
      matB ( 19 ,  3 ) =   -2.225672E-02
      matB ( 19 ,  4 ) =   -5.437889E-03
      matB ( 19 ,  5 ) =   -7.949124E-02
      matB ( 19 ,  6 ) =   -2.497762E-02
      matB ( 19 ,  7 ) =   -2.312986E-03
      matB ( 19 ,  8 ) =    6.102021E-02
      matB ( 19 ,  9 ) =    4.509505E-02
      matB ( 20 ,  1 ) =   -4.135806E-02
      matB ( 20 ,  2 ) =    2.745409E-02
      matB ( 20 ,  3 ) =    1.050289E-02
      matB ( 20 ,  4 ) =    6.061071E-02
      matB ( 20 ,  5 ) =   -4.135806E-02
      matB ( 20 ,  6 ) =   -3.176370E-02
      matB ( 20 ,  7 ) =   -5.400154E-02
      matB ( 20 ,  8 ) =    2.850700E-02
      matB ( 20 ,  9 ) =    4.509505E-02
      matB ( 21 ,  1 ) =   -7.949124E-02
      matB ( 21 ,  2 ) =   -3.859451E-02
      matB ( 21 ,  3 ) =    3.275961E-02
      matB ( 21 ,  4 ) =   -5.437889E-03
      matB ( 21 ,  5 ) =   -3.224887E-03
      matB ( 21 ,  6 ) =   -6.786080E-03
      matB ( 21 ,  7 ) =   -5.168856E-02
      matB ( 21 ,  8 ) =   -3.251321E-02
      matB ( 21 ,  9 ) =    4.509505E-02
      matB ( 22 ,  1 ) =   -3.224887E-03
      matB ( 22 ,  2 ) =    3.859451E-02
      matB ( 22 ,  3 ) =   -2.225672E-02
      matB ( 22 ,  4 ) =    5.437889E-03
      matB ( 22 ,  5 ) =   -7.949124E-02
      matB ( 22 ,  6 ) =    2.497762E-02
      matB ( 22 ,  7 ) =   -2.312986E-03
      matB ( 22 ,  8 ) =   -6.102021E-02
      matB ( 22 ,  9 ) =    4.509505E-02
      matB ( 23 ,  1 ) =   -7.949124E-02
      matB ( 23 ,  2 ) =    3.859451E-02
      matB ( 23 ,  3 ) =   -3.275961E-02
      matB ( 23 ,  4 ) =    5.437889E-03
      matB ( 23 ,  5 ) =   -3.224887E-03
      matB ( 23 ,  6 ) =   -6.786080E-03
      matB ( 23 ,  7 ) =    5.168856E-02
      matB ( 23 ,  8 ) =   -3.251321E-02
      matB ( 23 ,  9 ) =    4.509505E-02
      matB ( 24 ,  1 ) =   -4.135806E-02
      matB ( 24 ,  2 ) =   -2.745409E-02
      matB ( 24 ,  3 ) =   -1.050289E-02
      matB ( 24 ,  4 ) =   -6.061071E-02
      matB ( 24 ,  5 ) =   -4.135806E-02
      matB ( 24 ,  6 ) =   -3.176370E-02
      matB ( 24 ,  7 ) =    5.400154E-02
      matB ( 24 ,  8 ) =    2.850700E-02
      matB ( 24 ,  9 ) =    4.509505E-02
      matB ( 25 ,  1 ) =   -3.224887E-03
      matB ( 25 ,  2 ) =    3.859451E-02
      matB ( 25 ,  3 ) =    2.225672E-02
      matB ( 25 ,  4 ) =    5.437889E-03
      matB ( 25 ,  5 ) =   -7.949124E-02
      matB ( 25 ,  6 ) =   -2.497762E-02
      matB ( 25 ,  7 ) =    2.312986E-03
      matB ( 25 ,  8 ) =    6.102021E-02
      matB ( 25 ,  9 ) =    4.509505E-02
      matB ( 26 ,  1 ) =   -7.949124E-02
      matB ( 26 ,  2 ) =    3.859451E-02
      matB ( 26 ,  3 ) =    3.275961E-02
      matB ( 26 ,  4 ) =    5.437889E-03
      matB ( 26 ,  5 ) =   -3.224887E-03
      matB ( 26 ,  6 ) =    6.786080E-03
      matB ( 26 ,  7 ) =   -5.168856E-02
      matB ( 26 ,  8 ) =    3.251321E-02
      matB ( 26 ,  9 ) =    4.509505E-02
      matB ( 27 ,  1 ) =   -4.135806E-02
      matB ( 27 ,  2 ) =   -2.745409E-02
      matB ( 27 ,  3 ) =    1.050289E-02
      matB ( 27 ,  4 ) =   -6.061071E-02
      matB ( 27 ,  5 ) =   -4.135806E-02
      matB ( 27 ,  6 ) =    3.176370E-02
      matB ( 27 ,  7 ) =   -5.400154E-02
      matB ( 27 ,  8 ) =   -2.850700E-02
      matB ( 27 ,  9 ) =    4.509505E-02
      matB ( 28 ,  1 ) =   -2.481466E-02
      matB ( 28 ,  2 ) =    3.676235E-02
      matB ( 28 ,  3 ) =   -1.823385E-02
      matB ( 28 ,  4 ) =    3.676235E-02
      matB ( 28 ,  5 ) =   -6.726416E-02
      matB ( 28 ,  6 ) =    3.158195E-02
      matB ( 28 ,  7 ) =    2.876961E-02
      matB ( 28 ,  8 ) =   -4.983043E-02
      matB ( 28 ,  9 ) =    5.019939E-02
      matB ( 29 ,  1 ) =   -8.848891E-02
      matB ( 29 ,  2 ) =    0.000000E+00
      matB ( 29 ,  3 ) =   -3.646770E-02
      matB ( 29 ,  4 ) =    0.000000E+00
      matB ( 29 ,  5 ) =   -3.589914E-03
      matB ( 29 ,  6 ) =    0.000000E+00
      matB ( 29 ,  7 ) =    5.753922E-02
      matB ( 29 ,  8 ) =    0.000000E+00
      matB ( 29 ,  9 ) =    5.019939E-02
      matB ( 30 ,  1 ) =   -2.481466E-02
      matB ( 30 ,  2 ) =   -3.676235E-02
      matB ( 30 ,  3 ) =   -1.823385E-02
      matB ( 30 ,  4 ) =   -3.676235E-02
      matB ( 30 ,  5 ) =   -6.726416E-02
      matB ( 30 ,  6 ) =   -3.158195E-02
      matB ( 30 ,  7 ) =    2.876961E-02
      matB ( 30 ,  8 ) =    4.983043E-02
      matB ( 30 ,  9 ) =    5.019939E-02
      matB ( 31 ,  1 ) =   -2.481466E-02
      matB ( 31 ,  2 ) =    3.676235E-02
      matB ( 31 ,  3 ) =    1.823385E-02
      matB ( 31 ,  4 ) =    3.676235E-02
      matB ( 31 ,  5 ) =   -6.726416E-02
      matB ( 31 ,  6 ) =   -3.158195E-02
      matB ( 31 ,  7 ) =   -2.876961E-02
      matB ( 31 ,  8 ) =    4.983043E-02
      matB ( 31 ,  9 ) =    5.019939E-02
      matB ( 32 ,  1 ) =   -8.848891E-02
      matB ( 32 ,  2 ) =    0.000000E+00
      matB ( 32 ,  3 ) =    3.646770E-02
      matB ( 32 ,  4 ) =    0.000000E+00
      matB ( 32 ,  5 ) =   -3.589914E-03
      matB ( 32 ,  6 ) =    0.000000E+00
      matB ( 32 ,  7 ) =   -5.753922E-02
      matB ( 32 ,  8 ) =    0.000000E+00
      matB ( 32 ,  9 ) =    5.019939E-02
      matB ( 33 ,  1 ) =   -2.481466E-02
      matB ( 33 ,  2 ) =   -3.676235E-02
      matB ( 33 ,  3 ) =    1.823385E-02
      matB ( 33 ,  4 ) =   -3.676235E-02
      matB ( 33 ,  5 ) =   -6.726416E-02
      matB ( 33 ,  6 ) =    3.158195E-02
      matB ( 33 ,  7 ) =   -2.876961E-02
      matB ( 33 ,  8 ) =   -4.983043E-02
      matB ( 33 ,  9 ) =    5.019939E-02
      matB ( 34 ,  1 ) =    8.166213E-02
      matB ( 34 ,  2 ) =    2.385720E-02
      matB ( 34 ,  3 ) =    0.000000E+00
      matB ( 34 ,  4 ) =   -1.181525E-01
      matB ( 34 ,  5 ) =   -8.166213E-02
      matB ( 34 ,  6 ) =    0.000000E+00
      matB ( 34 ,  7 ) =    0.000000E+00
      matB ( 34 ,  8 ) =    0.000000E+00
      matB ( 34 ,  9 ) =    0.000000E+00
      matB ( 35 ,  1 ) =   -8.166213E-02
      matB ( 35 ,  2 ) =    2.385720E-02
      matB ( 35 ,  3 ) =    0.000000E+00
      matB ( 35 ,  4 ) =   -1.181525E-01
      matB ( 35 ,  5 ) =    8.166213E-02
      matB ( 35 ,  6 ) =    0.000000E+00
      matB ( 35 ,  7 ) =    0.000000E+00
      matB ( 35 ,  8 ) =    0.000000E+00
      matB ( 35 ,  9 ) =    0.000000E+00
      matB ( 36 ,  1 ) =    0.000000E+00
      matB ( 36 ,  2 ) =    1.653002E-01
      matB ( 36 ,  3 ) =    0.000000E+00
      matB ( 36 ,  4 ) =    2.329046E-02
      matB ( 36 ,  5 ) =    0.000000E+00
      matB ( 36 ,  6 ) =    0.000000E+00
      matB ( 36 ,  7 ) =    0.000000E+00
      matB ( 36 ,  8 ) =    0.000000E+00
      matB ( 36 ,  9 ) =    0.000000E+00
      matB ( 37 ,  1 ) =    0.000000E+00
      matB ( 37 ,  2 ) =    0.000000E+00
      matB ( 37 ,  3 ) =   -1.754553E-01
      matB ( 37 ,  4 ) =    0.000000E+00
      matB ( 37 ,  5 ) =    0.000000E+00
      matB ( 37 ,  6 ) =    0.000000E+00
      matB ( 37 ,  7 ) =   -1.823385E-02
      matB ( 37 ,  8 ) =    0.000000E+00
      matB ( 37 ,  9 ) =    0.000000E+00
      matB ( 38 ,  1 ) =    0.000000E+00
      matB ( 38 ,  2 ) =    0.000000E+00
      matB ( 38 ,  3 ) =   -8.772767E-02
      matB ( 38 ,  4 ) =    0.000000E+00
      matB ( 38 ,  5 ) =    0.000000E+00
      matB ( 38 ,  6 ) =   -1.519488E-01
      matB ( 38 ,  7 ) =   -9.116924E-03
      matB ( 38 ,  8 ) =   -1.579098E-02
      matB ( 38 ,  9 ) =    0.000000E+00
      matB ( 39 ,  1 ) =    0.000000E+00
      matB ( 39 ,  2 ) =    0.000000E+00
      matB ( 39 ,  3 ) =    8.772767E-02
      matB ( 39 ,  4 ) =    0.000000E+00
      matB ( 39 ,  5 ) =    0.000000E+00
      matB ( 39 ,  6 ) =   -1.519488E-01
      matB ( 39 ,  7 ) =    9.116924E-03
      matB ( 39 ,  8 ) =   -1.579098E-02
      matB ( 39 ,  9 ) =    0.000000E+00
      matB ( 40 ,  1 ) =    0.000000E+00
      matB ( 40 ,  2 ) =    7.891340E-02
      matB ( 40 ,  3 ) =    1.541707E-01
      matB ( 40 ,  4 ) =    1.111874E-02
      matB ( 40 ,  5 ) =    0.000000E+00
      matB ( 40 ,  6 ) =    0.000000E+00
      matB ( 40 ,  7 ) =    1.602188E-02
      matB ( 40 ,  8 ) =    0.000000E+00
      matB ( 40 ,  9 ) =    0.000000E+00
      matB ( 41 ,  1 ) =   -3.898506E-02
      matB ( 41 ,  2 ) =    1.138930E-02
      matB ( 41 ,  3 ) =    7.708534E-02
      matB ( 41 ,  4 ) =   -5.640537E-02
      matB ( 41 ,  5 ) =    3.898506E-02
      matB ( 41 ,  6 ) =    1.335157E-01
      matB ( 41 ,  7 ) =    8.010941E-03
      matB ( 41 ,  8 ) =    1.387536E-02
      matB ( 41 ,  9 ) =    0.000000E+00
      matB ( 42 ,  1 ) =    3.898506E-02
      matB ( 42 ,  2 ) =    1.138930E-02
      matB ( 42 ,  3 ) =   -7.708534E-02
      matB ( 42 ,  4 ) =   -5.640537E-02
      matB ( 42 ,  5 ) =   -3.898506E-02
      matB ( 42 ,  6 ) =    1.335157E-01
      matB ( 42 ,  7 ) =   -8.010941E-03
      matB ( 42 ,  8 ) =    1.387536E-02
      matB ( 42 ,  9 ) =    0.000000E+00
      matB ( 43 ,  1 ) =    0.000000E+00
      matB ( 43 ,  2 ) =    7.891340E-02
      matB ( 43 ,  3 ) =   -1.541707E-01
      matB ( 43 ,  4 ) =    1.111874E-02
      matB ( 43 ,  5 ) =    0.000000E+00
      matB ( 43 ,  6 ) =    0.000000E+00
      matB ( 43 ,  7 ) =   -1.602188E-02
      matB ( 43 ,  8 ) =    0.000000E+00
      matB ( 43 ,  9 ) =    0.000000E+00
      matB ( 44 ,  1 ) =   -3.898506E-02
      matB ( 44 ,  2 ) =    1.138930E-02
      matB ( 44 ,  3 ) =   -7.708534E-02
      matB ( 44 ,  4 ) =   -5.640537E-02
      matB ( 44 ,  5 ) =    3.898506E-02
      matB ( 44 ,  6 ) =   -1.335157E-01
      matB ( 44 ,  7 ) =   -8.010941E-03
      matB ( 44 ,  8 ) =   -1.387536E-02
      matB ( 44 ,  9 ) =    0.000000E+00
      matB ( 45 ,  1 ) =    3.898506E-02
      matB ( 45 ,  2 ) =    1.138930E-02
      matB ( 45 ,  3 ) =    7.708534E-02
      matB ( 45 ,  4 ) =   -5.640537E-02
      matB ( 45 ,  5 ) =   -3.898506E-02
      matB ( 45 ,  6 ) =   -1.335157E-01
      matB ( 45 ,  7 ) =    8.010941E-03
      matB ( 45 ,  8 ) =   -1.387536E-02
      matB ( 45 ,  9 ) =    0.000000E+00
      matB ( 46 ,  1 ) =   -4.988725E-02
      matB ( 46 ,  2 ) =   -4.762038E-03
      matB ( 46 ,  3 ) =   -4.005471E-03
      matB ( 46 ,  4 ) =   -7.255670E-02
      matB ( 46 ,  5 ) =    9.432920E-03
      matB ( 46 ,  6 ) =    2.081304E-02
      matB ( 46 ,  7 ) =   -7.021235E-02
      matB ( 46 ,  8 ) =    1.134682E-02
      matB ( 46 ,  9 ) =    2.205483E-02
      matB ( 47 ,  1 ) =    2.808287E-02
      matB ( 47 ,  2 ) =    2.754063E-02
      matB ( 47 ,  3 ) =   -2.002735E-02
      matB ( 47 ,  4 ) =   -4.025403E-02
      matB ( 47 ,  5 ) =   -6.853720E-02
      matB ( 47 ,  6 ) =    6.937678E-03
      matB ( 47 ,  7 ) =   -4.493281E-02
      matB ( 47 ,  8 ) =   -5.513227E-02
      matB ( 47 ,  9 ) =    2.205483E-02
      matB ( 48 ,  1 ) =   -3.887712E-02
      matB ( 48 ,  2 ) =    7.891340E-02
      matB ( 48 ,  3 ) =   -1.602188E-02
      matB ( 48 ,  4 ) =    1.111874E-02
      matB ( 48 ,  5 ) =   -1.577209E-03
      matB ( 48 ,  6 ) =   -1.387536E-02
      matB ( 48 ,  7 ) =    2.527954E-02
      matB ( 48 ,  8 ) =   -6.647909E-02
      matB ( 48 ,  9 ) =    2.205483E-02
      matB ( 49 ,  1 ) =   -4.988725E-02
      matB ( 49 ,  2 ) =   -4.762038E-03
      matB ( 49 ,  3 ) =    4.005471E-03
      matB ( 49 ,  4 ) =   -7.255670E-02
      matB ( 49 ,  5 ) =    9.432920E-03
      matB ( 49 ,  6 ) =   -2.081304E-02
      matB ( 49 ,  7 ) =    7.021235E-02
      matB ( 49 ,  8 ) =   -1.134682E-02
      matB ( 49 ,  9 ) =    2.205483E-02
      matB ( 50 ,  1 ) =    2.808287E-02
      matB ( 50 ,  2 ) =    2.754063E-02
      matB ( 50 ,  3 ) =    2.002735E-02
      matB ( 50 ,  4 ) =   -4.025403E-02
      matB ( 50 ,  5 ) =   -6.853720E-02
      matB ( 50 ,  6 ) =   -6.937678E-03
      matB ( 50 ,  7 ) =    4.493281E-02
      matB ( 50 ,  8 ) =    5.513227E-02
      matB ( 50 ,  9 ) =    2.205483E-02
      matB ( 51 ,  1 ) =   -3.887712E-02
      matB ( 51 ,  2 ) =    7.891340E-02
      matB ( 51 ,  3 ) =    1.602188E-02
      matB ( 51 ,  4 ) =    1.111874E-02
      matB ( 51 ,  5 ) =   -1.577209E-03
      matB ( 51 ,  6 ) =    1.387536E-02
      matB ( 51 ,  7 ) =   -2.527954E-02
      matB ( 51 ,  8 ) =    6.647909E-02
      matB ( 51 ,  9 ) =    2.205483E-02
      matB ( 52 ,  1 ) =    4.988725E-02
      matB ( 52 ,  2 ) =   -4.762038E-03
      matB ( 52 ,  3 ) =   -4.005471E-03
      matB ( 52 ,  4 ) =   -7.255670E-02
      matB ( 52 ,  5 ) =   -9.432920E-03
      matB ( 52 ,  6 ) =   -2.081304E-02
      matB ( 52 ,  7 ) =   -7.021235E-02
      matB ( 52 ,  8 ) =   -1.134682E-02
      matB ( 52 ,  9 ) =   -2.205483E-02
      matB ( 53 ,  1 ) =    3.887712E-02
      matB ( 53 ,  2 ) =    7.891340E-02
      matB ( 53 ,  3 ) =    1.602188E-02
      matB ( 53 ,  4 ) =    1.111874E-02
      matB ( 53 ,  5 ) =    1.577209E-03
      matB ( 53 ,  6 ) =   -1.387536E-02
      matB ( 53 ,  7 ) =   -2.527954E-02
      matB ( 53 ,  8 ) =   -6.647909E-02
      matB ( 53 ,  9 ) =   -2.205483E-02
      matB ( 54 ,  1 ) =   -2.808287E-02
      matB ( 54 ,  2 ) =    2.754063E-02
      matB ( 54 ,  3 ) =    2.002735E-02
      matB ( 54 ,  4 ) =   -4.025403E-02
      matB ( 54 ,  5 ) =    6.853720E-02
      matB ( 54 ,  6 ) =    6.937678E-03
      matB ( 54 ,  7 ) =    4.493281E-02
      matB ( 54 ,  8 ) =   -5.513227E-02
      matB ( 54 ,  9 ) =   -2.205483E-02
      matB ( 55 ,  1 ) =    4.988725E-02
      matB ( 55 ,  2 ) =   -4.762038E-03
      matB ( 55 ,  3 ) =    4.005471E-03
      matB ( 55 ,  4 ) =   -7.255670E-02
      matB ( 55 ,  5 ) =   -9.432920E-03
      matB ( 55 ,  6 ) =    2.081304E-02
      matB ( 55 ,  7 ) =    7.021235E-02
      matB ( 55 ,  8 ) =    1.134682E-02
      matB ( 55 ,  9 ) =   -2.205483E-02
      matB ( 56 ,  1 ) =    3.887712E-02
      matB ( 56 ,  2 ) =    7.891340E-02
      matB ( 56 ,  3 ) =   -1.602188E-02
      matB ( 56 ,  4 ) =    1.111874E-02
      matB ( 56 ,  5 ) =    1.577209E-03
      matB ( 56 ,  6 ) =    1.387536E-02
      matB ( 56 ,  7 ) =    2.527954E-02
      matB ( 56 ,  8 ) =    6.647909E-02
      matB ( 56 ,  9 ) =   -2.205483E-02
      matB ( 57 ,  1 ) =   -2.808287E-02
      matB ( 57 ,  2 ) =    2.754063E-02
      matB ( 57 ,  3 ) =   -2.002735E-02
      matB ( 57 ,  4 ) =   -4.025403E-02
      matB ( 57 ,  5 ) =    6.853720E-02
      matB ( 57 ,  6 ) =   -6.937678E-03
      matB ( 57 ,  7 ) =   -4.493281E-02
      matB ( 57 ,  8 ) =    5.513227E-02
      matB ( 57 ,  9 ) =   -2.205483E-02
      matB ( 58 ,  1 ) =    4.339781E-02
      matB ( 58 ,  2 ) =    1.267846E-02
      matB ( 58 ,  3 ) =   -1.337656E-02
      matB ( 58 ,  4 ) =   -6.278993E-02
      matB ( 58 ,  5 ) =   -4.339781E-02
      matB ( 58 ,  6 ) =   -7.722959E-03
      matB ( 58 ,  7 ) =   -6.408926E-02
      matB ( 58 ,  8 ) =   -3.700195E-02
      matB ( 58 ,  9 ) =    0.000000E+00
      matB ( 59 ,  1 ) =    0.000000E+00
      matB ( 59 ,  2 ) =    8.784567E-02
      matB ( 59 ,  3 ) =    0.000000E+00
      matB ( 59 ,  4 ) =    1.237728E-02
      matB ( 59 ,  5 ) =    0.000000E+00
      matB ( 59 ,  6 ) =   -1.544592E-02
      matB ( 59 ,  7 ) =    0.000000E+00
      matB ( 59 ,  8 ) =   -7.400390E-02
      matB ( 59 ,  9 ) =    0.000000E+00
      matB ( 60 ,  1 ) =   -4.339781E-02
      matB ( 60 ,  2 ) =    1.267846E-02
      matB ( 60 ,  3 ) =    1.337656E-02
      matB ( 60 ,  4 ) =   -6.278993E-02
      matB ( 60 ,  5 ) =    4.339781E-02
      matB ( 60 ,  6 ) =   -7.722959E-03
      matB ( 60 ,  7 ) =    6.408926E-02
      matB ( 60 ,  8 ) =   -3.700195E-02
      matB ( 60 ,  9 ) =    0.000000E+00
      matB ( 61 ,  1 ) =    4.339781E-02
      matB ( 61 ,  2 ) =    1.267846E-02
      matB ( 61 ,  3 ) =    1.337656E-02
      matB ( 61 ,  4 ) =   -6.278993E-02
      matB ( 61 ,  5 ) =   -4.339781E-02
      matB ( 61 ,  6 ) =    7.722959E-03
      matB ( 61 ,  7 ) =    6.408926E-02
      matB ( 61 ,  8 ) =    3.700195E-02
      matB ( 61 ,  9 ) =    0.000000E+00
      matB ( 62 ,  1 ) =    0.000000E+00
      matB ( 62 ,  2 ) =    8.784567E-02
      matB ( 62 ,  3 ) =    0.000000E+00
      matB ( 62 ,  4 ) =    1.237728E-02
      matB ( 62 ,  5 ) =    0.000000E+00
      matB ( 62 ,  6 ) =    1.544592E-02
      matB ( 62 ,  7 ) =    0.000000E+00
      matB ( 62 ,  8 ) =    7.400390E-02
      matB ( 62 ,  9 ) =    0.000000E+00
      matB ( 63 ,  1 ) =   -4.339781E-02
      matB ( 63 ,  2 ) =    1.267846E-02
      matB ( 63 ,  3 ) =   -1.337656E-02
      matB ( 63 ,  4 ) =   -6.278993E-02
      matB ( 63 ,  5 ) =    4.339781E-02
      matB ( 63 ,  6 ) =    7.722959E-03
      matB ( 63 ,  7 ) =   -6.408926E-02
      matB ( 63 ,  8 ) =    3.700195E-02
      matB ( 63 ,  9 ) =    0.000000E+00
      
      END SUBROUTINE

      ! adapted from Kourosh's code. 
      ! Modified input as vectorGND_b. This is the lattice incompatibility densities for each GND family (rho x b, 1/length)
      ! these values are divided by the relevant burgers vectors of the slip system to give GNDs as outputs (1/area)
      SUBROUTINE getGNDComponents_HCP(ScrwGND,EnGND,EtGND,vectorGND_b)
      implicit none
      real(8), intent(out) :: ScrwGND(30)
      real(8), intent(out) :: EnGND(30)
      real(8), intent(out) :: EtGND(30)
      real(8), intent(in)  :: vectorGND_b(63)
      
      integer :: isys

      !*************************************************************
      !     GND DENSITY FOR SOME SLIP SYSTEMS ARE EQUAL, INCLUDING
      !     				SCREW GND: 
      !				 1=6=9=12
      !				 2=5=8=11
      !				 3=4=7=10
      !				 13=24=30
      !				 14=19=25
      !				 15=20=26
      !				 16=21=27
      !				 17=22=28
      !				 18=23=29

      !    DENSITIES FOR 6 EDGE_T GNDS ARE EQUAL TO 6 EDGE_N GNDS
      !               EDGE_T GND         EDGE_N GND
      !                   1        =         6
      !                   2        =         5
      !                   3        =         4
      !                   4        =         3
      !                   5        =         2
      !                   6        =         1

      !    GND DENSITIES ARE EQUALLY DISTRIBUTED BETWEEN THESE SYSTEMS

      ScrwGND(1)=vectorGND_b(1)/4.D0
      ScrwGND(2)=vectorGND_b(2)/4.D0
      ScrwGND(3)=vectorGND_b(3)/4.D0
      ScrwGND(4)=ScrwGND(3)	
      ScrwGND(5)=ScrwGND(2)	
      ScrwGND(6)=ScrwGND(1)
      ScrwGND(7)=ScrwGND(3)
      ScrwGND(8)=ScrwGND(2)
      ScrwGND(9)=ScrwGND(1)
      ScrwGND(10)=ScrwGND(3)
      ScrwGND(11)=ScrwGND(2)	
      ScrwGND(12)=ScrwGND(1)
      ScrwGND(13)=vectorGND_b(4)/3.D0
      ScrwGND(14)=vectorGND_b(5)/3.D0	
      ScrwGND(15)=vectorGND_b(6)/3.D0
      ScrwGND(16)=vectorGND_b(7)/3.D0
      ScrwGND(17)=vectorGND_b(8)/3.D0
      ScrwGND(18)=vectorGND_b(9)/3.D0
      ScrwGND(19)=ScrwGND(14)
      ScrwGND(20)=ScrwGND(15)
      ScrwGND(21)=ScrwGND(16)
      ScrwGND(22)=ScrwGND(17)
      ScrwGND(23)=ScrwGND(18)
      ScrwGND(24)=ScrwGND(13)
      ScrwGND(25)=ScrwGND(14)
      ScrwGND(26)=ScrwGND(15)
      ScrwGND(27)=ScrwGND(16)
      ScrwGND(28)=ScrwGND(17)
      ScrwGND(29)=ScrwGND(18)
      ScrwGND(30)=ScrwGND(13)

      EnGND(7:30)=vectorGND_b(10:33)
      EnGND(6)=vectorGND_b(34)/2.D0
      EnGND(5)=vectorGND_b(35)/2.D0
      EnGND(4)=vectorGND_b(36)/2.D0
      EnGND(3)=vectorGND_b(37)/2.D0
      EnGND(2)=vectorGND_b(38)/2.D0
      EnGND(1)=vectorGND_b(39)/2.D0


      EtGND(1)=EnGND(6)
      EtGND(2)=EnGND(5)
      EtGND(3)=EnGND(4)
      EtGND(4)=EnGND(3)
      EtGND(5)=EnGND(2)
      EtGND(6)=EnGND(1)
      EtGND(7:30)=vectorGND_b(40:63)
      
      ! divide by respective burgers vectors to get the true dislocation densities (number per area)
      ! assumed units: vectorGND_b ~ curlxFp: [1/um]. burgers: [nm]
      ! conversion multiplier for burgers vectors [um/nm]=1E-3, to get rho_GND [1/um^2]
      ! materials module parameter : nm_to_um = 1.d-3
      do isys=1,30
         ScrwGND(isys) = ScrwGND(isys) / (burgers(isys) * nm_to_um)
         EnGND(isys) = EnGND(isys) / (burgers(isys) * nm_to_um)
         EtGND(isys) = EtGND(isys) / (burgers(isys) * nm_to_um)
      enddo
      
      END SUBROUTINE
      end module
