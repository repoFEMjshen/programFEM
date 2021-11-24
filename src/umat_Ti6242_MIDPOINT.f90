! INPUT:
!  STATEV         - state variables at time t
!  NPROPS
!  F_t         - deformation gradient at time t
!  F_tau         - deformation gradient at time tau
!  DTIME          - dt = time tau - time t
!  TEMP(2)        - ambient temperature at time t and tau (Kelvins)
!  TEMP_0         - initial (stress free) ambient temperature (Kelvins)
!  PROPS(NPROPS)  - material properties
!  KINC           - simulation time step
!  JELEM          - element ID
!
! OUTPUT:
!  stress_mix(6)  - cauchy stress at time tau (homogenized if within trans-Beta grain) 
!  STATEV(NSTATV) - state variables at time tau
!  DDSDDE_mix     - tangent elasticity matrix dS/de (homogenized if within trans-Beta grain) 
!  umat_Err       - error flag (0: success. otherwise, error)
!                 0: success
!                 1: tau/g > 2, or flow rule produced NaN
!                 2: g is negative
!                 3: delGamma > 0.02 (algorithm assumes delGamma << 1)
!                 4: 1st level iteration did not convergence
!                 5: 2nd level iteration did not convergence
!                 6: back-stress iteration did not converge
!                 7: matrix inversion error
!  PNEWDT         - suggested time cut-back ratio (error: 0.5, success: 1)

      SUBROUTINE UMAT_CPFEM_THERM(stress_mix,STATEV,PROPS,DDSDDE_mix,NTENS,  &
                  NSTATV,NPROPS,PNEWDT,F_t,F_tau,TEMP,TEMP_0,KINC,  &
                  DTIME,JELEM,umat_Err)
                  
      use options, only: USE_POWER_LAW, wt_BCC
      use material_Ti6242

      implicit none

      ! added by deniz
      ! check if the umat and the imported model is compatible.
      ! this umat has 351 state variables. update this number as you develop the material model.
      ! parameter (nstatv_umat=352) ! temperature added as a state variable
      integer, parameter :: nstatv_umat=386

      ! arguments
      real(8), intent(out)  :: stress_mix(ntens)         ! Cauchy stress of the phase (homogenized if trans-B)
      real(8), intent(inout):: statev(nstatv)
      real(8), intent(in)   :: props(nprops)
      real(8), intent(out)  :: ddsdde_mix(ntens,ntens)   ! Tangent elasticity matrix of the phase (homogenized if trans-B)
      integer, intent(in)   :: ntens
      integer, intent(in)   :: nstatv
      integer, intent(in)   :: nprops
      real(8), intent(out)  :: pnewdt
      real(8), intent(in)   :: F_t(3,3)
      real(8), intent(in)   :: F_tau(3,3)
      real(8), intent(in)   :: TEMP(2) ! ambient temperature: TEMP(1) at time t. TEMP(2) at time tau.
      real(8), intent(in)   :: TEMP_0  ! initial/stress-free temperature
      real(8), intent(in)   :: dtime
      integer, intent(in)   :: kinc
      integer, intent(in)   :: jelem
      integer, intent(out)  :: umat_Err              ! deniz - error flag

      ! local variables
      real(8) :: g_hp_famBCC(nslipfam,2)  ! HP contrib. to slip resistance (slip family, b1/others) (read from statev)
      real(8) :: delta_g_BCC(nslip_max)   ! work-hardening part of slip resistance (slip systems)   (read from statev)
      real(8) :: g_hp_famHCP(nslipfam,2)  ! HP contrib. to slip resistance (slip family, tens/comp) (read from statev)
      real(8) :: delta_g_HCP(nslip_max)   ! work-hardening part of slip resistance (slip systems)   (read from statev)
      
      real(8) :: re_tau_t(3,3)
      real(8) :: euler(3)
      real(8) :: g_alpha_t(nslip_max),g_alpha_tau(nslip_max)
      real(8) :: Fp_t(3,3),Fp_tau(3,3)
      real(8) :: dsrot(3,nslip_max),dmrot(3,nslip_max)
      real(8) :: s_tan(3,3,nslip_max)
      real(8) :: ddsdde_4d(3,3,3,3)
      real(8) :: de_gg(3,3,3,3),de_gg_2d(6,6)
      real(8) :: gst(3,nslip_max),gmt(3,nslip_max)
      real(8) :: delta_gamma(nslip_max),tot_gamma(nslip_max),dw_p_mix,dw_p,w_p_t
      real(8) :: dGamma_dTau(nslip_max)
      real(8) :: tpk_1d_t(6),tpk_1d_tau(6)
      real(8) :: c_alpha_1d(6,nslip_max)
      real(8) :: tau(nslip_max)
      real(8) :: Fe_tau(3,3)
      real(8) :: r_el_cry(3,3),u_el_cry(3,3)
      real(8) :: r_tot_cry(3,3)
      real(8) :: stress(6)                  ! Cauchy stress                     of the crystal HCP/BCC
      real(8) :: strain_1d(6)               ! Total strain (Green-Lagrange, 1D) of the crystal HCP/BCC
      real(8) :: crp_strain_1d(6)           ! Creep strain (Green-Lagrange, 1D) of the crystal HCP/BCC
      real(8) :: ddsdde(6,6)                ! Tangent stiffness                 of the crystal HCP/BCC
      real(8) :: strain_mix(6)              ! Total strain (Logarithmic, 1D)    of the phase (homogenized if trans-B)
      real(8) :: creep_mix(6)               ! Creep strain (Green-Lagrange, 1D) of the phase (homogenized if trans-B)
      real(8) :: wt(2)

      real(8) :: ps(3)
      real(8) :: cauchy_2d(3,3)
      real(8) :: Fe_tau_inv(3,3)
      real(8) :: cbasal(3), b_edge(3)
      real(8) :: bas_tract(3)
      real(8) :: g1_bcc(5),g2_bcc(5)
      real(8) :: s_gg_2d(6,6)
      ! crystal/phase related:
      integer :: nslip      ! # of slip systems
      integer :: iCrystal
      integer :: iPhase

      ! dummy variables
      real(8) :: bc1,bc2,bc3,bc4,bc5   ! orientation relation
      real(8) :: R3,x1,x2
      real(8) :: pnewdt1
      integer :: i,j,isys,ic,ICOUNT,index,index1,index2
      integer :: ikin_init_state
      integer :: INFO,isign
      
      ! added by deniz
      integer :: sys_basal
      real(8) :: stress_normal_basal, stress_shear_basal,   &
                                 stress_shearvec_basal(3)
      
      ! ----- added by deniz ---------
      real(8) :: F_thrm_tau(3,3)       ! thermal deformation gradient at tau
      real(8) :: F_thrm_tau_inv(3,3)   ! ... inverse of it
      real(8) :: F_thermoelastic_tau(3,3)      ! thermo-elastic deformation grad at tau,
      real(8) :: F_thermoelastic_tau_inv(3,3)  ! ... inverse of it
      real(8) :: det_thrm              ! determinant of the thermal deformation (volume change)
      real(8) :: Q(3,3)                ! coordinate (passive) transformation matrix from crystal frame to spatial frame. {v}_i = Q_ij {v_cry}_j
      real(8) :: rot(3,3)              ! an active transformation matrix takes vectors from crystal frame to spatial frame. {v'}_i = R_ij {v}_j
      real(8) :: rot_HCP(3,3),rot_BCCtoHCP(3,3) ! rotation matrices between HCP/BCC phases
      real(8) :: Nye(3,3)              ! Nye tensor
                  
      Q(:,:) = 0.D0
      rot(:,:) = 0.D0
     
      if(nstatv_umat.NE.NSTATV)then
        write(*,*) 'UMAT and the imported problem is incompatible.'
        write(*,*) 'This UMAT expects ',nstatv_umat,' state variables'
        write(*,*) 'model/input file has ',NSTATV,' state variables'
        STOP
      endif
      ! -------- deniz --------

      pnewdt1=1.0
      umat_Err = 0



      iPhase = statev(304) ! PHASE NUMBER OF THE ELEMENT: 1=PRIMARY-alpha, 2=TRANSFORMED-beta
      
      if(iPhase.NE.1.AND.iPhase.NE.2)then
         write(*,*) 'element phase number: ',iPhase
         write(*,*) 'allowed values: 1- primary-Alpha, 2- transformed-Beta'
         STOP
      endif

      if(iPhase==pri_A)then ! element is pure HCP
         wt(1)=1.d0
         wt(2)=0.d0
      elseif(iPhase==trans_B)then ! element is HCP+BCC. homogenize with rule of mixtures. weights:
         wt(1)= 1.00d0 - wt_BCC   ! prev:0.93 / 0.07 -- paper: 88/12
         wt(2)= wt_BCC
      endif
      ! homogenized variables:
      creep_mix = 0.d0
      strain_mix = 0.d0
      stress_mix = 0.d0
      ddsdde_mix = 0.d0
      dw_p_mix    = 0.d0
      w_p_t = statev(299)
      
      do iCrystal=1,iPhase     
         ! if iPhase=2: transformed beta grain
         !     iCrystal=1: HCP within transformed beta grains
         !     iCrystal=2: BCC within transformed beta grains
         !
         ! if iPhase=1: primary alpha grain
         !     iCrystal=1: HCP within primary alpha grains
         
         if(iCrystal==HCP) nslip = nslip_HCP ! 30
         if(iCrystal==BCC) nslip = nslip_BCC ! 48
         
         if(iCrystal==HCP)then   ! HCP in primary-Alpha or transformed-Beta
         
            ! assign reference slip rates
            a_dot(1:nslip) = a_dot_HCP(1:nslip)
            a_dot_therm(1:nslip) = a_dot_therm_HCP(1:nslip)
            rate_exp(1:nslip) = rate_exp_HCP(1:nslip)
            E_act(1:nslip) = E_act_HCP(1:nslip)
            burgers(1:nslip) = burgers_HCP(1:nslip)
                        
            ! elasticity matrix -- HCP
            CALL crys_ElastMat_HCP(ddsdde,c11_HCP,c12_HCP,c13_HCP,c33_HCP,c55_HCP) 
            
            call tr2to4(ddsdde,ddsdde_4d)
            
            euler(1)=STATEV(301)
            euler(2)=STATEV(302)
            euler(3)=STATEV(303)
            call euler_slip(euler,rot)
            call rotate_crys_vector(cst_HCP,rot,gst,nslip_HCP)
            call rotate_crys_vector(cmt_HCP,rot,gmt,nslip_HCP)
            call rotate4d(ddsdde_4d,rot,de_gg)
            call tr4to2(de_gg,de_gg_2d)      
            ! deniz - comment
            ! here:
            ! rot is an active transformation matrix that rotates global axes to iCrystal axes. {v'}_i = R_ij {v}_j
            ! or, a (passive) coordinate transformation matrix from iCrystal frame to global frame. {v_cry}_i = R_ij {v}_j
            ! cst, cmt are the slip directions and plane normals in the natural iCrystal frame
            ! gst(3,slip), gmt(3,slip) are slip directions and plane normals in the undeformed config of the grain
            
            
!----------------------------------------------    
         endif
         
         if(iCrystal==BCC)then    ! BCC in trans-Beta
         
            ! assign reference slip rates, slip sens. and activation energies
            ! ASSUME THESE PARAMETERS ARE THE SAME ALONG SOFT/HARD DIRECTIONS
            ! AS IN THE PAPERS
            ! OTHERWISE, THESE ASSIGNMENTS HAVE TO BE CARRIED TO update_implicit
            ! and Stress state/flow directions must be taken into account
            a_dot(1:nslip) = a_dot_BCC(1:nslip,1)
            a_dot_therm(1:nslip) = a_dot_therm_BCC(1:nslip,1)
            rate_exp(1:nslip) = rate_exp_BCC(1:nslip,1)
            E_act(1:nslip) = E_act_BCC(1:nslip,1)
            burgers(1:nslip) = burgers_BCC(1:nslip,1)
         
            ! create elasticity matrix -- BCC
            CALL crys_ElastMat_BCC(ddsdde,c1_BCC,c2_BCC,c3_BCC) 
            
            call tr2to4(ddsdde,ddsdde_4d)

            euler(1)=STATEV(301)
            euler(2)=STATEV(302)
            euler(3)=STATEV(303)
            call euler_slip(euler,rot_HCP)
            
!-----------------Orientation Relation------------------------        
            ! rot_BCCtoHCP is an active rotation matrix 
            ! from BCC to HCP within a single colony
            !(HCP/BCC related via Burger's orientation)            
            bc1=dsqrt(2.d0)
            bc2=2.d0*dsqrt(3.d0)
            bc3=bc2**2.d0+6.d0
            bc4=dsqrt(bc3)
            bc5=(bc4-bc2)/12.d0             
            rot_BCCtoHCP(1,1)=bc5
            rot_BCCtoHCP(1,2)=2.d0*bc5+dsqrt(3.d0)/2.d0
            rot_BCCtoHCP(1,3)=-bc5
            rot_BCCtoHCP(3,1)=1.d0/bc1
            rot_BCCtoHCP(3,2)=0.0d0
            rot_BCCtoHCP(3,3)=1.d0/bc1
            rot_BCCtoHCP(2,1)=rot_BCCtoHCP(3,2)*rot_BCCtoHCP(1,3)-  &
              rot_BCCtoHCP(1,2)*rot_BCCtoHCP(3,3)
            rot_BCCtoHCP(2,2)=rot_BCCtoHCP(3,3)*rot_BCCtoHCP(1,1)-  &
              rot_BCCtoHCP(3,1)*rot_BCCtoHCP(1,3)
            rot_BCCtoHCP(2,3)=rot_BCCtoHCP(3,1)*rot_BCCtoHCP(1,2)-  &
              rot_BCCtoHCP(1,1)*rot_BCCtoHCP(3,2)       
              
!              rot_BCCtoHCP = 0.0
              
!              rot_BCCtoHCP(1,1) = - dsin(0.103d0) / dsqrt(2.d0)
!              rot_BCCtoHCP(2,1) = - dcos(0.103d0) / dsqrt(2.d0)
!              rot_BCCtoHCP(3,1) = + 1.d0 / dsqrt(2.d0)

!              rot_BCCtoHCP(1,2) = + dcos(0.103d0)
!              rot_BCCtoHCP(2,2) = - dsin(0.103d0)
!              rot_BCCtoHCP(3,2) =   0.d0

!              rot_BCCtoHCP(1,3) = + dsin(0.103d0) / dsqrt(2.d0)
!              rot_BCCtoHCP(2,3) = + dcos(0.103d0) / dsqrt(2.d0)
!              rot_BCCtoHCP(3,3) = + 1.d0 / dsqrt(2.d0)
              
              
!--------------------------------------------------------------------------
            ! ********************************************************************************!
            rot = MATMUL(rot_HCP,rot_BCCtoHCP) ! THIS DOES NOT SEEM RIGHT - incorrect ordering? !
            ! ********************************************************************************!
            call rotate_crys_vector(cst_BCC,rot,gst,nslip_BCC) ! cst, cmt are for BCC here.
            call rotate_crys_vector(cmt_BCC,rot,gmt,nslip_BCC)
            call rotate4d(ddsdde_4d,rot,de_gg)
            call tr4to2(de_gg,de_gg_2d)
         ! deniz - comment
         ! here:
         ! rot is the rotation matrix that transforms vectors as the crystal orientation
         ! cst, cmt are the slip directions and plane normals in the natural crystal frame
         ! gst(3,slip), gmt(3,slip) are slip directions and plane normals in the undeformed config of the grain           
!--------------------------------------------------
         endif
         
         s_gg_2d(:,:)=de_gg_2d(:,:)
         
         call lapack_invert(6,s_gg_2d,INFO)
         if (INFO.NE.0) then
            ! matrix inversion error
            write(*,*) 'inversion error for de_gg_2d:'
            CALL printMatrix(de_gg_2d,6,6, & 
                  'singular matrix de_gg_2d:',0)
            umat_Err = 7
            stop
         endif
         
         
         ! NOW, the model parameters (flow rule and hardening evolution params) 
         ! have been imported according to the phase and deformation type (tens/compression).
         ! pass them over to UMAT routine
                     
         ! deniz - comment
         ! here the schmidt tensors are formed, using the slip direction and plane normal vectors
         s_tan=0.0d0      
         do isys=1,nslip
            do j=1,3
             do i=1,3
                s_tan(i,j,isys)=gst(i,isys)*gmt(j,isys)
             enddo
            enddo
         enddo
            
         ! HERE WE ARE INITIALIZING TO THE INIT. VALUES OF SOFT/HARD BCC slip Systems--------
         ! ORGANIZE THIS IN A BETTER WAY

         ! traditionally,
         ! code assigns g0_bcc(1,2,4) to all BCC families,
         ! g1_bcc(1,2,4) to all with b1
         ! g2_bcc(1,2,4) to all with b2
                  
         ! IF THIS IS NOT THE FIRST STEP (NSTEP>1)    
         ! retrieve internal variables:
         ! slip resistances, 2nd PK, Fp
         ! and tot_gamma from previous step statev.

         ! calculate initial and current slip system resistances
         ! from intrinsic slip resistance and Hall-Petch contribution
         if(iCrystal==HCP) then

            delta_g_HCP(1:nslip)= statev(1:30)     ! HCP, delta_g - work hardening at each slip system         
            g_hp_famHCP(1:8,1)  = statev(305:312)  ! HCP, hall petch contrib., tension
            g_hp_famHCP(1:8,2)  = statev(313:320)  ! HCP, hall petch contrib., compression
            ! intrinsic slip resistance is already stored in: g0_bar_HCP
            ! calculate current slip resistance g_alpha_t = g0_bar + g_hp + delta_g
            CALL calc_SlipSysStrengthsHCP(g0_HCP,g_HCP,iPhase,g_hp_famHCP,delta_g_HCP)
            
            tot_gamma(1:nslip)=statev(31:60)   ! tot gamma
            tpk_1d_t(1:6) =statev(61:66)       ! 2nd PK stress - HCP
            Fp_t(1,1:3)   =statev(67:69)       ! Fp - HCP - stored ROW-MAJOR in STATE VARIABLES
            Fp_t(2,1:3)   =statev(70:72)       ! Fp - HCP - stored ROW-MAJOR in STATE VARIABLES
            Fp_t(3,1:3)   =statev(73:75)       ! Fp - HCP - stored ROW-MAJOR in STATE VARIABLES
            chi_t(1:nslip)=statev(221:250)     ! Chi_t - HCP
            
         elseif(iCrystal==BCC) then

            delta_g_BCC(1:nslip)= statev(76:123)   ! BCC, delta_g - work hardening at each slip system
            g_hp_famBCC(1:5,1)  = statev(337:341)  ! BCC, hall petch contrib., b1[-111]
            g_hp_famBCC(1:5,2)  = statev(342:346)  ! BCC, hall petch contrib., all others
            ! intrinsic slip resistance is already stored in: g0_bar_BCC
            ! calculate current slip resistance g_alpha_t = g0_bar + g_hp + delta_g
            CALL calc_SlipSysStrengthsBCC(g0_BCC,g_BCC,g_hp_famBCC,delta_g_BCC)
            
            tot_gamma(1:nslip)=statev(124:171) ! tot gamma
            tpk_1d_t(1:6)=statev(172:177)      ! 2nd PK stress - BCC
            Fp_t(1,1:3)  =statev(178:180)      ! Fp - BCC - stored ROW-MAJOR in STATE VARIABLES
            Fp_t(2,1:3)  =statev(181:183)      ! Fp - BCC - stored ROW-MAJOR in STATE VARIABLES
            Fp_t(3,1:3)  =statev(184:186)      ! Fp - BCC - stored ROW-MAJOR in STATE VARIABLES
            chi_t(1:nslip)=statev(251:298)     ! Chi_t - BCC
            
         endif
         
         if (iCrystal==HCP) then
            wpBackStress_t(1:4) = statev(378:381)      ! plastic work dissipated on the back-stress on the (1) basal and (3) prismatic planes
         endif

         ! Calculate lambda_beta_t, lattice plane incompatibilities
         ! lambda_beta's store the burgers vector area densities
         ! they will be used update_implicit() to calculate the
         ! GND contribution to hardening (See: Ghosh, Anahid, 2013)
         Nye(1,1:3) = statev(357:359)  ! read Nye from state vars. (row major)
         Nye(2,1:3) = statev(360:362)
         Nye(3,1:3) = statev(363:365)
         do isys=1,nslip
            ! burgers vector on the slip plane isys.
            b_edge(:) = MATMUL(Nye(:,:),gmt(:,isys))
            ! subtract the screw modes, 
            ! so that burgers() measures the edge modes only.
            b_edge(:) = b_edge(:) -   &
               DOT_PRODUCT(b_edge(:),gmt(:,isys))*gmt(:,isys)  &
               /DOT_PRODUCT(gmt(:,isys),gmt(:,isys))
            lambda_beta_t(isys) =   &
               SQRT(DOT_PRODUCT(b_edge,b_edge))  ! REF: Ghosh, Anahid, 2013
         enddo
         
         ! added by deniz
         ! calculate thermal deformation gradient using
         ! thermal expansion coeffs. a1,a2,a3
         ! difference between current and initial temp.
         ! rotation matrix, rot, calculated from euler angls, is an active transformation from iCrystal to spatial frame
         CALL euler_slip(euler,rot)
         CALL calc_F_Thermal(F_thrm_tau,F_thrm_tau_inv,det_thrm,  &
                             thrm_C(:,iCrystal),TEMP(2),TEMP_0,rot)
                             
!         if (iCrystal==HCP) then
!            write(*,*) 'HCP basal', gmt(:,1)
!            write(*,*) 'HCP sys 1 - b', gst(:,1)
!            write(*,*) 'HCP sys 2 - b', gst(:,2)
!            write(*,*) 'HCP sys 3 - b', gst(:,3)
!            write(*,*) 'schmidt factors'
!            write(*,*) 'HCP sys 1', s_tan(3,3,1)
!            write(*,*) 'HCP sys 2', s_tan(3,3,2)
!            write(*,*) 'HCP sys 3', s_tan(3,3,3)
!         endif
!         if (iCrystal==HCP) then
!            write(*,*) 'HCP basal', gmt(:,1)
!            write(*,*) 'HCP sys 4 - b', gst(:,4)
!            write(*,*) 'HCP sys 5 - b', gst(:,5)
!            write(*,*) 'HCP sys 6 - b', gst(:,6)
!            write(*,*) 'schmidt factors'
!            write(*,*) 'HCP sys 4', s_tan(3,3,4)
!            write(*,*) 'HCP sys 5', s_tan(3,3,5)
!            write(*,*) 'HCP sys 6', s_tan(3,3,6)
!         endif

!         if (iCrystal==BCC) then
!            write(*,*) 'BCC 101', gmt(:,5)
!            write(*,*) 'BCC sys 5 - b', gst(:,5)
!            write(*,*) 'BCC sys 6 - b', gst(:,6)
!         endif
         
         ! solve the material model, calculate
         ! 2nd PK, delta gamma, hardening and back-stress at tau
         call ssmatc3_therm(F_tau,Fp_t,F_thrm_tau,F_thrm_tau_inv,  &
                 TEMP(2),dtime,iCrystal,iPhase,gst,gmt,s_tan,nslip, &
                 tpk_1d_t,tpk_1d_tau,tau,de_gg_2d,s_gg_2d,  &
                 c_alpha_1d,tot_gamma,delta_gamma,dGamma_dTau, &
                 pnewdt1,umat_Err)
!        converged on 2nd PK at tau.
!        OUTPUTS:
!           tpk_1d_tau  - converged 2nd PK STRESS at tau
!           delta_gamma - converged slips for the time step
!           dd_g        - (in MAT.MODULE) increase in the work-hardening contribution delta_g(t) over this time step (tens/comp/soft/hard independent)
!           chi_tau     - (in MAT.MODULE) back stress at tau
!           umat_Err:
!                 0: success
!                 1: tau/g > 2, or flow rule produced NaN
!                 2: g is negative
!                 3: delGamma > 0.02 (algorithm assumes delGamma << 1)
!                 4: 1st level iteration did not convergence
!                 5: 2nd level iteration did not convergence
!                 6: back-stress iteration did not converge

         ! actually, fem_uel should decide how to change the time step in case of umat error
         if(umat_Err.ne.0)then  
            pnewdt=pnewdt1
            return
         endif

         ! now that we have 2nd PK tau (tpk_1d_tau)
         ! state variables, orientations and Fp at tau will be calculated
         call dstress_therm(F_tau,F_thrm_tau,F_thrm_tau_inv,det_thrm,  &
                            Fp_t,Fp_tau,Fe_tau,Fe_tau_inv, &
                            gst,gmt,s_tan,nslip,tpk_1d_tau,delta_gamma,  &
                            stress,cauchy_2d,strain_1d,crp_strain_1d,dw_p)
        
         ! deniz calculate thermoelastic deform. grad. at time TAU:
         F_thermoelastic_tau = MATMUL(Fe_tau,F_thrm_tau)
         F_thermoelastic_tau_inv = MATMUL(F_thrm_tau_inv,Fe_tau_inv)

         ! OUTPUT: ddsdde
         ! Kourosh says this has to be updated accord. to therm deform.
         call jacobian(F_t,F_tau,Fp_t,Fp_tau,s_tan,c_alpha_1d,  &
           tau,delta_gamma,dGamma_dTau,dtime,  &
           de_gg,tpk_1d_tau,g_alpha_tau, &
           ddsdde,ntens,nslip,tot_gamma,iCrystal,INFO)
           
         if (INFO.NE.0) then
            write(*,*) 'inversion error at jacobian()'
            umat_Err = 7
            pnewdt = 0.5
            return
         endif

!******************************************************************************

         ! UPDATE THE STATE VAR VECTOR with CALCULATED VARIABLES at time TAU 
         if(iCrystal==HCP) then
            statev(1:30)   = delta_g_HCP(1:nslip) + dd_g(1:nslip)           ! change in work-hardening delta_g(t) during this time step
            statev(31:60)  = tot_gamma(1:nslip)+dabs(delta_gamma(1:nslip))  ! update Total plastic slips
            statev(61:66)  = tpk_1d_tau(1:6)                                ! 2nd PK resovled at [????] configuration
            statev(67:69)  = Fp_tau(1,1:3)                                  ! Fp - (ROW MAJOR)
            statev(70:72)  = Fp_tau(2,1:3)                                  !
            statev(73:75)  = Fp_tau(3,1:3)                                  !
            statev(221:250)= chi_tau(1:nslip)                               ! back stresses 
         elseif(iCrystal==BCC) then
            statev(76:123) = delta_g_BCC(1:nslip) + dd_g(1:nslip)         ! change in work-hardening delta_g(t) during this time step
            statev(124:171)= tot_gamma(1:nslip)+dabs(delta_gamma(1:nslip))! update Total plastic slips
            statev(172:177)= tpk_1d_tau(1:6)                              ! 2nd PK resovled at [????] configuration
            statev(178:180)= Fp_tau(1,1:3)                                ! Fp - (ROW MAJOR)
            statev(181:183)= Fp_tau(2,1:3)                                !
            statev(184:186)= Fp_tau(3,1:3)                                !
            statev(251:298)= chi_tau(1:nslip)                             ! back stresses 
         endif
         
         ! additional output:
         statev(368) = SQRT(DOT_PRODUCT(delta_gamma(1:nslip),delta_gamma(1:nslip))) ! norm of delta_gamma
                  
         if (iCrystal==HCP) then
            statev(378:381)= wpBackStress_tau(1:4)                        ! plastic work on the back stress on (1) basal and (3) prism planes
            statev(382:385)= nucleationFactor(1:4)                        ! thermally activated void nucleation factor
         endif


         ! HCP - c axis, stresses on the basal plane
         if(iCrystal==HCP)then
            ! modified - deniz 
            sys_basal = 1 ! plane = basal
            cbasal(:)=MATMUL(TRANSPOSE(F_thermoelastic_tau_inv),  &
                                                  gmt(:,sys_basal))
            ! normalize
            cbasal(:) = cbasal(:)/SQRT(DOT_PRODUCT(cbasal,cbasal))

            ! traction on basal plane
            bas_tract(:)=MATMUL(cauchy_2d(:,:),cbasal(:))

            ! normal and shear components
            stress_normal_basal=DOT_PRODUCT(bas_tract,cbasal)
            stress_shearvec_basal(:) =   &
                     bas_tract(:)-stress_normal_basal*cbasal(:)
            stress_shear_basal =   &
         SQRT(DOT_PRODUCT(stress_shearvec_basal,stress_shearvec_basal))
     
            statev(352) = stress_normal_basal      ! save normal stress on basal
            statev(353) = stress_shear_basal       ! save resolved shear stress on basal
            statev(354:356) = cbasal(:)            ! save c-axis (basal normal) direction             
         endif

         ! save total crystal rotation for later
         ! grain orientation/texture calculations
         CALL rudcmp(Fe_tau,R_el_cry,U_el_cry)
         R_tot_cry = MATMUL(R_el_cry,rot)
         statev(369:371) = R_tot_cry(1,1:3)  ! save total iCrystal rotation to 369-377
         statev(372:374) = R_tot_cry(2,1:3)  ! (ROW-MAJOR)
         statev(375:377) = R_tot_cry(3,1:3)
         
         ! Homogenization w/ Rule of Mixtures
         ! sum up the stress, strain and tangent elasticity of the constituent phases
         creep_mix  = creep_mix  + wt(iCrystal)*crp_strain_1d ! Creep strain (Green-Lagrange)
         strain_mix = strain_mix + wt(iCrystal)*strain_1d     ! Total strain (Logarithmic)
         stress_mix = stress_mix + wt(iCrystal)*stress        ! Cauchy stress 
         ddsdde_mix = ddsdde_mix + wt(iCrystal)*ddsdde
         dw_p_mix    = dw_p_mix    + wt(iCrystal)*dw_p          ! plastic work 

         pnewdt=pnewdt1

      enddo ! loop over the crystals within a phase

      statev(211:216)=stress_mix(1:6)       ! Cauchy stress time tau
      statev(187:192)=strain_mix(1:6)       ! Logarithmic Strain time Tau
      statev(321:326)=creep_mix(1:6)        ! Creep strain (Green-Lagrange)
      statev(299) = w_p_t + dw_p_mix        ! plastic work       
      
      RETURN
      end
                
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!  CALCULATES THE INVERSE OF A 3*3 MATRIX

      subroutine matinv3(a,ai,det)
!      implicit double precision*8(a-h,o-z)
       implicit double precision (a-h,o-z)

!      include 'aba_param.inc'

      dimension a(3,3), ai(3,3)
!
      det=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)  &
           *a(3,3)+a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)  &
           *a(1,3)*a(2,2))
      ai(1,1) =  ( a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
      ai(1,2) = -( a(1,2)*a(3,3)-a(1,3)*a(3,2))/det
      ai(1,3) = -(-a(1,2)*a(2,3)+a(1,3)*a(2,2))/det
      ai(2,1) = -( a(2,1)*a(3,3)-a(2,3)*a(3,1))/det
      ai(2,2) =  ( a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
      ai(2,3) = -( a(1,1)*a(2,3)-a(1,3)*a(2,1))/det
      ai(3,1) =  ( a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
      ai(3,2) = -( a(1,1)*a(3,2)-a(1,2)*a(3,1))/det
      ai(3,3) =  ( a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
      return
      end

!-----------------------------------------------------

!    THIS SUBROUTINE CALCULATES THE TRANSFORMATION MATRIX
!----------------------------------------------------------
!     phi   - euler(1)
!     theta - euler(2)
!     omega - euler(3)
!
! deniz - comment
!
! output, tlgt,:
!   an active transformation matrix corresponding to the actual rotation due to euler angles {v}_i = TLGT_ij {v0}_j
!   as a (passive) coordinate transformation matrix, it transforms coordinates from iCrystal to global spatial frame. {v}_i = TLGT_ij {v_cry}_j
!---------------------------------------------------
      subroutine euler_slip(euler,tlg)   ! Kourosh: tlgt --- RESTORED to DEKA: tlg (OIM-convention) !
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),tlgt(3,3)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)

  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      ! tlg is a passive rotation correspt to euler angles
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct
      
      ! tlgt is an active rotation correspt to euler angles
      call trans(tlg,tlgt)

      return
      end   

!---------------------------------------------------------
!---------------------------------------------------------

   ! modified by deniz for thermal expansion

!   REFERING TO THE THESIS THIS SUBROUTINE CALCULATES A,T*tr,B,C
!   AND IT CALLS THE UPDATE IMPLICIT SUBROUTINE WHICH USES 
!   NEWTON RAPHSON FOR CALCULATING T* 
!-----------------------------------------------------------------
!  INPUTS
!        F_tau       - total   deform grad at tau
!        F_thrm_tau  - thermal deform grad at tau
!        Fp_t        - plastic deform grad at t
!        tpk_1d_t    - 2nd PK stress at t
!        s_tan       - Schmidt tensor at t, wrt plastically def conf
!        g0_BCC,g_BCC- initial and current slip systems resistances
!        g0_HCP,g_HCP
!
!  LOCALS
!         a_tan    - A (from thesis Pg 175) --- now C_thermelast_0tau
!         t_tr_1d  - T*tr (in thesis). zeroth order approximation to STRESS at tau.
!         b_alpha  - B (Pg 175)
!         c_alpha_1d - C (Pg 175)
!
!  OUTPUTS
!        tpk_1d_tau  - converged 2nd PK STRESS at tau
!        delta_gamma - converged slips for the time step
!        dd_g        - (in MAT.MODULE) increase in the work-hardening contribution delta_g(t) over this time step (tens/comp/soft/hard independent)
!----------------------------------------------------------------
      subroutine ssmatc3_therm(F_tau,Fp_t,F_thrm_tau,F_thrm_tau_inv,  &
                        TEMP,dtime,iCrystal,iPhase,gst,gmt,s_tan,nslip, &
                        tpk_1d_t,tpk_1d_tau,tau,de_gg_2d,s_gg_2d,  &
                        c_alpha_1d,tot_gamma,delta_gamma,dGamma_dTau, &
                        pnewdt1,umat_Err)
            
      use crystal_Titanium
      
      implicit none
      
      real(8), intent(in) :: TEMP,dtime              ! ambient temperature (*K)
      integer, intent(in) :: iCrystal                ! 1=HCP, 2=BCC
      integer, intent(in) :: iPhase                  ! 1=priA,2=traB
      real(8), intent(in) :: gst(3,nslip_max),gmt(3,nslip_max)
      real(8), intent(in) :: s_tan(3,3,nslip_max)
      integer, intent(in) :: nslip
      
      real(8), intent(in) :: F_tau(3,3),Fp_t(3,3)
      real(8), intent(in) :: F_thrm_tau(3,3)     ! thermal deformation gradient at tau
      real(8), intent(in) :: F_thrm_tau_inv(3,3) ! inv of thermal deformation gradient at tau
      real(8), intent(in) :: tpk_1d_t(6)
      real(8), intent(out):: tpk_1d_tau(6)
      real(8), intent(inout):: tau(nslip_max)
      real(8), intent(in) :: de_gg_2d(6,6),s_gg_2d(6,6)
      real(8), intent(out):: c_alpha_1d(6,nslip_max)
      real(8), intent(in) :: tot_gamma(nslip_max)
      real(8), intent(out):: delta_gamma(nslip_max)
      real(8), intent(out):: dGamma_dTau(nslip_max)  ! derivative of the flow rule wrt tau, for jacobian/tang. stiffness calculation
      
      real(8), intent(out):: pnewdt1
      integer, intent(out):: umat_Err
      
      ! local variables
      real(8) :: C_elast_0tau(3,3), C_thermelast_0tau(3,3)
      real(8) :: F_thrm_invtra_tau(3,3)
      
      
      real(8) :: fp_t_inv(3,3),fp_t_inv_t(3,3),F_tau_TR(3,3),det_fp_t
      real(8) :: C_tau(3,3),a1_tan(3,3),a_tan(3,3),a_tan_1d(6),  &
                 b_alpha(3,3,nslip_max),b_alpha_1d(6,nslip_max), &
                 t_tr_1d(6),ckrone_1d(6)
                 
      ! dummy
      integer :: i,j,isys
                 
      ! begin statements
           
      ! deniz - thermal vars
      C_elast_0tau(:,:) = 0.D0
      C_thermelast_0tau(:,:) = 0.D0
      F_thrm_invtra_tau = TRANSPOSE(F_thrm_tau_inv)

      call matinv3(Fp_t,fp_t_inv,det_fp_t) 

      if(det_fp_t.eq.0d0) then
        write(*,*) '0 divide'
        stop
      endif    

      fp_t_inv_t = TRANSPOSE(fp_t_inv)
      F_tau_TR = TRANSPOSE(F_tau)        
      C_tau = MATMUL(F_tau_TR,F_tau)     ! C_tau
      ! C_thermelast_0tau = Fp_t_INVTRA*F_tau_TRA*F_tau*Fp_t_inv
      C_thermelast_0tau = MATMUL(MATMUL(fp_t_inv_t,C_tau),fp_t_inv)  
      
      ! modified by deniz
      ! multiply from sides: F_thr_tau^-T ( * ) T_thr_tau^-1
      ! to include thermal deformation in calculation of the
      ! zeroth order approx to elastic deform tensor at tau
      C_elast_0tau =   &
       MATMUL(MATMUL(F_thrm_invtra_tau,C_thermelast_0tau),F_thrm_tau_inv)
           
      call tr2to1(C_elast_0tau,a_tan_1d)
      call cal_ckrone_1d(ckrone_1d)

      t_tr_1d = 0.d0
      t_tr_1d = 0.5d0*MATMUL(de_gg_2d(:,:),(a_tan_1d(:)-ckrone_1d(:)))
      ! comment - deniz
      ! Here,
      !     de_gg_2d is the elasticity matrix reduced to rank 2 in Voigt notation
      !     a_tan_1d is (in my notation) C_e0_tau = Fthr_tau^-T Fp_t^-T F_tau^T F_tau Fp_t^-1 Fthr_tau^-1, zeroth order approx to elastic deform tensor at tau
         
      ! now we got t_tr_1d. this is the zeroth order approx to STRESS, the first term in the residual
      
      do isys=1,nslip
         b_alpha(:,:,isys) =   &
            MATMUL(F_thrm_invtra_tau(:,:),C_thermelast_0tau(:,:))
         b_alpha(:,:,isys) =   &
            MATMUL(b_alpha(:,:,isys),s_tan(:,:,isys))
         b_alpha(:,:,isys) =   &
            MATMUL(b_alpha(:,:,isys),F_thrm_tau_inv(:,:))
            
         ! symmetrize -- factor 0.5 will be applied later (see below)
         b_alpha(:,:,isys) =  &
            b_alpha(:,:,isys) + TRANSPOSE(b_alpha(:,:,isys))
      enddo
      ! now, b_alpha, has the thermal deformation included

      call tr2to1crys(b_alpha,b_alpha_1d,nslip)

      do isys=1,nslip
        do i=1,6
          c_alpha_1d(i,isys)=0.d0
          do j=1,6
            c_alpha_1d(i,isys)=c_alpha_1d(i,isys)  &
            +0.5d0*de_gg_2d(i,j)*b_alpha_1d(j,isys)
          enddo
        enddo
      enddo      
      
      ! now we got c_alpha_1d(slips). this is the matrix multiplier of the second term in STRESS.
      ! below, it will be multiplied with trial slip rates in a summation, to give the second term in stress.
      ! below subroutine runs NR method and converges on the 2nd PK at tau: 
      ! outputs: tpk_1d_tau, delta_gamma, dd_g, chi_tau
      
      ! deniz - modified : added arguments: dGamma_dTau
      call update_implicit(TEMP,dtime,iCrystal,iPhase,gst,gmt,s_tan,nslip, &
                           t_tr_1d,tpk_1d_t,tpk_1d_tau,tau,s_gg_2d,  &
                           c_alpha_1d,tot_gamma,delta_gamma,dGamma_dTau, &
                           pnewdt1,umat_Err)

      return
      end




!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine sigmat(x,xx)    
!                                                                     
!  arranges into a symmetric array xx, the six commponents of a vector x
!                                                                     
       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension x(6),xx(3,3)
!
      xx(1,1)=x(1)                                                    
      xx(2,2)=x(2)                                                    
      xx(3,3)=x(3)                                                    
      xx(2,3)=x(6)                                                    
      xx(3,2)=x(6)                                                    
      xx(3,1)=x(5)                                                    
      xx(1,3)=x(5)                                                    
      xx(1,2)=x(4)                                                    
      xx(2,1)=x(4)                                                    
      return
      end




!---------------------------------------------------------

!   THIS USES NEWTON RAPHSON TO CALCULATE T*
! ---------------------------------------------
!      xjn_2d - Jacobian of Newton Raphson for T*
!      g_n    - Gn (Residual) (Pg 175)
!      xjn_gn - Jn_inv*[Gn] (Pg 175)
!
!  INPUTS
!       TEMP          - ambient temperature (*K)
!       lambda_beta_t - (from MODULE) burgers vectors on slip planes (due to Nye at time t)
!
!  OUTPUTS - converged values for time tau:
!       tpk_1d_tau  - 2nd PK at tau
!       delta_gamma - plastic slip incr.
!       dd_g        - (in MAT.MODULE) increase in the work-hardening contribution delta_g(t) over this time step (tens/comp/soft/hard independent)
!       chi_tau     - (in MAT.MODULE) back stress at time tau
!--------------------------------------------
      subroutine update_implicit(TEMP,dtime,iCrystal,iPhase,gst,gmt,s_tan,nslip, &
                                 t_tr_1d,tpk_1d_t,tpk_1d_tau,tau,s_gg_2d,  &
                                 c_alpha_1d,tot_gamma,delta_gamma,dGamma_dTau, &
                                 pnewdt1,umat_Err)

      ! options module - import UMAT options
      use options, only: ELASTIC_ONLY, PRISMATIC_ONLY, BASAL_ONLY, &
                         HARDENING_DISABLED, BACKSTRESS_DISABLED, GND_DISABLED, &
                         damage_enabled, &
                         DEBUG_RESID, DEBUG_UMAT, UMAT_halfStep, UMAT_stabilize, &
                         USE_POWER_LAW, p_umat, q_umat, g0_add
                         
      use material_Ti6242

      implicit none
            
      real(8), intent(in) :: TEMP,dtime              ! ambient temperature (*K)
      integer, intent(in) :: iCrystal                ! 1=HCP, 2=BCC
      integer, intent(in) :: iPhase                  ! 1=priA,2=traB
      real(8), intent(in) :: gst(3,nslip_max),gmt(3,nslip_max)
      real(8), intent(in) :: s_tan(3,3,nslip_max)
      integer, intent(in) :: nslip
      
      real(8), intent(in) :: t_tr_1d(6)              ! zeroth order (trial) stress for time tau
      real(8), intent(in) :: tpk_1d_t(6)
      real(8), intent(out):: tpk_1d_tau(6)
      real(8), intent(inout):: tau(nslip_max)
      real(8), intent(in) :: s_gg_2d(6,6)
      real(8), intent(in) :: c_alpha_1d(6,nslip_max)
      real(8), intent(in) :: tot_gamma(nslip_max)
      real(8), intent(out):: delta_gamma(nslip_max)
      real(8), intent(out):: dGamma_dTau(nslip_max)  ! derivative of the flow rule wrt tau, for jacobian/tang. stiffness calculation
      
      real(8), intent(out):: pnewdt1
      integer, intent(out):: umat_Err
      
                         
      ! locals
      
      real(8), parameter :: g0=8.0d0
      
      ! LAPACK error flag:
      integer :: INFO
      ! residual calc error flag
      integer :: iatg

      real(8) :: g_n(6),xjn_2d(6,6),xjn_gn(6)
      real(8) :: xjn_2d_inv(6,6),xjn_2d_inv_temp(6,6)
      real(8) :: xhard(nslip_max,nslip_max)
      real(8) :: xhard_predictor(nslip_max,nslip_max),xhard_corrector(nslip_max,nslip_max)

      integer :: softHard(nslip_max)
      integer :: tensComp(nslip_max)
      real(8) :: dd_g_prev(nslip_max),dd_g_predictor(nslip_max)
      real(8) :: g_tau(nslip_max),g_tau_predictor(nslip_max)
      real(8) :: g_t(nslip_max),g0_bar(nslip_max)
      
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)     ! 1... CRSS, 1000... gamma dot, 2000... g_sat
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT     
      ! ---------- GND hardening variables ---------- !
      ! see REF Anahid, Ghosh 2013
      real(8) :: kGND                ! multiplier for GND hardening
      ! GND hardening term parameters
      real(8), parameter :: k0=2.d0
      real(8), parameter :: alpha_hat=1.d0/3.d0
      real(8), parameter :: G=48.d0    ! GPa
      
      real(8) :: hardMOB, hardGND
      
      ! PNEWDT Multipliers for UMAT errors and adaptive increase
      real(8), parameter :: PNEWDT_DECR_MULTIPLIER = 0.5
      real(8), parameter :: PNEWDT_INCR_MULTIPLIER = 1.5
      
      integer :: i,j,isys,jsys,ni
      integer :: nIter_g,nIter_S
      logical :: tensComp_Converged, softHard_Converged
      real(8), parameter :: g_tol=1.0d0
      real(8), parameter :: xtol=1.d-10
      real(8), parameter :: chi_tol=1.d0
      real(8), parameter :: gamma_tol=0.02d0
      integer, parameter :: iter_max_S=300
      integer, parameter :: iter_max_g=30
      real(8) :: g_max, chi_max, chi_prev, tensComp_prev(nslip_max), softHard_prev(nslip_max)
      real(8) :: rnorm,xeta
      
      iatg = 0
      tpk_1d_tau(1:6)=tpk_1d_t(1:6)

      ! for the initial iter loop, we assume the hardening params and back stress stay constant at time t values
      dd_g(1:nslip)    = 0.d0       ! increment in work-hardening contribution during this time step t -> tau
      chi_tau(1:nslip) = chi_t(1:nslip)
      if (iCrystal==HCP) wpBackStress_tau(1:4) = wpBackStress_t(1:4)

      umat_Err=0

!     rnorm, g_max defined to enter the iteration loops
      g_max=2*g_tol
      chi_max=2*chi_tol
      softHard_Converged = .true.
      tensComp_Converged = .true.
      
      nIter_g = 0
      do while( ((g_max > g_tol.or.chi_max > chi_tol).or. &
                 (.not.tensComp_Converged.or..not.softHard_Converged)).and. &
                 nIter_g < iter_max_g )
                 
         if (UMAT_stabilize.and.nIter_g > 10) then  ! if not converging, try to stabilize by damping/averaging
            dd_g = 0.5d0 * (dd_g + dd_g_prev)
         endif
               
         if (iCrystal==HCP) then
            ! determine tension/compression state      
            CALL getTensionCompressionState(tensComp,tpk_1d_t,gmt,nslip)   ! STRESS at TIME t ! ! REVERT THIS
            ! assign CRSSs
            do isys = 1, nslip
               g_tau(isys) = g_HCP(isys,tensComp(isys)) + dd_g(isys) ! at time tau
               g_t(isys) = g_HCP(isys,tensComp(isys))                ! at time t
               g0_bar(isys) = g0_bar_HCP(isys,tensComp(isys),iPhase) ! intrinsic slip resistance (without HP contrib)
            enddo
         elseif(iCrystal==BCC) then
            ! determine if slip is along soft/hard direction
            CALL getSoftHardDirection(softHard,tpk_1d_t,s_tan,nslip)
            ! assign CRSSs
            do isys = 1, nslip
               g_tau(isys) = g_BCC(isys,softHard(isys)) + dd_g(isys) ! at time tau
               g_t(isys)   = g_BCC(isys,softHard(isys))                ! at time t
               g0_bar(isys) = g0_bar_BCC(isys,softHard(isys))
            enddo                        
         endif
         
         ! entering slip-iteration (inner loop)
         ! assuming CRSS/back-stress: chi_tau(1:nslip) and g_tau(1:nslip)
         ! record the assumed back-stresses for debugging purposes:
         if(DEBUG_UMAT.and.iCrystal==HCP) then
            realData(nIter_g*6+1:nIter_g*6+6) = g_tau(1:6)        ! ! ! ! stabilization attempt
         endif
            
         nIter_S=1
         rnorm=2*xtol
         do while(nIter_S < iter_max_S)  ! convergence criteria below, right after the residual calculation
         
            !     Calculating the residual
            !       tpk_1d_tau  (in)  the proposed/trial STRESS at tau.
            !       t_tr_1d     (in)  zeroth order approximation to the stress
            !       TEMP        (in)  ambient temp
            !       tau         (out) resolved shear stress due to trial stress
            !       g_n         (out) residual = trial stress - stress calculated via slips: (t_tr_1d + slips*c_alpha_1d)
            !       xjn_2d      (out) Jacobian for Newton Raphson            
            call stress_residue(TEMP,dtime,iCrystal,s_tan,nslip, &
                                0.5d0*(g_t+g_tau),       &              ! ! ! ! stabilization attempt: use g at midstep
                                t_tr_1d,tpk_1d_tau,tau,s_gg_2d,  &
                                g_n,xjn_2d,1,c_alpha_1d,         &
                                tot_gamma,delta_gamma,dGamma_dTau, &
                                iatg)  
                                
            ! check UMAT errors - iatg:
            !  0: success
            !  1: tau/g > 2  --> should I tolerate this during iterations, or check only after convergence?
            !  2: g is negative
            !  3: flow rule has produced NaN            
            ! * should I tolerate this during iterations, or check only after convergence?
            ! * I can try to recover from a tau/g > 2, or avoid the steep the region
            !  - if tau/g > 1 for any system during the first few iterations, 
            !  take a half step back, and recalculate residual
            !     if(nIter_S <= 2 .and. ANY(tau(:)/g_tau(:) > 1.d0) ) then
            !        tpk_1d_tau = tpk_1d_tau + 0.5d0*xeta*xjn_gn
            !        cycle
            !     endif
            if(iatg.ne.0) goto 701  ! error! jump to error handler and exit. 

            rnorm = dsqrt(DOT_PRODUCT(g_n,g_n))
            
            if (rnorm <= xtol) exit  ! stress converged, exit NR loop
            
            ! residual still large. invert jacobian and update stress.
            
            ! tweak step size.
            xeta=1.d0
            ! if this is the first iteration, 
            ! dampen the NR update to avoid over-shooting the flow rule
            ! (avoids NaNs and tau/g>2 errors for large strain increments)
            if (nIter_S == 1.and.UMAT_halfStep) xeta=0.50d0
            
            do isys=1,nslip
               if(dabs(delta_gamma(isys)).gt.gamma_tol) xeta=0.25d0
            enddo
            
            ! NEW:
            ! use LAPACK to invert matrix:
            xjn_2d_inv(:,:) = xjn_2d(:,:)
            call lapack_invert(6,xjn_2d_inv,INFO)
            if (INFO.NE.0) then
               ! matrix inversion error
               write(*,*) 'singular umat jac xjn_2d. g_n(1:6)',g_n(1:6)
               pnewdt1=PNEWDT_DECR_MULTIPLIER
               umat_Err = 7
               return
            endif
            
            ! NR-update of the trial/proposed 2nd PK STRESS at tau
            xjn_gn = MATMUL(xjn_2d_inv,g_n)
            tpk_1d_tau = tpk_1d_tau - xeta*xjn_gn
            
            nIter_S=nIter_S+1
            
         enddo ! N-R loop for STRESS, assuming constant: g_alpha, chi, tens/Comp, soft/hard
         
         if (DEBUG_UMAT.and.iCrystal==HCP) then     ! update UMAT CPU counter for debugging/optiization purposes: 1000: gammadot, 
            intData(nIter_g+1) = nIter_S            ! number of delta-gamma iterations within this back-stress iteration
            intData(4999) = intData(4999) + nIter_S ! total number of delta-gamma iterations
            intData(5000) = nIter_g+1               ! total number of back-stress iterations
            realData(nIter_g*6+1001:nIter_g*6+1006) = delta_gamma(1:6)/dtime

         endif
         
         if(rnorm.gt.xtol)then
            write(*,*)'The First Level (STRESS) Iteration did not converge'
            write(*,*)'rnorm=',rnorm
            umat_Err=4
            pnewdt1=PNEWDT_DECR_MULTIPLIER
            return
         endif

!        Hardness Iteration
!        updates the hardening rate according to the converged stress and slip

         ! INPUTS
         !     delta_gamma
         !     g_tau
         ! OUTPUTS
         !     xhard(a,b)  - hardening rate due to self & latent: h_ab
         !     kGND        - GND hardening multiplier
         
         !---PREDICTOR
         call make_hard(g_tau,g_t,delta_gamma,tpk_1d_t,TEMP,  &
         dtime,xhard_predictor,tau,nslip,tot_gamma,iCrystal,iPhase,softHard,tensComp)
         
         do isys=1,nslip
            
            ! calc kGND for this slip system
            kGND = 0.5*k0*(alpha_hat**2.0)*(G**2.0)*burgers(isys)  &
                      /(20.0 + g_tau(isys)-g0_bar(isys))
            kGND = 1.D3*kGND  ! transform to proper units. G: GPa, b: nm, g_alpha: MPa, lambda_beta: 1/um.
            if(GND_DISABLED) kGND=0.D0 ! GND disable switch for debug purposes
            
            dd_g_predictor(isys) = 0.d0
            do jsys=1,nslip

               dd_g_predictor(isys)=dd_g_predictor(isys)  &
                  + xhard_predictor(isys,jsys)*dabs(delta_gamma(jsys)) &  ! xhard(a,b) hardening rate due to self & latent: h_ab
                  + kGND*lambda_beta_t(jsys)*dabs(delta_gamma(jsys))   ! GND hardening REF Anahid, Ghosh 2013
                              
            enddo
            
            ! to disable hardening:
            if(HARDENING_DISABLED) dd_g_predictor(isys)=0.d0

         enddo
         
         !---CORRECTOR
         g_tau_predictor(1:nslip) = g_t(1:nslip) + dd_g_predictor(1:nslip)
                  
         call make_hard(g_tau_predictor,g_t,delta_gamma,tpk_1d_t,TEMP,  &
         dtime,xhard_corrector,tau,nslip,tot_gamma,iCrystal,iPhase,softHard,tensComp)
         
         xhard = 0.5d0*(xhard_predictor + xhard_corrector)                    ! ! ! ! ! stabilization attempt: use a predictor-corrector method for converging on g

!         call make_damage(g_tau,g_t,delta_gamma,wpBackStress_t, &
!         dtime,xhard,tau,nslip,tot_gamma,iCrystal,iPhase,softHard,tensComp)
         
         ! ---------- GND HARDENING TERM ----------- !
         ! deniz - fix this
         ! G actually depends on the dislocation line direction
         ! for transv. isotropic, shear modulus sub-matrix is NOT diagonal
         ! for non-principal coord frames.
         ! G = 1.0/3.0*(de_gg_2d(4,4)+de_gg_2d(5,5)+de_gg_2d(6,6))*1.D6 ! in Pa         
         ! ----------------------------------------- !
         g_max=0.d0
         dd_g_prev(1:nslip) = dd_g(1:nslip)
         do isys=1,nslip
            
            ! calc kGND for this slip system
            kGND = 0.5*k0*(alpha_hat**2.0)*(G**2.0)*burgers(isys)  &
                      /(20.0 + g_tau(isys)-g0_bar(isys))
            kGND = 1.D3*kGND  ! transform to proper units. G: GPa, b: nm, g_alpha: MPa, lambda_beta: 1/um.
            if(GND_DISABLED) kGND=0.D0 ! GND disable switch for debug purposes
            
            dd_g(isys) = 0.d0
            do jsys=1,nslip

               dd_g(isys)=dd_g(isys)  &
                  + xhard(isys,jsys)*dabs(delta_gamma(jsys)) &  ! xhard(a,b) hardening rate due to self & latent: h_ab
                  + kGND*lambda_beta_t(jsys)*dabs(delta_gamma(jsys))   ! GND hardening REF Anahid, Ghosh 2013
                              
            enddo
            
            ! to disable hardening:
            if(HARDENING_DISABLED) dd_g(isys)=0.d0

            if(dabs(dd_g(isys)-dd_g_prev(isys)) > g_max) &
               g_max = dabs(dd_g(isys)-dd_g_prev(isys))
         enddo
         
         ! update hardening parameters and check if tens/comp/soft/hard has changed.
         if (iCrystal==HCP) then
            ! update the CRSSs with the calculated hardening increment
            do isys = 1, nslip
               g_tau(isys) = g_HCP(isys,tensComp(isys)) + dd_g(isys)
               g0_bar(isys) = g0_bar_HCP(isys,tensComp(isys),iPhase)
            enddo
            tensComp_prev(1:nslip) = tensComp(1:nslip)
            ! determine tension/compression state      
            CALL getTensionCompressionState(tensComp,tpk_1d_t,gmt,nslip)         ! STRESS at TIME t ! ! REVERT THIS
            tensComp_Converged = (ALL(tensComp_prev(1:nslip) == tensComp(1:nslip)))
         elseif(iCrystal==BCC) then
            ! update the CRSSs with the calculated hardening increment
            do isys = 1, nslip
               g_tau(isys) = g_BCC(isys,softHard(isys)) + dd_g(isys)
               g0_bar(isys) = g0_bar_BCC(isys,softHard(isys))
            enddo           
            softHard_prev(1:nslip) = softHard(1:nslip)
            ! determine if slip is along soft/hard direction
            CALL getSoftHardDirection(softHard,tpk_1d_t,s_tan,nslip)
            softHard_Converged = (ALL(softHard_prev(1:nslip) == softHard(1:nslip)))
         endif
                  
!        Back Stress Iteration
         chi_max=0.d0
         do isys=1,nslip
            chi_prev=chi_tau(isys)
            
            chi_tau(isys)=chi_t(isys)+xc*delta_gamma(isys)
            chi_tau(isys)=chi_tau(isys)/(1.d0+xd*  &
                          dabs(delta_gamma(isys)))

            if (BACKSTRESS_DISABLED) chi_tau(isys) = 0.d0
            
            if(dabs(chi_tau(isys)-chi_prev).gt.chi_max)chi_max=  &
              dabs(chi_tau(isys)-chi_prev)
         enddo
         
         ! update energy dissipated on back stress
         if (iCrystal==HCP) then
            wpBackStress_tau(1) = wpBackStress_t(1) & 
                        + delta_gamma(1)*(chi_t(1) + chi_tau(1)) / 2.d0 &
                        + delta_gamma(2)*(chi_t(2) + chi_tau(2)) / 2.d0 &
                        + delta_gamma(3)*(chi_t(3) + chi_tau(3)) / 2.d0
            wpBackStress_tau(2) = wpBackStress_t(2) & 
                        + delta_gamma(4)*(chi_t(4) + chi_tau(4)) / 2.d0
            wpBackStress_tau(3) = wpBackStress_t(3) & 
                        + delta_gamma(5)*(chi_t(5) + chi_tau(5)) / 2.d0
            wpBackStress_tau(4) = wpBackStress_t(4) & 
                        + delta_gamma(6)*(chi_t(6) + chi_tau(6)) / 2.d0
         endif
                  
         nIter_g=nIter_g+1

      enddo ! do while((g_max.gt.g_tol.or.chi_max.gt....
      
      ! update converged CRSS value
      if(DEBUG_UMAT.and.iCrystal==HCP) then
         realData(nIter_g*6+1:nIter_g*6+6) = g_tau(1:6)
      endif
      
      call stress_residue(TEMP,dtime,iCrystal,s_tan,nslip, &
                          0.5d0*(g_t+g_tau),       &           ! ! ! ! stabilization attempt: use g at midstep
                          t_tr_1d,tpk_1d_tau,tau,s_gg_2d,  &
                          g_n,xjn_2d,0,c_alpha_1d,         &
                          tot_gamma,delta_gamma,dGamma_dTau, &
                          iatg)
      ! umat_Err:
      !  0: success
      !  1: tau/g > 2
      !  2: g is negative
      !  3: delGamma > 0.02 (algorithm assumes its much smaller than 1)
      !  4: 1st level iteration did not convergence
      !  5: 2nd level iteration did not convergence
      !  6: back-stress iteration did not converge
      !(-): matrix inversion/numerical error
      
701   continue ! error handling

      if(iatg.eq.1)then
         umat_Err=1
         pnewdt1=PNEWDT_DECR_MULTIPLIER
         ! write(*,*)'tau/g exceed 2 @elem',curElemUMAT
         return
      elseif (iatg.eq.2)then
         umat_Err=2
         pnewdt1=PNEWDT_DECR_MULTIPLIER
         write(*,*)'g is negative, umat line 2129 @elem',curElemUMAT
         return
      elseif (iatg.eq.3)then
         umat_Err=1
         pnewdt1=PNEWDT_DECR_MULTIPLIER
         write(*,*)'flow rule has produced a NaN (stress_residue) @elem',curElemUMAT
         return
      endif

      ! residual calculations were successful. 
      ! check other potential errors:

      ni=0      
      do isys=1,nslip
         if(dabs(delta_gamma(isys)).gt.gamma_tol) then
            ! the mathematical assumption, delGamma << 1, is not valid
            !write(*,*)'assumption delGamma << 1, is not valid',curElemUMAT
            umat_Err = 3
            pnewdt1=PNEWDT_DECR_MULTIPLIER
            return 
         elseif(dabs(delta_gamma(isys)).le.0.01d0*gamma_tol) then
            ni=ni+1
         endif
      enddo
      if(ni.eq.nslip)pnewdt1=PNEWDT_INCR_MULTIPLIER
      
      if(g_max.gt.g_tol)then
         !write(*,*)'Second Level Iteration did not converge',curElemUMAT
         !write(*,*)'g_max=',g_max,'# iterations:',nIter_g
         umat_Err=5
         pnewdt1=PNEWDT_DECR_MULTIPLIER
         
         return 
      endif

      if(chi_max.gt.chi_tol)then
         write(*,*)'Back Stress iteration did not converge',curElemUMAT
         write(*,*)'chi_max=',chi_max,'# iterations:',nIter_g
         umat_Err=6
         pnewdt1=PNEWDT_DECR_MULTIPLIER
         return
      endif
      
      return
      end  


!-----------------------------------------------------------------------
      subroutine make_hard(g_alpha_tau,g_alpha_t,delta_gamma,tpk_1d,TEMP,  &
        dtime,xhard,tau,nslip,tot_gamma,iCrystal,iPhase,softHard,tensComp)
        
      use material_Ti6242
      
      use options, only: DEBUG_UMAT, UMAT_stabilize, damage_gsat_enabled
      
      implicit none
      
      integer, intent(in) :: softHard(nslip_max)   ! slip is along 1: soft direction, 2: hard direction
      integer, intent(in) :: tensComp(nslip_max)   ! stress state is 1: tension, 2: compression
      integer, intent(in) :: nslip
      real(8), intent(in) :: g_alpha_tau(nslip_max)
      real(8), intent(in) :: g_alpha_t(nslip_max)
      real(8), intent(in) :: delta_gamma(nslip_max)
      real(8), intent(in) :: tpk_1d(6)             ! 2nd PK stress for calculating damage progression/void nucleation rate
      real(8), intent(in) :: TEMP                  ! ambient temperature
      real(8), intent(in) :: dtime
      real(8), intent(out):: xhard(nslip_max,nslip_max)
      real(8), intent(in) :: tot_gamma(nslip_max)
      real(8), intent(in) :: tau(nslip_max)
      integer, intent(in) :: iCrystal,iPhase
      
      real(8) :: qab(nslip_max,nslip_max),hbeta(nslip_max),gamma_dot(nslip_max)
      real(8) :: tg,xharg,xhsech,xtemp
      integer :: i,j,isys_a,isys_b,isys,iPlane
      ! values of parameters in tens/comp, soft/hard, HCP/BCC etc.
      real(8) :: xs_bar(nslip_max),xh0(nslip_max),xn(nslip_max),xr(nslip_max)
      real(8) :: xs(nslip_max)
      real(8) :: xg_0,xg_inf,xh_0,xh_inf
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT   
      
      ! deniz - thermally activated void nucleation/damage
      real(8), parameter :: PI=3.14159265359d0 ! pi
      real(8), parameter :: kB=1.3806488d-23   ! boltzmann constant (J/K)
      real(8), parameter :: surfaceEnergy=2.d0 ! J/m^2
      real(8), parameter :: MPatoSI = 1.0d6    ! x10^6
      real(8) :: hydrostaticTension, p, R_crit

      xhard(1:nslip,1:nslip)=0.d0

      if(dabs(dtime) <= 1d-10) return
      
      
      do isys_a=1,nslip
         do isys_b=1,nslip
            if(isys_a.eq.isys_b)then
             qab(isys_a,isys_b)= hardeningSelf  ! self hardening
            else
             qab(isys_a,isys_b)= hardeningLatent  ! latent hardening
            endif
         enddo
      enddo
      
      if(iCrystal==BCC)then
         
         tg=0.d0
         do i=1,nslip
            tg=tg+tot_gamma(i)+dabs(delta_gamma(i))
         enddo    
     
         do isys_b=1,nslip
            ! CHECK values
            xh_0   =   xh_0_BCC(isys_b,softHard(isys_b))
            xh_inf = xh_inf_BCC(isys_b,softHard(isys_b))
            xg_0   =   xg_0_BCC(isys_b,softHard(isys_b))
            xg_inf = xg_inf_BCC(isys_b,softHard(isys_b))
         
            ! CHECK formulation
            xharg=(xh_0-xh_inf)/(xg_inf-xg_0)*tg
            xhsech=1.d0/dcosh(xharg)
            hbeta(isys_b)=xh_inf+(xhsech**2.d0)*(xh_0-xh_inf)
         enddo
         
         ! EVOLUTION OF HARDENING RATES. xhard here = q_AB*h_B in the papers         
         do isys_a=1,nslip
            do isys_b=1,nslip        

            
               ! CHECK FORMULATION
               xhard(isys_a,isys_b)=qab(isys_a,isys_b)*hbeta(isys_b)
            enddo
         enddo        
         
      endif
      
      if(iCrystal==HCP)then   ! HCP

         do isys=1,nslip
            gamma_dot(isys)=delta_gamma(isys)/dtime
         enddo
       
         do isys=1,nslip
            ! assign hardening parameter values according to tens/Comp/priA/transB
            xs_bar(isys) = xs_bar_HCP(isys,tensComp(isys),iPhase)
            xh0(isys)    = xh0_HCP(isys,tensComp(isys),iPhase)
            xn(isys)     = xn_HCP(isys,tensComp(isys),iPhase)  
            xr(isys)     = xr_HCP(isys,tensComp(isys),iPhase)  
         enddo
         
         ! ----- deniz damage -------
         hydrostaticTension = (tpk_1d(1)+tpk_1d(2)+tpk_1d(3))/3.d0
         do iPlane = 1,4
            p = 1.0d4*( wpBackStress_t(iPlane) + abs(hydrostaticTension) ) * MPatoSI
            if (p /= 0.d0) then
               R_crit = 2.d0*surfaceEnergy/p
               nucleationFactor(iPlane) = dexp( -(+surfaceEnergy*4.d0*PI*R_crit**2.d0 &
                                                     -p*4.d0/3.d0*PI*R_crit**3.d0)/ (kB*TEMP))
            else
               nucleationFactor(iPlane) = 0.d0
            endif
            ! apply the damage coefficient to the saturation stress
         enddo
         
         if (damage_gsat_enabled) then            
            xs_bar(1:3) = xs_bar(1:3)*(1.d0 - nucleationFactor(1))   ! basal
            xs_bar(4) = xs_bar(4)*(1.d0 - nucleationFactor(2))       ! prismatics
            xs_bar(5) = xs_bar(5)*(1.d0 - nucleationFactor(3))       ! prismatics
            xs_bar(6) = xs_bar(6)*(1.d0 - nucleationFactor(4))       ! prismatics
         endif
         
         ! HCP - evolution of slip system deformation resistance
         do isys_b=1,nslip
            hbeta(isys_b)=0.d0
            if(gamma_dot(isys_b).ne.0.d0)then
               
               xs(isys_b) = xs_bar(isys_b)* &
                   (dabs(gamma_dot(isys_b)/a_dot(isys_b)))**xn(isys_b)
               
!  use g_t for pulling g towards g_sat
!               xtemp = 1.d0 - (g_alpha_t(isys_b))/xs(isys_b) 
!  use g_tau for pulling g towards g_sat
!               xtemp = 1.d0 - (g_alpha_tau(isys_b))/xs(isys_b)
!  use 0.5(g_t + g_tau) for pulling g towards g_sat (mid-point method)
               xtemp = 1.d0 - 0.5d0*(g_alpha_t(isys_b)+g_alpha_tau(isys_b))/xs(isys_b)    ! ! ! ! stabilization attempt: use g at midstep --- all these changes do not completely solve the problem. still got to use hbeta-damping below..

               ! saturation stress stabilizer !
!               if (UMAT_stabilize.and.dabs(xtemp) < 0.03d0) then   ! g is too close to g_sat !
!                  xh0(isys_b) = xh0(isys_b) * (dabs(xtemp) / 0.03d0)
!               endif

               hbeta(isys_b) = &
                     xh0(isys_b)*((dabs(xtemp))**xr(isys_b))*dsign(1.d0,xtemp)
            endif
         enddo

         ! set up the h_ab matrix
         do isys_a=1,nslip
            do isys_b=1,nslip
               if(dabs(gamma_dot(isys_b)).gt.0.d0)then
                  xhard(isys_a,isys_b) = qab(isys_a,isys_b)*hbeta(isys_b)
               else
                  xhard(isys_a,isys_b) = 0.d0
               endif
            enddo
         enddo
         
         ! record the saturation stresses for debugging purposes - REVERT THIS ! ! !
         if (DEBUG_UMAT.and.iCrystal==HCP) then
            realData((intData(5000)-1)*6+2001:(intData(5000)-1)*6+2006) = xs(1:6)
         endif
         
         
      endif ! HCP

      end subroutine

      subroutine make_damage(g_alpha_tau,g_alpha_t,delta_gamma,wpBackStress,  &
        dtime,xhard,tau,nslip,tot_gamma,iCrystal,iPhase,softHard,tensComp)
        
      use material_Ti6242
      use options, only: damage_enabled
      
      implicit none
      
      integer, intent(in) :: softHard(nslip_max)   ! slip is along 1: soft direction, 2: hard direction
      integer, intent(in) :: tensComp(nslip_max)   ! stress state is 1: tension, 2: compression
      integer, intent(in) :: nslip
      real(8), intent(in) :: g_alpha_tau(nslip_max)
      real(8), intent(in) :: g_alpha_t(nslip_max)
      real(8), intent(in) :: delta_gamma(nslip_max)
      real(8), intent(in) :: wpBackStress(4)
      real(8), intent(in) :: dtime
      real(8), intent(out):: xhard(nslip_max,nslip_max)
      real(8), intent(in) :: tot_gamma(nslip_max)
      real(8), intent(in) :: tau(nslip_max)
      integer, intent(in) :: iCrystal,iPhase
      
      real(8) :: slipSystemDamage(4)
      real(8) :: gamma_dot(nslip_max)
      integer :: i,j,isys_a,isys_b,isys
      ! values of parameters in tens/comp, soft/hard, HCP/BCC etc.
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT   

      if(dabs(dtime) <= 1d-10) return
      
      do isys=1,nslip
         gamma_dot(isys)=delta_gamma(isys)/dtime
      enddo
      
      if(iCrystal==HCP)then   ! HCP
      
         !-- plastic slip based
         !slipSystemDamage(1) = (dabs(tot_gamma(1))/slipSystemLife + &
         !                       dabs(tot_gamma(2))/slipSystemLife + &
         !                       dabs(tot_gamma(3))/slipSystemLife)**slipSystemExhaustionExp
         !-- dislocation interaction energy based
         slipSystemDamage(1) = (abs(wpBackStress(1))/slipSystemLife)**slipSystemExhaustionExp

         !slipSystemDamage(1) = MAX(slipSystemDamage(1),1.d0)   ! can do up to 1
         
         if (damage_enabled) then
         
            do isys=1,3
            xhard(isys,isys) = xhard(isys,isys) &
                     - g_alpha_tau(isys) * slipSystemDamage(1) * &
                     (1.d0 / slipSystemDecohesionLife)
            enddo
            
         endif
         
         do isys=4,6
            !-- plastic slip based
            !slipSystemDamage(isys-2) = (dabs(tot_gamma(isys))/slipSystemLife)**slipSystemExhaustionExp
            !-- dislocation interaction energy based
            slipSystemDamage(isys-2) = (wpBackStress(isys-2)/slipSystemLife)**slipSystemExhaustionExp
            
            !slipSystemDamage(isys-2) = MAX(slipSystemDamage(isys-2),1.d0)   ! can do up to 1
            
            if (damage_enabled) then
            
               xhard(isys,isys) = xhard(isys,isys) &
                        - g_alpha_tau(isys) * slipSystemDamage(isys-2) * &
                        (1.d0 / slipSystemDecohesionLife)
                        
            endif
            
         enddo

         
      endif ! HCP

      end subroutine



      subroutine rclear (a,max)                                       
!                                                                     
!  fills a real array a with zeros
!                                                                     
       implicit double precision (a-h,o-z)                             

      dimension a(max)
!
      do i=1,max                                                 
      a(i) = 0.d0                                                     
      enddo
      return
      end         
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine prod(a,b,c)                                          
!                                                                     
!  computes the matrix product c=a*b  all matricies are 3x3
!                                                                     
      implicit double precision (a-h,o-z)                             

      dimension a(3,3),b(3,3),c(3,3)                                  
!                                                                     
      do 200 j=1,3
         do 200 i=1,3
            s=0.0
            do 100 k=1,3
               s=s+a(i,k)*b(k,j)
 100  continue
        c(i,j)=s
 200  continue
      return
      end

!---------------------------------------------------------------------
!      THIS SUBROUTINE CALCULATES THE VALUE OF Fp,NORMALIZES IT AND
!      CALCULATES THE CAUCHY STRESS AND THE NEW TEXTURE
!      * modified by deniz for thermal deformations
!-----------------------------------------------------------------
!  INPUTS
!     F_tau      - total deform grad at tau
!     F_thrm_tau  - thermal deformation gradient at tau
!     Fp_t        - Fp at t
!     s_tan       - Schmidt tensor at t
!     delta_gamma - slips during time step
!     tpk_1d      - converged 2nd PK at tau (in Voigt form)
!     gst, gmt    - slip direction and plane normal vectors at t
!
!  OUTPUTS
!     Fp_tau      - plastic deform grad at tau
!     Fe_tau      - elastic deform grad at tau
!     stress      - Cauchy stress at tau (in Voigt form)
!     cauchy_2d   - Cauchy stress at tau (in 3x3 matrix form)
!     strain_1d   - Logarithmic strain at tau (in Viogt form)
!     dw_p        - plastic work density over this increment within this phase      
!     
!  LOCALS
!      gamma_s_2d - delta_gamma*So (used in calculation of Fp_tau)
!----------------------------------------------------------------
      subroutine dstress_therm(F_tau,F_thrm_tau,F_thrm_tau_inv,det_thermal,  &
                               Fp_t,Fp_tau,Fe_tau,Fe_tau_inv, &
                               gst,gmt,s_tan,nslip,tpk_1d,delta_gamma,  &
                               stress,cauchy_2d,strain_1d,crp_strain_1d,dw_p)             
      use crystal_Titanium
      
      implicit none

            
      ! locals
      real(8) :: F_thermoelastic(3,3)        ! thermoelastic deformation gradient tensor
      real(8) :: F_thermoelastic_t(3,3)      ! its transpose
      real(8) :: F_thermoelastic_inv_t(3,3)  ! its inverse transpose
      real(8) :: F_thrm_tau_t(3,3)           ! transpose of the therm. grad. tensor
      real(8) :: F_thrm_tau_inv_t(3,3)       ! inv transpose of the therm. grad. tensor
      real(8) :: det_thermoelastic           ! det of thermoelastic deformation
      
      real(8) :: C_inv(3,3), detC
      real(8) :: lapack_work(27)
      integer :: INFO
      real(8), intent(in) :: gst(3,nslip_max),gmt(3,nslip_max),s_tan(3,3,nslip_max)
      integer, intent(in) :: nslip
      
      real(8), intent(in) :: F_tau(3,3),Fp_t(3,3)
      real(8), intent(out):: Fp_tau(3,3), Fe_tau(3,3), Fe_tau_inv(3,3)
      real(8), intent(in) :: F_thrm_tau(3,3)       ! thermal deformation gradient at tau
      real(8), intent(in) :: F_thrm_tau_inv(3,3)   ! inv of thermal deformation gradient at tau
      real(8), intent(in) :: det_thermal           ! determinant of the thermal deformation
      real(8), intent(in) :: tpk_1d(6)             ! 2nd PK at time t
      real(8), intent(in) :: delta_gamma(nslip_max)! plastic slip increments
      real(8), intent(out):: cauchy_2d(3,3)        ! cauchy at time tau - 2d
      real(8), intent(out):: stress(6)             ! cauchy at time tau - Voigt form
      real(8), intent(out):: strain_1d(6)          ! Total strain (Logarithmic)
      real(8), intent(out):: crp_strain_1d(6)      ! Creep strain (Green-Lagrange)
      real(8), intent(out):: dw_p                  ! plastic work density over this increment within this phase      

      real(8) :: strain_2d_green_lagrange(3,3),strain_2d_logarithmic(3,3)
      real(8) :: gamma_s_2d(3,3),gamma_s_fp(3,3),fe_tau_inv_t(3,3),det_fp_tau
      real(8) :: fp_tau_inv(3,3),Fe_tau_t(3,3),fe_tpk_2d_fe_t(3,3),det_fe_tau
      real(8) :: Fp_tau_TR(3,3),fp_t_fp(3,3),f_t_f(3,3),F_tau_TR(3,3)
      real(8) :: tpk_2d_g(3,3)
      real(8) :: eigvalue(3),eigvector(3,3),ckrone_2d(3,3)
           
      real(8) :: dmrot(3,nslip_max),dsrot(3,nslip_max)
      
      integer :: i,j,k,isys,nrot
     
      F_thermoelastic(:,:) =0.D0
      F_thrm_tau_t(:,:) = TRANSPOSE(F_thrm_tau)
      F_thrm_tau_inv_t(:,:) = TRANSPOSE(F_thrm_tau_inv)
      F_thermoelastic_inv_t(:,:)=0.D0
      
      ! calc the plastic deform generator (SIGMA delta_gamma * Schmidt)
      do i=1,3
        do j=1,3
          gamma_s_2d(i,j)=0.d0
          do isys=1,nslip
            gamma_s_2d(i,j)=gamma_s_2d(i,j)  &
               +delta_gamma(isys)*s_tan(i,j,isys)
          enddo
        enddo 
      enddo

      call mat33(gamma_s_fp,gamma_s_2d,Fp_t,3)


      ! Fp_tau, calculated from slips, to the first order. 
      ! deniz - FIX THIS - try to come up with a closed form
      !     expression for finite plastic deformation tensors,
      !     by projecting the (slip*schmidt) to a properly-chosen
      !     basis of the Lie algebra for the group of isochoric deformations
      !     and exploiting the exponential maps of these basis matrices.
      !     (See Glen-Mann matrices)
      do i=1,3
        do j=1,3
          Fp_tau(i,j)=Fp_t(i,j)+gamma_s_fp(i,j)
        enddo
      enddo
      
      call matinv3(Fp_tau,fp_tau_inv,det_fp_tau)
      if(det_fp_tau.eq.0.d0) then
        write(*,*) '0 divide'        
        stop
      endif

      ! 'NORMALIZE' Fp_tau to be volume-preserving (det=1)
      if(det_fp_tau.ne.1.d0) then
        do i=1,3
          do j=1,3
            Fp_tau(i,j)=(det_fp_tau**(-1.d0/3.d0))*Fp_tau(i,j)
          enddo
        enddo
        call matinv3(Fp_tau,fp_tau_inv,det_fp_tau)
      endif
     
      ! below, we will calculate Fe_tau,
      ! given F_tau, Fp_tau, F_therm_tau
      call sigmat(tpk_1d,tpk_2d_g)
      Fe_tau(:,:) = MATMUL(MATMUL(F_tau,fp_tau_inv),F_thrm_tau_inv)
      Fe_tau_t = TRANSPOSE(Fe_tau)
      call matinv3(Fe_tau,Fe_tau_inv,det_fe_tau)
      fe_tau_inv_t = TRANSPOSE(Fe_tau_inv)

      if(det_fe_tau.eq.0.d0) then
        write(*,*) '0 divide'
        stop
      endif
      
      ! form the thermoelastic deformation grad tensor.
      ! this transforms tangents from plastically deformed config. to current config.
      F_thermoelastic(:,:) = MATMUL(Fe_tau,F_thrm_tau)
      F_thermoelastic_t(:,:) = MATMUL(F_thrm_tau_t,Fe_tau_t)
      F_thermoelastic_inv_t(:,:) = MATMUL(fe_tau_inv_t,F_thrm_tau_inv_t)

      ! now we have Fe_tau, we will calculate
      ! the Cauchy stress (push-back operation)
      ! NOTE: 2nd PK is w.r.t thermally deformed configuration !
      det_thermoelastic = det_fe_tau * det_thermal
      cauchy_2d = (1.0/det_thermoelastic) *   &
                       MATMUL(MATMUL(F_thermoelastic,tpk_2d_g),  &
                                             F_thermoelastic_t)

      ! in Voigt notation
      ! NOTE THE ORDERING CONVENTION HERE !
      !        xx,yy,zz,xy,xz,yz
      stress(1)=cauchy_2d(1,1)
      stress(2)=cauchy_2d(2,2)
      stress(3)=cauchy_2d(3,3)
      stress(4)=0.5d0*(cauchy_2d(1,2)+cauchy_2d(2,1))
      stress(5)=0.5d0*(cauchy_2d(1,3)+cauchy_2d(3,1))
      stress(6)=0.5d0*(cauchy_2d(2,3)+cauchy_2d(3,2))

      ! calculate slip directions and plane normals in current config.
      ! use F_thermoelastic: maps from plastically deformed config to current config.
      ! dmrot=0.d0
      ! dsrot=0.d0
      ! do isys=1,nslip
      !   ! slip direction
      !   dsrot(:,isys)=MATMUL(F_thermoelastic(:,:),gst(:,isys))
      !   ! slip plane normal
      !   dmrot(:,isys)=1.d0/det_thermoelastic* &
      !                 MATMUL(F_thermoelastic_inv_t(:,:),gmt(:,isys))
      !enddo

!---------------creep strain-------------------------
!   Green-Lagrange strain for plastic deformation
      Fp_tau_TR = TRANSPOSE(Fp_tau)
      call mat33(fp_t_fp,Fp_tau_TR,Fp_tau,3)
      call cal_ckrone_2d(ckrone_2d)
      do i=1,3
         do j=1,3
            strain_2d_green_lagrange(i,j)=0.5d0*(fp_t_fp(i,j)-ckrone_2d(i,j))
         enddo
      enddo
      crp_strain_1d(1)=strain_2d_green_lagrange(1,1)
      crp_strain_1d(2)=strain_2d_green_lagrange(2,2)
      crp_strain_1d(3)=strain_2d_green_lagrange(3,3)
      crp_strain_1d(4)=strain_2d_green_lagrange(1,2)
      crp_strain_1d(5)=strain_2d_green_lagrange(1,3)
      crp_strain_1d(6)=strain_2d_green_lagrange(2,3)


!----- Kourosh ----------
!----- Strain  ----------
      call trans(F_tau,F_tau_TR)
      call mat33(f_t_f,F_tau_TR,F_tau,3)
      ! f_t_f = MATMUL(TRANSPOSE(F_tau),F_tau)
      
      CALL matinv3(f_t_f,C_inv,detC)

      !call jacobi(f_t_f,3,3,eigvalue,eigvector,nrot)
      eigvector = f_t_f
      call dsyev('V','U',3,eigvector,3,eigvalue,lapack_work,27,INFO)
      if (INFO.gt.0) then
         write(*,*) 'Eigenvalue calculation failed'
         stop
      endif
      
      do i=1,3
         do j=1,3
            strain_2d_logarithmic(i,j)=0.d0
            do k=1,3
               if (eigvalue(k).le.0.d0) then
                  write(*,*) 'Negative Eigenvalue @ Strain Cal.'
                  stop
               endif
         
               ! logarithmic strain
               strain_2d_logarithmic(i,j)=strain_2d_logarithmic(i,j) &
                        + 0.5d0*dlog(eigvalue(k))*eigvector(i,k)*eigvector(j,k)
            enddo
         enddo
      enddo
      
      call tr2to1(strain_2d_logarithmic,strain_1d) ! logarithmic strain is exported in 1D to state variables vector.      

      ! plastic work over this increment, within this crystal/phase
      dw_p = 0.0d0
      do i=1,3
        do j=1,3
           dw_p=dw_p+0.5d0*cauchy_2d(i,j)*  &
                  (gamma_s_2d(i,j)+gamma_s_2d(j,i))
        enddo
      enddo
      
      return
      end  



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!    MULTIPLICATION OF TWO 3*3 MATRICES

      subroutine mat33(x,y,z,l)
!       
       implicit double precision(a-h,o-z)
!       include 'aba_param.inc'
       dimension x(3,l),y(3,3),z(3,l)
!
        do i=1,3
          do j=1,l
             x(i,j)=0d0
             do  k=1,3
                x(i,j)=x(i,j)+y(i,k)*z(k,j)
             enddo
          enddo
        enddo
        return
        end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--               
      subroutine rotate(tch,r,t)                                      
!                                                                     
!  computes the matrix t=r*tch*(rtrans)  (all matrices are 3x3)
!                                                                     
       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension t(3,3),r(3,3),tch(3,3)                                
!                                                                     
      do 200 j=1,3                                                    
      do 200 i=1,3                                                    
      s=0.0                                                           
      do 100 l=1,3                                                    
        rjl=r(j,l)                                                    
        do 100 k=1,3                                                  
        s=s+r(i,k)*rjl*tch(k,l)                                       
 100  continue                                                      
      t(i,j)=s                                                        
 200  continue                                                        
      return
      end 


!------------------------------------------------------------------------

!    CALCULATES THE JACOBIAN FOR UMAT
!-------------------------------------------------------------------------
!      REFERING TO THESIS Pg 185
!      cl_4d - Lijkl
!      de_g  - Cijkl
!      d_tang_4d -Dijkl
!      g_tan_4d -Gmnkl
!      cj_tan_4d - Jijkl
!      b_tan_2d  - B
!      ck_tan_4d - Kijkl
!      q_tan     -Qijkl
!      r_tan     -Rij
!      Sijkl_tan - Sijkl
!      w_tan     - Wijkl 
!
! OUTPUT: ddsdde -- tangent elasticity matrix, dS/de
!
!  modified - deniz 
!      TEMP         - ambient temperature (in Kelvins)
!------------------------------------------------------------------------
      subroutine jacobian(F_t,F_tau,Fp_t,Fp_tau,s_tan,c_alpha_1d,  &
      tau,delta_gamma,dGamma_dTau,dtime,  &
      de_g,tpk_1d_tau,g_alpha_tau,  &
      ddsdde,ntens,nslip,tot_gamma,iCrystal,  &
      error)
      
      ! import UMAT options
      use options, only: ELASTIC_ONLY, PRISMATIC_ONLY, BASAL_ONLY, HARDENING_DISABLED, GND_DISABLED, &
                         USE_POWER_LAW, p_umat, q_umat, g0_add
      
      use crystal_Titanium
      
      implicit none
      
      real(8), intent(in) :: F_t(3,3),F_tau(3,3),Fp_t(3,3),Fp_tau(3,3)
      real(8), intent(in) :: s_tan(3,3,nslip_max)
      real(8), intent(in) :: c_alpha_1d(6,nslip_max)
      real(8), intent(in) :: tau(nslip_max)
      real(8), intent(in) :: delta_gamma(nslip_max)
      real(8), intent(in) :: dGamma_dTau(nslip_max)
      real(8), intent(in) :: dtime
      real(8), intent(in) :: de_g(3,3,3,3)            ! rotated elasticity tensor 
      real(8), intent(in) :: tpk_1d_tau(6)  
      real(8), intent(in) :: g_alpha_tau(nslip_max)
      real(8), intent(out):: ddsdde(ntens,ntens)
      integer, intent(in) :: nslip, ntens
      real(8), intent(in) :: tot_gamma(nslip_max)
      integer, intent(in) :: iCrystal
      integer, intent(out):: error
      integer :: INFO
      
      real(8) :: Fe_t_t(3,3),Fe_t(3,3),Fe_tau_inv(3,3),Fp_t_inv(3,3),  &
           cl1(3,3),cl2(3,3),fe_gamma_s(3,3), &
           cl_4d(3,3,3,3), &
           d_tang_4d(3,3,3,3),g_tan_4d(3,3,3,3,nslip_max),tpk_2d_g(3,3),  &
           cj_tan_4d(3,3,3,3,nslip_max),cj1_tan_4d(3,3,3,3,nslip_max),   &
           b_tan_2d(3,3,nslip_max),ck_tan_4d(3,3,3,3),ck_tan_2d(6,6),  &
           ckrone_symm(3,3,3,3),cb_tan(3,3,3,3),c_alpha(3,3,nslip_max),  &
           gamma_j(3,3,3,3),q_tan(3,3,3,3),r_tan(3,3,nslip_max),  &
      ck_tan_2d_inv(6,6),ck_tan_2d_inv_temp(6,6),ck_tan_4d_inv(3,3,3,3),  &
           q1_tan(3,3,3,3),  &
           sijkl_tan_1st(3,3,3,3),sijkl_tan_2nd(3,3,3,3),r_s(3,3,3,3),  &
         t_fe_t(3,3),w_tan_1st(3,3,3,3),fe_q(3,3,3,3),F_t_inv(3,3),  &
           w_tan_2nd(3,3,3,3),w_tan(3,3,3,3),w_tan_3rd(3,3,3,3),  &
           w_tan_4th(3,3,3,3),fe_stress_2d(3,3),fe_stress_2d_fe_t(3,3),  &
         sijkl_tan_3rd(3,3,3,3),sijkl_tan(3,3,3,3),  &
           sijkl_tan_fe_inv(3,3),rt_u_fe_t(3,3),  &
           w_tan_2d(6,6),det,det_fe_tau

      real(8) :: q2_tan(3,3,3,3) 
      real(8) :: Fe_tau(3,3),Fe_tau_t(3,3),  &
                 fp_tau_inv(3,3),gamma_s_2d(3,3), &
                 dfgrdt(3,3),rt_tau(3,3),ut_tau(3,3)
      integer :: i,j,k,l,m,n,ip,iq,isys
      
      real(8), parameter :: kB=1.3806488d-23   ! boltzmann constant (J/K)

      error = 0

      call matinv3(F_t,F_t_inv,det)
      
      call mat33(dfgrdt,F_tau,F_t_inv,3)

      call rudcmp(dfgrdt,rt_tau,ut_tau)

      call matinv3(Fp_t,fp_t_inv,det)

      call mat33(Fe_t,F_t,fp_t_inv,3)
      call mat33(cl1,ut_tau,Fe_t,3)
      call trans(Fe_t,fe_t_t)
      call mat33(cl2,fe_t_t,ut_tau,3)
      call sigmat(tpk_1d_tau,tpk_2d_g)


      call matinv3(Fp_tau,fp_tau_inv,det)

      call mat33(Fe_tau,F_tau,fp_tau_inv,3)
      call trans(Fe_tau,Fe_tau_t)
      call matinv3(Fe_tau,Fe_tau_inv,det_fe_tau)

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              cl_4d(i,j,k,l)=fe_t_t(i,k)*cl1(l,j)+cl2(i,k)*Fe_t(l,j)
            enddo
          enddo
        enddo
      enddo

      call mat44(d_tang_4d,de_g,cl_4d) 
      d_tang_4d = 0.5d0*d_tang_4d

      do isys=1,nslip
        do m=1,3
          do n=1,3
            do k=1,3
              do l=1,3
                g_tan_4d(m,n,k,l,isys)=0d0
                do ip=1,3
                  g_tan_4d(m,n,k,l,isys)=g_tan_4d(m,n,k,l,isys)  &
          +cl_4d(m,ip,k,l)*s_tan(ip,n,isys)  &
          +s_tan(ip,m,isys)*cl_4d(ip,n,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo 
      enddo



      call mat44crys(cj1_tan_4d,de_g,g_tan_4d,nslip)
      do isys=1,nslip
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                cj_tan_4d(i,j,k,l,isys)=0.5d0  &
       *cj1_tan_4d(i,j,k,l,isys)
              enddo
            enddo
          enddo
        enddo
      enddo
 

      do isys=1,nslip
         do i=1,3
            do j=1,3
               b_tan_2d(i,j,isys)=0.5d0*dGamma_dTau(isys)  &
               *(s_tan(i,j,isys)+s_tan(j,i,isys))
            enddo
         enddo
      enddo

      call sigmat_crys(c_alpha_1d,c_alpha,nslip)  
      call mat33crys(cb_tan,c_alpha,b_tan_2d,nslip)
      call cal_ckrone_symm(ckrone_symm)  ! this has to be the symm 4th order tensor. see derivation.
!
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              ck_tan_4d(i,j,k,l)=ckrone_symm(i,j,k,l)+cb_tan(i,j,k,l)
            enddo
          enddo
        enddo
      enddo          

!--- (b) qijkl=
      call tr4to2(ck_tan_4d,ck_tan_2d)

      ! use Lapack:
      ck_tan_2d_inv(:,:) = ck_tan_2d(:,:)
      call lapack_invert(6,ck_tan_2d_inv,INFO)
      if (INFO.NE.0) then
         ! matrix inversion error
         error = -1
         return
      endif

      call tr2to4(ck_tan_2d_inv,ck_tan_4d_inv)

      do m=1,3
        do n=1,3
          do k=1,3
            do l=1,3
              gamma_j(m,n,k,l)=0d0
              do isys=1,nslip 
                gamma_j(m,n,k,l)=gamma_j(m,n,k,l)  &
            +delta_gamma(isys)*cj_tan_4d(m,n,k,l,isys)
              enddo
            enddo
          enddo
        enddo
      enddo

      call mat44(q1_tan,ck_tan_4d_inv,d_tang_4d)
      call mat44(q2_tan,ck_tan_4d_inv,gamma_j)
      
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              q_tan(i,j,k,l)=q1_tan(i,j,k,l)-q2_tan(i,j,k,l)
            enddo
          enddo
        enddo
      enddo


      do isys=1,nslip
        do i=1,3
          do j=1,3
            r_tan(i,j,isys)=0d0
            do k=1,3 
              do l=1,3
                r_tan(i,j,isys)=r_tan(i,j,isys)  &
          +b_tan_2d(k,l,isys)*q_tan(k,l,i,j)      
              enddo
            enddo
          enddo
        enddo
      enddo 


      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_1st(i,j,k,l)=rt_tau(i,k)*Fe_t(l,j)
            enddo
          enddo
        enddo
      enddo
      
      
      do i=1,3
         do j=1,3
            gamma_s_2d(i,j)=0.d0
            do isys=1,nslip
               
               gamma_s_2d(i,j)=gamma_s_2d(i,j)+delta_gamma(isys)  &
                    *s_tan(i,j,isys)
            enddo
         enddo
      enddo

      call mat33(fe_gamma_s,Fe_t,gamma_s_2d,3)


      do i=1,3 
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_2nd(i,j,k,l)=rt_tau(i,k)*fe_gamma_s(l,j)
            enddo
          enddo
        enddo
      enddo

      call mat33crys(r_s,r_tan,s_tan,nslip)
      call mat33(rt_u_fe_t,dfgrdt,Fe_t,3)

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan_3rd(i,j,k,l)=0d0
              do ip=1,3
                sijkl_tan_3rd(i,j,k,l)=sijkl_tan_3rd(i,j,k,l)  &
          +rt_u_fe_t(i,ip)*r_s(k,l,ip,j) 
              enddo
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              sijkl_tan(i,j,k,l)=sijkl_tan_1st(i,j,k,l)  &
        -sijkl_tan_2nd(i,j,k,l)-sijkl_tan_3rd(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

      call mat33(t_fe_t,tpk_2d_g,Fe_tau_t,3) 
      call mat42(w_tan_1st,sijkl_tan,t_fe_t)


      do i=1,3
        do n=1,3
          do k=1,3
            do l=1,3
              fe_q(i,n,k,l)=0d0
              do m=1,3
                fe_q(i,n,k,l)=fe_q(i,n,k,l)+Fe_tau(i,m)  &
          *q_tan(m,n,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo

      call mat42(w_tan_2nd,fe_q,Fe_tau_t)


      call mat33(fe_stress_2d,Fe_tau,tpk_2d_g,3)

      

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan_3rd(i,j,k,l)=0d0
              do n=1,3
                w_tan_3rd(i,j,k,l)=w_tan_3rd(i,j,k,l)  &
        +fe_stress_2d(i,n)*sijkl_tan(j,n,k,l)
              enddo
            enddo
          enddo
        enddo
      enddo


      call mat33(fe_stress_2d_fe_t,fe_stress_2d,Fe_tau_t,3)
      do k=1,3
        do l=1,3
          sijkl_tan_fe_inv(k,l)=0d0 
          do ip=1,3
            do iq=1,3
              sijkl_tan_fe_inv(k,l)=sijkl_tan_fe_inv(k,l)  &
        +sijkl_tan(ip,iq,k,l)*Fe_tau_inv(iq,ip)   
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan_4th(i,j,k,l)=fe_stress_2d_fe_t(i,j)  &
        *sijkl_tan_fe_inv(k,l)
            enddo
          enddo
        enddo
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              w_tan(i,j,k,l)=1.0/det_fe_tau*(w_tan_1st(i,j,k,l)  &
        +w_tan_2nd(i,j,k,l)  &
        +w_tan_3rd(i,j,k,l)-w_tan_4th(i,j,k,l))
            enddo
          enddo
        enddo
      enddo

      call tr4to2(w_tan,w_tan_2d)

       do i=1,6
       do j=1,3
          ddsdde(i,j)=w_tan_2d(i,j)
       enddo
         do j=4,6
            ddsdde(i,j)=w_tan_2d(i,j)/2
         enddo
      enddo



    
      return
      end

      subroutine mat44crys(x,y,z,nslip)

      implicit double precision(a-h,o-z)
      include 'slip_sys_param.inc'    ! nslipfam, nslip_max, HCP/BCC, pri_A, trans_B
       
      dimension x(3,3,3,3,nslip_max),y(3,3,3,3),z(3,3,3,3,nslip_max) 

      do isys=1,nslip
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                x(i,j,k,l,isys)=0d0
                do m=1,3
                  do n=1,3
                    x(i,j,k,l,isys)=x(i,j,k,l,isys)  &
                      +y(i,j,m,n)*z(m,n,k,l,isys)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      return
      end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!    CONVERTS FOURTH ORDER TENSOR TO EQUIVALENT SECOND ORDER TENSOR

      subroutine tr4to2(a_4d,a_2d)
      implicit real*8(a-h,o-z)
!     include 'aba_param.inc'
      dimension a_4d(3,3,3,3),a_2d(6,6)

      call rclear66(a_2d)

      do i=1,3
        a_2d(i,1)=a_4d(i,i,1,1)
        a_2d(i,2)=a_4d(i,i,2,2)
        a_2d(i,3)=a_4d(i,i,3,3)
        a_2d(i,4)=a_4d(i,i,1,2)+a_4d(i,i,2,1)
        a_2d(i,5)=a_4d(i,i,1,3)+a_4d(i,i,3,1)
        a_2d(i,6)=a_4d(i,i,2,3)+a_4d(i,i,3,2)
      enddo

      do i=1,3
        a_2d(4,i)=a_4d(1,2,i,i)
        a_2d(5,i)=a_4d(1,3,i,i)
        a_2d(6,i)=a_4d(2,3,i,i)
      enddo

      a_2d(4,4)=a_4d(1,2,1,2)+a_4d(1,2,2,1)
      a_2d(4,5)=a_4d(1,2,1,3)+a_4d(1,2,3,1)
      a_2d(4,6)=a_4d(1,2,2,3)+a_4d(1,2,3,2)

      a_2d(5,4)=a_4d(1,3,1,2)+a_4d(1,3,2,1)
      a_2d(5,5)=a_4d(1,3,1,3)+a_4d(1,3,3,1)
      a_2d(5,6)=a_4d(1,3,2,3)+a_4d(1,3,3,2)

      a_2d(6,4)=a_4d(2,3,1,2)+a_4d(2,3,2,1)
      a_2d(6,5)=a_4d(2,3,1,3)+a_4d(2,3,3,1)
      a_2d(6,6)=a_4d(2,3,2,3)+a_4d(2,3,3,2)

      return
      end


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!     CONVERTS A SECOND ORDER TENSOR TO A EQUIVALENT FOURTH ORDER TENSOR

      subroutine tr2to4(a_2d,a_4d)
      implicit real*8(a-h,o-z)

      dimension a_4d(3,3,3,3),a_2d(6,6)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  a_4d(i,j,k,l)=0d0
               enddo
            enddo
         enddo
      enddo
!
!---1-3
      do i=1,3
         a_4d(i,i,1,1)=a_2d(i,1)
         a_4d(i,i,2,2)=a_2d(i,2)
         a_4d(i,i,3,3)=a_2d(i,3)
         a_4d(i,i,1,2)=a_2d(i,4)
         a_4d(i,i,1,3)=a_2d(i,5)
         a_4d(i,i,2,3)=a_2d(i,6)
      enddo

      a_4d(1,2,1,1)=a_2d(4,1)
      a_4d(1,2,2,2)=a_2d(4,2)
      a_4d(1,2,3,3)=a_2d(4,3)
      a_4d(1,2,1,2)=a_2d(4,4)
      a_4d(1,2,1,3)=a_2d(4,5)
      a_4d(1,2,2,3)=a_2d(4,6)

      a_4d(1,3,1,1)=a_2d(5,1)
      a_4d(1,3,2,2)=a_2d(5,2)
      a_4d(1,3,3,3)=a_2d(5,3)
      a_4d(1,3,1,2)=a_2d(5,4)
      a_4d(1,3,1,3)=a_2d(5,5)
      a_4d(1,3,2,3)=a_2d(5,6)

      a_4d(2,3,1,1)=a_2d(6,1)
      a_4d(2,3,2,2)=a_2d(6,2)
      a_4d(2,3,3,3)=a_2d(6,3)
      a_4d(2,3,1,2)=a_2d(6,4)
      a_4d(2,3,1,3)=a_2d(6,5)
      a_4d(2,3,2,3)=a_2d(6,6)

      a_4d(2,1,1,1)=a_2d(4,1)
      a_4d(2,1,2,2)=a_2d(4,2)
      a_4d(2,1,3,3)=a_2d(4,3)
      a_4d(2,1,1,2)=a_2d(4,4)
      a_4d(2,1,1,3)=a_2d(4,5)
      a_4d(2,1,2,3)=a_2d(4,6)

      a_4d(3,1,1,1)=a_2d(5,1)
      a_4d(3,1,2,2)=a_2d(5,2)
      a_4d(3,1,3,3)=a_2d(5,3)
      a_4d(3,1,1,2)=a_2d(5,4)
      a_4d(3,1,1,3)=a_2d(5,5)
      a_4d(3,1,2,3)=a_2d(5,6)

      a_4d(3,2,1,1)=a_2d(6,1)
      a_4d(3,2,2,2)=a_2d(6,2)
      a_4d(3,2,3,3)=a_2d(6,3)
      a_4d(3,2,1,2)=a_2d(6,4)
      a_4d(3,2,1,3)=a_2d(6,5)
      a_4d(3,2,2,3)=a_2d(6,6)
      
      return
      end

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine mat44(x,y,z)
!
       implicit double precision(a-h,o-z)
!      include 'aba_param.inc'
      dimension x(3,3,3,3),y(3,3,3,3),z(3,3,3,3) 
!
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
              do m=1,3
                do n=1,3
                  x(i,j,k,l)=x(i,j,k,l)+y(i,j,m,n)*z(m,n,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      return
      end           
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!    CALCULATES THE TRANSPOSE OF A MATRIX

       subroutine trans(x,x_t)
        implicit double precision(a-h,o-z)
!       include 'aba_param.inc'
       dimension x(3,3),x_t(3,3) 
!
       do i=1,3
         do j=1,3
           x_t(j,i)=x(i,j)
         enddo
       enddo
       return
       end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine mat33crys(x,y,z,nslip)

      implicit double precision (a-h,o-z)
      include 'slip_sys_param.inc'    ! nslipfam, nslip_max, HCP/BCC, pri_A, trans_B

      dimension x(3,3,3,3),y(3,3,nslip),z(3,3,nslip_max)

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
              do isys=1,nslip
                x(i,j,k,l)=x(i,j,k,l)  &
          +y(i,j,isys)*z(k,l,isys)
              enddo
            enddo
          enddo
        enddo
      enddo
!
      return
      end





      subroutine mat42(x,y,z)
       implicit real*8 (a-h,o-z)
!      include 'aba_param.inc'
      dimension x(3,3,3,3),y(3,3,3,3),z(3,3)
!
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              x(i,j,k,l)=0d0
                do n=1,3
                  x(i,j,k,l)=x(i,j,k,l)  &
                +y(i,n,k,l)*z(n,j)
                enddo
              enddo
            enddo
          enddo
        enddo
!
        return
        end


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-
!     FORMS A FOURTH ORDER IDENTITY TENSOR

      subroutine cal_ckrone(ckrone)

       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension ckrone(3,3,3,3)

      ckrone=0.d0
     
      ckrone(1,1,1,1)=1.d0
      ckrone(2,2,2,2)=1.d0
      ckrone(3,3,3,3)=1.d0
      ckrone(1,2,1,2)=1.d0
      ckrone(1,3,1,3)=1.d0
      ckrone(2,3,2,3)=1.d0
      ckrone(2,1,2,1)=1.d0
      ckrone(3,1,3,1)=1.d0
      ckrone(3,2,3,2)=1.d0
      
      return
      end
      
! forms a symmetric 4th order unit tensor
      subroutine cal_ckrone_symm(ckrone_symm)

       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension ckrone_symm(3,3,3,3)

      ckrone_symm = 0.D0
      
      ckrone_symm(1,1,1,1)=1.d0
      ckrone_symm(2,2,2,2)=1.d0
      ckrone_symm(3,3,3,3)=1.d0
      
      ckrone_symm(1,2,1,2)=0.5d0
      ckrone_symm(1,2,2,1)=0.5d0
      ckrone_symm(2,1,1,2)=0.5d0
      ckrone_symm(2,1,2,1)=0.5d0
      
      ckrone_symm(1,3,1,3)=0.5d0
      ckrone_symm(1,3,3,1)=0.5d0
      ckrone_symm(3,1,1,3)=0.5d0
      ckrone_symm(3,1,3,1)=0.5d0

      ckrone_symm(2,3,2,3)=0.5d0
      ckrone_symm(2,3,3,2)=0.5d0
      ckrone_symm(3,2,2,3)=0.5d0
      ckrone_symm(3,2,3,2)=0.5d0
      
      return
      end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!      FORMS A 1-D IDENTITY VECTOR

      subroutine cal_ckrone_1d(ckrone_1d)
       implicit double precision (a-h,o-z)
      dimension ckrone_1d(6)

      ckrone_1d = 0.d0
      ckrone_1d(1) = 1.d0
      ckrone_1d(2) = 1.d0
      ckrone_1d(3) = 1.d0

      return
      end
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!    FORMS A SECOND ORDER IDENTITY TENSOR

      subroutine cal_ckrone_2d(ckrone_2d)

       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
       dimension ckrone_2d(3,3)
!
       do i=1,3
          do j=1,3
             if(i.eq.j) then
                ckrone_2d(i,j)=1.0d0
             else
                ckrone_2d(i,j)=0.0d0
             endif
          enddo
       enddo
!     
       return
       end
!---  +----1----+----2----+----3----+----4----+----5----+----6----+----7--
!     SETS ALL THE COMPONENTS OF A 6*6 MATRIX TO ZERO

      subroutine rclear66(x)

       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension x(6,6)
!
      do i=1,6
        do j=1,6
          x(i,j)=0.d0
        enddo
      enddo
!
      return
      end
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      subroutine sigmat_crys(x,xx,nslip)    
!                                                                     
!  arranges into a symmetric array xx, the six commponents of a vector x
!                                                                     
       implicit double precision (a-h,o-z)
!      include 'aba_param.inc'
      dimension x(6,nslip),xx(3,3,nslip)
!
      do isys=1,nslip
        xx(1,1,isys)=x(1,isys)
        xx(2,2,isys)=x(2,isys)
        xx(3,3,isys)=x(3,isys)
        xx(2,3,isys)=x(6,isys)
        xx(3,2,isys)=x(6,isys)
        xx(3,1,isys)=x(5,isys)
        xx(1,3,isys)=x(5,isys)
        xx(1,2,isys)=x(4,isys)
        xx(2,1,isys)=x(4,isys)
      enddo

      return
      end


      subroutine lapack_invert(N,A,INFO)
      implicit none
      integer, intent(in) :: N
      real(8), intent(inout):: A(N,N)
      integer, intent(out) :: INFO
      real(8) :: IPIV(N), WORK(N*N)
      
      CALL DGETRF( N, N, A, N, IPIV, INFO )
      if(INFO.NE.0) return
      CALL DGETRI( N, A, N, IPIV, WORK, N*N, INFO)
      end subroutine
!******************************************************
       subroutine inverseMat(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0


      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine
!----------------------------------------------------
!     CALCULATES THE SIGN 
             
      double precision function dsignf(x)
      double precision x
      if(x.ge.0)then
         dsignf=+1.d0
      else
         dsignf=-1.d0
      endif
      return
      end
      
!-----------------------------------------------------
      subroutine rotate_crys_vector(a,q,a_prime,nslip)
      implicit real*8(a-h,o-z)
      include 'slip_sys_param.inc'      
      dimension a(3,nslip_max),a_prime(3,nslip_max),q(3,3)
      do isys=1,nslip
         do i=1,3
            a_prime(i,isys)=0.d0
            do j=1,3
               a_prime(i,isys)=a_prime(i,isys)+q(i,j)*a(j,isys)
            enddo
         enddo
      enddo
      return
      end


!------------------------------------------------------------
!   TRANSFORMS THE FOURTH ORDER ELASTICITY TENSOR FROM 
!   CRYSTAL SYSTEM TO GLOBAL 

      subroutine rotate4d(ddsdde_4d,tr,de_gg)
      implicit real*8(a-h,o-z)
!     transforms a fourth order tensor
      dimension ddsdde_4d(3,3,3,3),de_gg(3,3,3,3),tr(3,3)
      dimension de_g11(3,3,3,3),de_g22(3,3,3,3),de_g33(3,3,3,3)

      do i=1,3
         do n=1,3
            do io=1,3
               do ip=1,3
                  de_g11(i,n,io,ip)=0d0
                  do m=1,3
                     de_g11(i,n,io,ip)=de_g11(i,n,io,ip)  &
                          +tr(i,m)*ddsdde_4d(m,n,io,ip)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      do i=1,3
         do j=1,3
            do io=1,3
               do ip=1,3
              de_g22(i,j,io,ip)=0d0
              do n=1,3
                 de_g22(i,j,io,ip)=de_g22(i,j,io,ip)  &
                      +tr(j,n)*de_g11(i,n,io,ip)
              enddo
           enddo
          enddo
       enddo
      enddo
      
      do i=1,3
         do j=1,3
            do k=1,3
               do ip=1,3
                  de_g33(i,j,k,ip)=0d0
                  do io=1,3
                     de_g33(i,j,k,ip)=de_g33(i,j,k,ip)  &
                          +tr(k,io)*de_g22(i,j,io,ip)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      do i=1,3
         do j=1,3
          do k=1,3
             do l=1,3
              de_gg(i,j,k,l)=0d0
              do ip=1,3
                 de_gg(i,j,k,l)=de_gg(i,j,k,l)  &
                      +tr(l,ip)*de_g33(i,j,k,ip)
              enddo
           enddo
        enddo
      enddo
      enddo
      
      return
      end


      subroutine tr2to1(a_2d,a_1d)
      implicit real*8(a-h,o-z)
      dimension a_2d(3,3),a_1d(6)
      a_1d(1)=a_2d(1,1)
      a_1d(2)=a_2d(2,2)
      a_1d(3)=a_2d(3,3)
      a_1d(4)=0.5d0*(a_2d(1,2)+a_2d(2,1))
      a_1d(5)=0.5d0*(a_2d(1,3)+a_2d(3,1))
      a_1d(6)=0.5d0*(a_2d(2,3)+a_2d(3,2))
      return
      end

      subroutine tr2to1crys(a_2d,a_1d,nslip)
      implicit real*8(a-h,o-z)
      include 'slip_sys_param.inc' ! nslip_max
      dimension a_2d(3,3,nslip_max),a_1d(6,nslip_max)
      do isys=1,nslip
         a_1d(1,isys)=a_2d(1,1,isys)
         a_1d(2,isys)=a_2d(2,2,isys)
         a_1d(3,isys)=a_2d(3,3,isys)
         a_1d(4,isys)=0.5d0*(a_2d(1,2,isys)+a_2d(2,1,isys))
         a_1d(5,isys)=0.5d0*(a_2d(1,3,isys)+a_2d(3,1,isys))
         a_1d(6,isys)=0.5d0*(a_2d(2,3,isys)+a_2d(3,2,isys))
      enddo
      return
      end
      
!----------------------------------------------------------------
!       CALCULATES THE RESIDUE FOR NEWTON RAPHSON FOR T*
!------------------------------------------------------
!       g_n             (out) residual = trial stress - stress calculated via slips: (t_tr_1d + slips*c_alpha_1d)
!       xjn_2d          (out) Jacobian for Newton Raphson
!       tpk_1d_tau      (in)  the trial STRESS at tau.
!       chi_tau         (in)  (in MAT.MODULE) back stress
!       t_tr_1d         (in)  zeroth order approximation to the stress
!       s_tan           (in)  Schmidt tensor in plastically deformed config.
!       TEMP            (in)  ambient temp
!       tau(slip)       (out) resolved shear stress due to trial stress
!       delta_gamma(slip)(out)plastic slip increments at each slip system, due to proposed/trial STRESS t_tr_1d
!----------------------------------------------------------------
      
      subroutine stress_residue(TEMP,dtime,iCrystal,s_tan,nslip, &
                                g_alpha_tau,       &
                                t_tr_1d,tpk_1d_tau,tau,s_gg_2d,  &
                                g_n,xjn_2d,iflag,c_alpha_1d,     &
                                tot_gamma,delta_gamma,dGamma_dTau, &
                                iatg)
      
      use material_Ti6242
      
      ! import UMAT options
      use options, only: ELASTIC_ONLY, PRISMATIC_ONLY, BASAL_ONLY, HARDENING_DISABLED, GND_DISABLED, &
                         DEBUG_RESID, &
                         USE_POWER_LAW, p_umat, q_umat, g0_add

      implicit none
      
      real(8), intent(in) :: TEMP,dtime              ! ambient temperature (*K)
      integer, intent(in) :: iCrystal                ! 1=HCP, 2=BCC
      real(8), intent(in) :: s_tan(3,3,nslip_max)
      integer, intent(in) :: nslip
      real(8), intent(in) :: g_alpha_tau(nslip_max)  ! hardening parameters at time tau (generic/crystal agnostic)
      real(8), intent(in) :: t_tr_1d(6)              ! zeroth order (trial) stress for time tau
      real(8), intent(in) :: tpk_1d_tau(6)
      real(8), intent(out):: tau(nslip_max)
      real(8), intent(in) :: s_gg_2d(6,6)
      real(8), intent(in) :: c_alpha_1d(6,nslip_max)
      real(8), intent(in) :: tot_gamma(nslip_max)
      real(8), intent(out):: delta_gamma(nslip_max)
      real(8), intent(out):: dGamma_dTau(nslip_max)  ! derivative of the flow rule wrt tau, for jacobian/tang. stiffness calculation

      integer, intent(in) :: iflag
      integer, intent(out):: iatg
      
      real(8), intent(out) :: g_n(6)
      
      real(8) :: tpk_2d(3,3),gamma_c_alpha(6)
      real(8) :: b_alpha_der(3,3,nslip_max), &
                 xjn(3,3,3,3),xjn_2d(6,6),c_alpha_2dd(3,3,nslip_max)
      real(8) :: e_strain(6,1),e_strain_2d(3,3)
      real(8) :: fe_t_fe_tpk(3,3),fe_t_fe(3,3)
      integer :: isys,i,j,k,l
      real(8) :: g_inv,rate_inv,sc,atg,atg_m
      
      real(8), parameter :: kB=1.3806488d-23  ! boltzmann constant (J/K)
      
      iatg = 0 !reset error flag:
      ! iatg:
      !  0: success
      !  1: tau/g > 2
      !  2: g is negative
      
!----------------------------------------------------------
!
!       for feTfe
!
      CALL matmul_flex(s_gg_2d,tpk_1d_tau,e_strain,6,6,1)
      call sigmat(e_strain,e_strain_2d)
      do i=1,3
         do j=1,3
            ! deniz - FIX THIS
            ! here, we are calculating fe_t_fe to give us the resolved shear stress wrt thermally deformed config.
            ! should this be f_thermoelast_t_f_thermoelast so res shear str is wrt plastically deformed config?
            if(i.eq.j)then
               fe_t_fe(i,j)=1.d0+(2.d0*e_strain_2d(i,j))
            else
               fe_t_fe(i,j)=2.d0*e_strain_2d(i,j)
            endif
         enddo
      enddo        
      call sigmat(tpk_1d_tau,tpk_2d)
      call mat33(fe_t_fe_tpk,fe_t_fe,tpk_2d,3)
               
      ! deniz -- comment
      ! here:
      !  tpk_1d_tau  (in)  the trial STRESS at tau.
      !  t_tr_1d     (in)  zeroth order approximation to the stress
      !  g_n         (out) residual = trial stress - stress calculated via slips: (t_tr_1d + slips*c_alpha_1d)
!
!-------------------------------------------------------------

      do isys=1,nslip
         tau(isys)=0.0d0
         do i=1,3
          do j=1,3
          
            ! deniz - FIX THIS -- CAUTION
            ! here, we are calculating resolved shear stress wrt thermally deformed config, not plastically deformed config.
            ! however, s_tan is wrt plastically deformed config
             tau(isys)=tau(isys)+fe_t_fe_tpk(i,j)*s_tan(i,j,isys) 
       ! here: 
       !    s_tan: schmidt tensor for slip system isys, wrt undeformed (or equiv, plastically deformed) config
       !    fe_t_fe_tpk: (C_e)(S*), mandell stress tensor at intermediate config.
          enddo
         enddo
           

         tau(isys)=tau(isys)-chi_tau(isys)
           
           
      enddo
      
      do isys=1,nslip
       
       
         sc=tau(isys)/(g_alpha_tau(isys))
         atg=dabs(sc)
         
         ! deniz -- elastic switch
         if(ELASTIC_ONLY) atg = 0.D0

         ! DISABLE ALL SLIPS BUT PRISMATIC
         if(PRISMATIC_ONLY.and.iCrystal==HCP.and. &
           (isys<=3.or.isys>=7)) then
            atg = 0.d0
         endif
         ! DISABLE ALL SLIPS BUT BASAL
         if(BASAL_ONLY.and.iCrystal==HCP.and. &
            isys>3) then
            atg = 0.d0
         endif         

         ! check if power law would blow up
         if(atg.gt.2.d0) then               
            iatg=1   ! tau/g > 2
            return
         endif    
         
         ! check if the shear res is physical.
         if(g_alpha_tau(isys).LT.0.0)then
            iatg=2
            return
         endif
                     
         if(USE_POWER_LAW) then
            ! POWER LAW
            delta_gamma(isys)=a_dot(isys)*dtime*  &
                              atg**(1.d0/rate_exp(isys))*  &
                              dsign(1.d0,tau(isys))
         else
            delta_gamma(isys)=a_dot_therm(isys)*dtime*  &
                        exp(-E_act(isys)/(kB*TEMP)*(1-atg**p_umat)**q_umat)*  &
                        dsign(1.d0,tau(isys))
            !write(*,*) 'using therm law!!'
         endif   
      enddo
      
      do i=1,6
         gamma_c_alpha(i)=DOT_PRODUCT(delta_gamma(1:nslip),c_alpha_1d(i,1:nslip))
      enddo
      
      g_n(:)=tpk_1d_tau(:)-t_tr_1d(:)+gamma_c_alpha(:)
      
      do i=1,6
         if ( isNaN(g_n(i)) ) then   ! NaN
            iatg=1
            return
         endif
      enddo
      ! here
      !     t_tr_1d: C:E_e0_tau where C is elastic stiffness tensor, and E_e0_tau = 0.5*(C_e0_tau-I), where C_e0_tau = Fp_t^-T F_tau^T F_tau Fp_t^-1 , zeroth order approx to elastic deform tensor
         
      ! Calculating the jacobian for nr iteration
      ! ----------------------------------------------
      if(iflag.eq.1)then
         do isys=1,nslip

         !  g_alpha_tau is an input. if negative, problem is somewhere else 
            ! check if the shear res is physical.
            if(g_alpha_tau(isys).LT.0.0)then
               iatg=2
               return
            endif
            
            g_inv=1.d0/dabs(g_alpha_tau(isys))
            rate_inv=1.d0/rate_exp(isys)
            sc=tau(isys)/g_alpha_tau(isys)
            atg=dabs(sc)
            

            ! deniz -- elastic switch
            if(ELASTIC_ONLY) atg = 0.D0

            ! DISABLE ALL SLIPS BUT PRISMATIC
            if(PRISMATIC_ONLY.and.iCrystal==HCP.and. &
               (isys<=3.or.isys>=7)) then
               atg = 0.d0
            endif
            ! DISABLE ALL SLIPS BUT BASAL
            if(BASAL_ONLY.and.iCrystal==HCP.and. &
               isys>3) then
               atg = 0.d0
            endif         
            
            if(atg.gt.2.d0) then
               iatg=1
               return
            endif       
            
            ! derivative of flow rule -- for Jacobian calculation
            if(USE_POWER_LAW) then
               ! POWER LAW
               dGamma_dTau(isys)=a_dot(isys)*dtime*rate_inv*atg**(rate_inv-1.d0)  &
                        *g_inv
            else
               ! Thermally activated flow
               atg_m = atg
               if (atg.le.0.05.and.q_umat.gt.1.0D0) atg_m = 0.05
               dGamma_dTau(isys)=a_dot_therm(isys)*dtime*  &
                  exp(-E_act(isys)/(kB*TEMP)*((1-atg_m**p_umat)**q_umat))*  &
                  E_act(isys)/(kB*TEMP)*q_umat*p_umat*((1-atg_m**p_umat)**(q_umat-1))*(atg_m**(p_umat-1))*g_inv
               !write(*,*) 'using therm law!!'
            end if

            
            do i=1,3
             do j=1,3
                b_alpha_der(i,j,isys)=0.5d0*dGamma_dTau(isys)  &
                       *(s_tan(i,j,isys)+s_tan(j,i,isys))   ! here, assuming C_elastic equal IDENTITY
             enddo
            enddo
         enddo
         
         
         call sigmat_crys(c_alpha_1d,c_alpha_2dd,nslip)
         
         do i=1,3
            do j=1,3
             do k=1,3
                do l=1,3
                   xjn(i,j,k,l)=0.0d0
                   do isys=1,nslip
                    xjn(i,j,k,l)=xjn(i,j,k,l)  &
                          +c_alpha_2dd(i,j,isys)*b_alpha_der(k,l,isys)
                   enddo
                   
                enddo
             enddo
            enddo
         enddo
         
         
         
         call tr4to2(xjn,xjn_2d)      
         
         do i=1,6
            do j=1,6
             if(i.eq.j)xjn_2d(i,j)=xjn_2d(i,j)+1.d0
            enddo
         enddo
      endif
      
      return
      end

!------------------------------------------
!    CALCULATES THE NORM OF THE RESIDUAL
      subroutine magn_vect(a,r,n)
      implicit real*8(a-h,o-z)
      dimension a(n)
      r=0.d0
      do i=1,n
         r=r+a(i)**2
      enddo
      r=dsqrt(r)
      return
      end
!-------------------------------------------
!    RU DECOMPOSITON          

      subroutine rudcmp(f,rd,ud) 
      implicit real*8 (a-h,o-z)                                        
               
      dimension rd(3,3),f(3,3),ud(3,3) 

!    
!      write(6,*)'entering ru'
      o3=1.0d0/3.0
      root3=dsqrt(3.d0)
      c11=f(1,1)*f(1,1)+f(2,1)*f(2,1)+f(3,1)*f(3,1)
      c12=f(1,1)*f(1,2)+f(2,1)*f(2,2)+f(3,1)*f(3,2)
      c13=f(1,1)*f(1,3)+f(2,1)*f(2,3)+f(3,1)*f(3,3)
      c23=f(1,2)*f(1,3)+f(2,2)*f(2,3)+f(3,2)*f(3,3)
      c22=f(1,2)*f(1,2)+f(2,2)*f(2,2)+f(3,2)*f(3,2)
      c33=f(1,3)*f(1,3)+f(2,3)*f(2,3)+f(3,3)*f(3,3)
      c1212=c12*c12
      c1313=c13*c13
      c2323=c23*c23
      c2313=c23*c13
      c1223=c12*c23
      c1213=c12*c13
      s11=c22*c33-c2323
      ui1=o3*(c11+c22+c33)
      ui2=s11+c11*c22+c33*c11-c1212-c1313
      ui3=c11*s11+c12*(c2313-c12*c33)  &
                 +c13*(c1223-c22*c13)
      ui1s=ui1*ui1
      q    =dsqrt(-dmin1(o3*ui2-ui1s,0.d0))                     
      r    =0.5*(ui3-ui1*ui2)+ui1*ui1s
      xmod =q*q*q
      scl1 =.5d0+dsign(.5d0,xmod-1.d-30)                          
      scl2 =.5d0+dsign(.5d0,xmod-dabs(r))                   
      scl0 =dmin1(scl1,scl2)                               
      scl1 =1.-scl0

      if(scl1.eq.0)then
        xmodscl1=xmod
      else
        xmodscl1=xmod+scl1
      endif

      sdetm=dacos(r/(xmodscl1))*o3

      q  =scl0*q
      ct3=q*dcos(sdetm)
      st3=q*root3*dsin(sdetm)
      sdetm=scl1*dsqrt(dmax1(0.0d0,r))
      aa=2.000*(ct3+sdetm)+ui1 
      bb=-ct3+st3-sdetm+ui1
      cc=-ct3-st3-sdetm+ui1                         
      xlamda1=dsqrt(dmax1(aa,0.d0))
      xlamda2=dsqrt(dmax1(bb,0.d0))
      xlamda3=dsqrt(dmax1(cc,0.d0))
      sdetm=xlamda1*xlamda2
      xli1=xlamda1+xlamda2+xlamda3
      xli2= sdetm+xlamda2*xlamda3+xlamda3*xlamda1
      xli3= sdetm*xlamda3/xli1
      s11= c11+xli3
      s22= c22+xli3
      s33= c33+xli3
      s12= c2313-c12*s33
      s13= c1223-s22*c13
      s23=-c2323+s22*s33
      sdetm=1./(xli1*(s11*s23+c12*s12+c13*s13))
      c11=c11+xli2
      c22=c22+xli2
      c33=c33+xli2
      si11=sdetm*s23
      si12=sdetm*s12
      si13=sdetm*s13
      si22=sdetm*( s11*s33-c1313)
      si23=sdetm*(-s11*c23+c1213)
      si33=sdetm*( s11*s22-c1212)
      s12=c12*si12
      s13=c13*si13
      s23=c23*si23
      ui11=c11*si11+s12+s13
      ui22=s12+c22*si22+s23
      ui33=s13+s23+c33*si33
      ui12=c11*si12+c12*si22+c13*si23
      ui13=c11*si13+c12*si23+c13*si33
      ui23=c12*si13+c22*si23+c23*si33
      rd(1,1)=f(1,1)*ui11+f(1,2)*ui12+f(1,3)*ui13
      rd(1,2)=f(1,1)*ui12+f(1,2)*ui22+f(1,3)*ui23
      rd(1,3)=f(1,1)*ui13+f(1,2)*ui23+f(1,3)*ui33
      rd(2,1)=f(2,1)*ui11+f(2,2)*ui12+f(2,3)*ui13
      rd(2,2)=f(2,1)*ui12+f(2,2)*ui22+f(2,3)*ui23
      rd(2,3)=f(2,1)*ui13+f(2,2)*ui23+f(2,3)*ui33
      rd(3,1)=f(3,1)*ui11+f(3,2)*ui12+f(3,3)*ui13
      rd(3,2)=f(3,1)*ui12+f(3,2)*ui22+f(3,3)*ui23
      rd(3,3)=f(3,1)*ui13+f(3,2)*ui23+f(3,3)*ui33 



      do i=1,3
         do j=1,3
            ud(i,j)=0.d0
            do k=1,3
               ud(i,j)=ud(i,j)+rd(k,i)*f(k,j)
            enddo
         enddo
      enddo
!      write(6,*)'exiting ru'
      return
      end


      subroutine stereograph(xl,xm,xn,x,y)
      implicit real*8 (a-h,o-z)
      x=dsign(1.d0,xn)*xl/(xn*dsign(1.d0,xn)+1.d0)
      y=dsign(1.d0,xn)*xm/(xn*dsign(1.d0,xn)+1.d0)
      return
      end

      integer function ngrain(nel,m,n)
      implicit real*8 (a-h,o-z)
      n3=(nel-1)/(m*n)**2+1
      nz=(n3-1)/m+1
      nel_2d=nel-(n3-1)*(m*n)**2
      nt1=nel_2d-1
      nt2=m*n
      n1=nt1-(nt1/nt2)*nt2+1
      nx=(n1-1)/m+1
      n2=(nel_2d-n1)/(m*n)+1
      ny=(n2-1)/m+1
      ngrain=n**2*(nz-1)+n*(ny-1)+nx
      return
      end

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
!    MATRIX MULTIPLICATION SUBROUTINE 

      subroutine calc_F_Thermal(F_therm,F_therm_inv,det_thrm,  &
                                thrm_C,temp,temp0,rot)
      implicit none
      real(8), intent(out) :: F_therm(3,3),F_therm_inv(3,3) ! thermal deformation tensor, and its inverse
      real(8), intent(out) :: det_thrm    ! det of thermal deform grad (volume change)
      real(8), intent(in) :: thrm_C(3)    ! thermal expansion coefficients along principal iCrystal directions
      real(8), intent(in) :: temp,temp0   ! current ambient, and stress free temperature
      real(8), intent(in) :: rot(3,3)     ! rotation tensor principal configuration to undeformed crystal configuration
      real(8) :: delTemp
      real(8) :: rot_inv(3,3)
      
      delTemp = temp-temp0
      
      rot_inv = TRANSPOSE(rot)
      
      det_thrm = exp((thrm_C(1)+thrm_C(2)+thrm_C(3))*delTemp)
      
      F_therm(:,:)=0.D0
      F_therm(1,1)=exp(thrm_C(1)*delTemp)
      F_therm(2,2)=exp(thrm_C(2)*delTemp)
      F_therm(3,3)=exp(thrm_C(3)*delTemp)
      F_therm(:,:) = MATMUL(MATMUL(rot,F_therm),rot_inv)

      F_therm_inv(:,:)=0.D0
      F_therm_inv(1,1)=1.d0/F_therm(1,1)
      F_therm_inv(2,2)=1.d0/F_therm(2,2)
      F_therm_inv(3,3)=1.d0/F_therm(3,3)
      F_therm_inv(:,:) = MATMUL(MATMUL(rot,F_therm_inv),rot_inv)
      
      return
      end
   
      subroutine matmul_flex(aa,bb,cc,n,m,l)
      
!      include 'aba_param.inc'
      implicit none

!   this subroutine returns martix [cc] = matrix [aa] * matrix [bb]
!   n=no. of rows of [aa]      =no. of rows of [cc]
!   m=no. of columns of [aa]   =no. of rows of [bb]
!   l=no. of columns of [bb]   =no. of columns of [cc]

      real(8), intent(in) :: aa(n,m),bb(m,l)
      real(8), intent(out) :: cc(n,l)
      integer, intent(in) :: n,m,l
      real(8) :: temp_out(n,l)
      integer :: i,j,k

      temp_out(:,:)=0.D0
      do i=1,n
         do j=1,l            
            do k=1,m
      temp_out(i,j)=temp_out(i,j)+aa(i,k)*bb(k,j)
            enddo
         enddo
      enddo
      cc(:,:)=temp_out(:,:)
      return 
      end

      FUNCTION TRAP_TERM(SIG_NUM)
      INTEGER TRAP_TERM
      INTEGER SIG_NUM
          write(*,*) 'trapped:',SIG_NUM
          TRAP_TERM = 1
      END

      
      SUBROUTINE getTensionCompressionState(tensComp,stress_1d,gmt,nslip)
      
      use options, only: tensOnly 
      
      implicit none
      
      real(8), intent(in) :: stress_1d(6)
      integer, intent(out):: tensComp(nslip)
      real(8), intent(in) :: gmt(3,nslip)
      integer, intent(in) :: nslip
      
      real(8) :: trPK, tpk_2d(3,3),eigvalue(3),eigvector(3,3),T_slipPlane
      integer :: nrot
      
      ! method 1: hydrostatic stress
      trPK = stress_1d(1)+stress_1d(2)+stress_1d(3)
      if(trPK.GE.0.D0) then
         tensComp = 1 ! tension
      else 
         tensComp = 2 ! compression
      endif
      
      ! overwrite tension/compression state
      if (tensOnly==1) tensComp = 1
      if (tensOnly==-1) tensComp = 2
      
!         method 2: principal stresses (as suggested the paper)
!         eigenvalue analysis
!         ------------------------------------------------------
!         tpk_2d(1,1) = stress_1d(1)
!         tpk_2d(2,2) = stress_1d(2)
!         tpk_2d(3,3) = stress_1d(3)
!         tpk_2d(1,2) = stress_1d(4)
!         tpk_2d(2,1) = stress_1d(4)
!         tpk_2d(1,3) = stress_1d(5)
!         tpk_2d(3,1) = stress_1d(5)
!         tpk_2d(2,3) = stress_1d(6)
!         tpk_2d(3,2) = stress_1d(6)
!         call jacobi(tpk_2d,3,3,eigvalue,eigvector,nrot)
!         if(MAX(abs(eigvalue(1:3))).ge.0) then
!            tensComp = 1 ! tension
!         else 
!            tensComp = 2 ! compression
!         endif
!
!         method 3:
!         determine tens/Comp state on individual slip planes
!         (better if integrated into/called from the implicit material solver)
!      do isys=1,nslip
!         T_slipPlane = &
!         gmt(1,isys)*(stress_1d(1)*gmt(1,isys)+stress_1d(4)*gmt(2,isys)+stress_1d(5)*gmt(3,isys)) + &
!         gmt(2,isys)*(stress_1d(4)*gmt(1,isys)+stress_1d(2)*gmt(2,isys)+stress_1d(6)*gmt(3,isys)) +
!         gmt(3,isys)*(stress_1d(5)*gmt(1,isys)+stress_1d(6)*gmt(2,isys)+stress_1d(3)*gmt(3,isys))
!         if (T_slipPlane>0.d0) then
!            tensComp(isys) = 1   ! tension
!         else
!            tensComp(isys) = 2   ! compression
!         endif
!      enddo
      
      END SUBROUTINE
      
      SUBROUTINE getSoftHardDirection(softHard,stress_1d,s_tan,nslip)
      
      use options, only: softOnly
      implicit none
      integer, intent(out):: softHard(nslip)
      real(8), intent(in) :: stress_1d(6)
      real(8), intent(in) :: s_tan(3,3,nslip)
      integer, intent(in) :: nslip
      
      real(8) :: tau, stress_2d(3,3), trPK
      integer :: isys,i,j
      
      call sigmat(stress_1d,stress_2d)
      
      trPK = stress_1d(1)+stress_1d(2)+stress_1d(3)
      if (trPK.gt.0.d0) then
         softHard(1:nslip) = 2
      else
         softHard(1:nslip) = 1
      endif
      return
      
      do isys=1,nslip
         tau = 0.d0
         do i=1,3
            do j=1,3
               tau = tau + stress_2d(i,j)*s_tan(i,j,isys)
            enddo
         enddo
         if(dsign(1.d0,tau) > 0.d0) then
            softHard(isys) = 1   ! SOFT direction
         else
            softHard(isys) = 2   ! HARD direction
         endif
      enddo
      
      if (softOnly==1) softHard = 1
      if (softOnly==-1) softHard = 2
      
      END SUBROUTINE
