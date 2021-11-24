      SUBROUTINE runUmatSubroutine(pathInput)

      use task_Timing
      use options
      use material_Ti6242
      use crystal_Titanium

      IMPLICIT NONE
      
      
      character(100), intent(in) :: pathInput
      character(100) :: path
      COMMON/pathVar/path

      logical, parameter :: Ti6242 = .TRUE.
      integer, parameter :: MAXPROPS = 400
      integer, parameter :: MAXSV = 400
      
      
      integer, parameter :: nPoints = 1
      integer, parameter :: NELX = 1
      integer, parameter :: elemID = 1
      integer :: grainIDelem(1)
      real(8) :: pointTexture(3), pointGrainSize(1)
      integer :: pointPhase(1)
      real(8), allocatable :: svarsOld(:),svarsNew(:)
      real(8) :: props(MAXPROPS)
      integer :: matID,nSvars,nprops
      logical :: importedMaterial
      
      logical :: elasticModulusDetected,yieldStressDetected
      real(8) :: elasticModulus,plasticModulus,plasticStrainApprox
      
      real(8) :: stressPK(6),stressCauchy(6),strainLog(6),stressCauchy_prev(6),strainLog_prev(6)
      real(8) :: stressZ, strainZ, stressZ_prev, strainZ_prev
      real(8) :: stressYield, strainYield
      real(8) :: R_stress(6)  ! stress residual for SIGxx=0 SIGyy=0
      real(8) :: maxR
      real(8) :: nye
      integer :: nIter, totalIter
      integer, parameter :: nIter_max = 200
      real(8) :: R_rxn
      real(8), parameter :: tolR = 0.01d0      ! 1% of the reaction force is fine
      real(8) :: DDSDDE(6,6), coeff
      integer, parameter :: NTENS = 6
      real(8) :: timeIncr,maxTimeIncr,timeIncrNew,timeLast
      real(8) :: timePrev, timeOld, timeNew, times(2)
      real(8) :: F_New(3,3), F_Old(3,3), E_New(6), E_Old(6)
      real(8) :: TEMP(2), TEMP_0
      integer :: umat_Err, ierr, nCutBack, nCutBackTotal
      integer, parameter :: maxCutBack = 5
      logical :: cutBack, forceBalance
      
      integer :: i,ii,j,jj,k,kk,iDOF, iDir
      integer :: INFO,error
      integer :: nStep
      integer :: iflag_over
      
      ! loading input variables
      integer, parameter :: MAX_HISTORY_N = 1000
      ! thermal loading information      
      integer :: temp_loading_type ! 1: temperature time history given
                                                ! 2: cyclic thermal loading around initial temp
                                                ! 3: cyclic thermal loading ABOVE initial temp
      real(8) :: temp_stress_free  ! initial/stress-free temperature
      real(8) :: temp_cyclic_amplitude   ! in Kelvins
      real(8) :: temp_cyclic_period      ! in seconds
      real(8) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
      integer :: temp_history_N                          ! for integration time steps

      ! info=1: temp loading imported, info=0: error, or no temp. loading
      
      integer(4) :: paddingTemp
      common/temperatureHistory/paddingTemp,temp_loading_type,temp_stress_free, &
            temp_cyclic_amplitude,temp_cyclic_period,temp_history_values, &
            temp_history_times,temp_history_N 
      
      ! loading input variables
      ! thermal loading information      
      integer :: load_type  ! 1: F-controlled
                            ! 2: eps_z-controlled (like CSR)
                            ! 3: S-controlled (like creep)
      
      integer :: deform_history_repeat    ! 1: repeat provided F-time history, 0: keep u at the last provided value
      real(8) :: deform_history_values(9,MAX_HISTORY_N) ! provide F-time points. 
      real(8) :: deform_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: deform_history_N    ! number of F-time points
      
      real(8) :: epsZ_cyclic_max  ! in MPa
      real(8) :: epsZ_cyclic_min  ! in MPa
      real(8) :: epsZ_cyclic_period ! in seconds
      integer :: epsZ_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: epsZ_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: epsZ_history_repeat
      real(8) :: epsZ_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: epsZ_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: epsZ_history_N    ! number of load-time points
      
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
      
      common/deformationHistory/load_type,deform_history_repeat, &
            deform_history_values,deform_history_times,deform_history_N, &
            epsZ_cyclic_max,epsZ_cyclic_min,epsZ_cyclic_period,epsZ_dwell_load, &
            epsZ_dwell_ramptime,epsZ_history_repeat,epsZ_history_values, &
            epsZ_history_times,epsZ_history_N,P_cyclic_max,P_cyclic_min, &
            P_cyclic_period,P_dwell_load,P_dwell_ramptime,P_history_repeat, &
            P_history_values,P_history_times,P_history_N
      
      
      ! blocks below transfer data between FEM_UEL and UMAT for debugging purposes only
      real(8) :: realData(5000)
      integer :: intData(5000)
      integer :: curElemUMAT
      COMMON/DEBUGINFO/realData,intData,curElemUMAT
      
      !modified by deniz: new FNAME size: 30
      CHARACTER(len=50) :: grainFileName
      
   ! export variables: total plastic slip
      real(8) :: tot_gamma(100)
      real(8) :: max_gamma,norm_tot_gamma
      
   ! other
      character(10) :: strNSTEP
      character(10) :: strNITER
      character(10) :: strPOSQ
      real(8) :: edgeLength
      real(8) :: xNodes(3,4), xNodes0(3,4)
      real(8) :: elemJacMat(3,3)
      real(8) :: elemJacInv(3,3)
      real(8) :: elemJac, elemVol, elemVol0, totVol
      real(8) :: totStress(6),totStrain(6),totFpZZ,totCrZZ,Ep_Equiv
      integer, parameter :: NGRAINAVG = 7
      real(8), allocatable :: grainAVG(:,:)   ! grain-averaged quantities are stored here
      real(8), allocatable :: grainRot(:,:)   ! grain-averaged quantities are stored here
      real(8) :: grainSchmidtBasalAll(6), grainSchmidtPrismAll(6)
      real(8) :: gst(3,nslip_max),gmt(3,nslip_max)
      real(8) :: R_cry(9), R_grain(9)
      integer :: totBndArea, grgrBndArea
      real(8) :: misorient, misorient_ave
      real(8) :: euler(3), rot(3,3), isliptmp
      integer :: n_HCP_planes
      parameter(n_HCP_planes=8)
      integer :: nslips(n_HCP_planes), ics(n_HCP_planes)
      integer :: isys, icrys ! icrys - a pointless variable
      real(8) :: hard(47)
      

      ! take the path provided by the user into the common block
      path = pathInput 
      
      ! read options
      call readOptions(trim(path)//'options.inp')
      
      call summarizeOptions()
         
      OPEN(981,FILE=trim(path)//'stress_strain_RVEavg.out')
      write(981,*) 'incr# time stress strain creepZZ'
      
      open(900,file=trim(path)//'debug.txt')
      open(901,file=trim(path)//'forceBalance.txt')
      open(902,file=trim(path)//'umatConverge.txt')
      open(904,file=trim(path)//'evol_g.txt')
      open(132,file=trim(path)//'time.out')
      open(504,file=trim(path)//'gamma.out')
      open(506,file=trim(path)//'fp.out')
      open(503,file=trim(path)//'F.out')
      open(507,file=trim(path)//'creep.out')
      
      if(Ti6242) nSvars = mat_nstate
      
      
      ! allocate state variables vector
      allocate(svarsOld(nSvars))
      allocate(svarsNew(nSvars))
      svarsOld(:) = 0.D0
      svarsNew(:) = 0.D0
      
      !read temperature
      ! deniz - read thermal loading info
      CALL readTemp(INFO)
                    
      !read deformation history
      CALL readDeformationHistory(timeIncr, maxTimeIncr, timeLast)

      !read texture and phase
      open(101,file=trim(path)//'pointMorphology.inp')
      read(101,*) pointTexture(1:3), pointPhase(1), pointGrainSize(1)
      read(101,*,iostat=error) nye
      if (error/=0) then
         nye=0.d0
      else
         write(*,*) 'assuming isotropic Nye: ',nye
      endif
      close(101)
   
      ! export Schmid factors for this orientation
      call euler_slip(pointTexture,rot)
      
      call init_crys_HCP() ! initialize HCP slip systems
      do isys=1,nslip_HCP
         gst(:,isys) = MATMUL(rot,cst_HCP(:,isys))
         gmt(:,isys) = MATMUL(rot,cmt_HCP(:,isys))
      enddo
      
      open(106,file='schmid.out')
      
      write(106,*) 'basal-a'
      do isys=1,3   
         write(106,'(F10.3)',advance='no') dabs(gst(3,isys)*gmt(3,isys))
      enddo      
      write(106,*) ''
      write(106,*) 'prism-a'
      do isys=4,6
         write(106,'(F10.3)',advance='no') dabs(gst(3,isys)*gmt(3,isys))
      enddo
      write(106,*) ''
      write(106,*) 'pyramidal-a'
      do isys=7,12
         write(106,'(F10.3)',advance='no') dabs(gst(3,isys)*gmt(3,isys))
      enddo
      write(106,*) ''
      write(106,*) '1st pyramidal{10-11}-c+a'
      do isys=13,24
         write(106,'(F10.3)',advance='no') dabs(gst(3,isys)*gmt(3,isys))
      enddo
      write(106,*) ''
      write(106,*) '2nd pyramidal{11-22}-c+a'
      do isys=25,30
         write(106,'(F10.3)',advance='no') dabs(gst(3,isys)*gmt(3,isys))
      enddo
      close(106)
      
      ! construct grainID-element map for a single point
      grainIDelem(1) = 1
      
      ! read material properties      
      write(*,*) 'Material: Ti-6242'
      call readMaterial(props,nprops,trim(path)//'props_Ti6242.dat',importedMaterial)
      if(.not.importedMaterial) then
         write(*,*) 'Can not import material parameters: props_Ti6242.dat corrupt/not found'
         STOP
      endif      
      ! initialize material parameters
      CALL initialize_Parameters(props,nprops,'Ti6242')
      ! initialize state variables
      CALL initialize_StateVariables(svarsOld,nSvars,elemID,pointGrainSize,pointTexture,pointPhase,nPoints,grainIDelem,NELX)
      
      ! set GND content/Nye tensor, if provided:
      svarsOld(357:365) = nye
      
      write(*,*) '# properties: ',nprops
      write(*,*) '# state variables: ',nSvars
!*************************************************************************

      OPEN(132,FILE=trim(path)//'time.out')

      ! initialize timer      
      CALL task_Timing_Initialize()


      ! *******************************************************************************!
		! the main loop over time steps
      ! *******************************************************************************!
      
      nStep= 1
      nCutBack = 0
      nCutBackTotal = 0
      totalIter = 0
      IFLAG_OVER=0
      
      ! initialize data block UMAT debugging
      curElemUMAT = elemID
      realData = 0.d0
      intData = 0
      
      
      ! initialize deformation gradients
      F_Old = 0.d0
      F_Old(1,1) = 1.d0
      F_Old(2,2) = 1.d0
      F_Old(3,3) = 1.d0
      F_New = 0.d0
      F_New(1,1) = 1.d0
      F_New(2,2) = 1.d0
      F_New(3,3) = 1.d0
      E_Old = 0.d0
      E_New = 0.d0
      timeOld = 0.d0
      
      stressCauchy = 0.d0
      strainLog = 0.d0
      elasticModulusDetected = .false.
      yieldStressDetected = .false.

      
		! IFLAG_OVER = 1 time integration is complete
      DO WHILE(IFLAG_OVER.EQ.0)
      
         cutBack = .false.
      
         ! update time
         timeNew = timeOld + timeIncr
      
         ! update temperature according to the prescribed thermal loading
         times(1) = timeOld
         times(2) = timeNew
         CALL tempInc(times(:), TEMP(:), TEMP_0)
                      
         ! update F or E
         F_new = F_Old
         E_new = E_Old
         CALL deformInc(timeNew,F_new,E_new)
      
         WRITE(*,*) nStep,'tNew=',timeNew,'dt=',timeIncr,'T=',TEMP(2)
         write(900,*) 'INCREMENT',nStep
         write(900,*) nStep,'tNew=',timeNew,'dt=',timeIncr,'T=',TEMP(2)
         write(901,*) 'INCREMENT',nStep
         write(901,*) nStep,'tNew=',timeNew,'dt=',timeIncr,'T=',TEMP(2)
         write(902,*) 'INCREMENT',nStep
         write(902,*) nStep,'tNew=',timeNew,'dt=',timeIncr,'T=',TEMP(2)
         
         if(load_type==1) then ! F-controlled
         
            ! reset UMAT interation counters
            intData = 0
            ! initialize state variables for time tau
            svarsNew = svarsOld
            ! F_Old: deform gradient at time t
            ! F_New: trial deform gradient at time tau
            CALL UMAT_CPFEM_THERM(stressPK,svarsNew,props,DDSDDE,NTENS,nSvars, & 
            nProps,timeIncrNew,F_Old,F_New,TEMP(:),TEMP_0,nStep,timeIncr,elemID,umat_Err) ! Deniz - thermal
            
            ! accept the stresses. this is F-controlled
            forceBalance = .true.
            
         elseif(load_type==2.or. &
                load_type==3) then     ! eps_Z-controlled -- need to solve implicit problem: Fzz=Fzz*, SIGxx=0, SIGyy=0, SIGxy=0, SIGxz=0, SIGyz=0
                  
            forceBalance = .false.
            nIter = 0
            
            
            
            do while(.not.forceBalance .and. nIter < nIter_max)
               
               ! reset UMAT interation counters
               intData = 0
               ! initialize state variables for time tau
               svarsNew = svarsOld
               CALL getFfromE(F_new,E_new)
               write(900,*) 'try E', E_new(:)
               CALL printMatrix(F_new,3,3,'try F',900)
               flush(900)
               write(902,*) 'old E:',E_Old(:)
               write(902,*) 'try E:',E_new(:)
               flush(902)
               
               CALL UMAT_CPFEM_THERM(stressPK,svarsNew,props,DDSDDE,NTENS,nSvars, & 
               nProps,timeIncrNew,F_Old,F_New,TEMP(:),TEMP_0,nStep,timeIncr,elemID,umat_Err) ! Deniz - thermal
               write(900,*) 'result S:', stressPK(:)
               if(DEBUG_UMAT) then
                  write(902,*) '# back-stress iterations:',intData(5000)
                  totalIter = totalIter + intData(5000)
                  do ii=1,intData(5000)
                     write(902,'(A,I5,A,6F10.1,A,6ES10.2,A,6F10.1)') '    ', intData(ii), &
                                  ' g_t', realData((ii-1)*6+1:(ii-1)*6+6), &      ! CRSS
                                  ' delt_gamma', realData((ii-1)*6+1001:(ii-1)*6+1006), &      ! delta gamma
                                  ' g_sat', realData((ii-1)*6+2001:(ii-1)*6+2006)   ! saturation stress
                  enddo
                  flush(902)
               endif
               
               if(umat_Err /= 0 .or. timeIncrNew < 1.d0) then
                  cutBack = .true.
                  nCutBackTotal = nCutBackTotal + 1
                  exit
               endif
               
               R_stress = stressPK
               !R_stress(1) = 0.d0   ! SIGxx=0
               !R_stress(2) = 0.d0   ! SIGyy=0
               R_rxn = R_stress(3)
               R_stress(3) = 0.d0   ! SIGzz= ?
               !R_stress(4) = 0.d0   ! SIGxy=0
               !R_stress(5) = 0.d0   ! SIGxz=0
               !R_stress(6) = 0.d0   ! SIGyz=0
               
               iDOF = 3
               coeff = DDSDDE(iDOF,iDOF)
               DDSDDE(:,iDOF) = 0.d0
               DDSDDE(iDOF,:) = 0.d0
               DDSDDE(iDOF,iDOF) = coeff
               
               DDSDDE(:,4) = DDSDDE(:,4) * 2.d0 ! strain definition -- mathematical vs. engineering strain
               DDSDDE(:,5) = DDSDDE(:,5) * 2.d0
               DDSDDE(:,6) = DDSDDE(:,6) * 2.d0
               
               DDSDDE(4,:) = DDSDDE(4,:) * 2.d0
               DDSDDE(5,:) = DDSDDE(5,:) * 2.d0
               DDSDDE(6,:) = DDSDDE(6,:) * 2.d0

!               iDOF = 5
!               coeff = DDSDDE(iDOF,iDOF)
!               DDSDDE(:,iDOF) = 0.d0
!               DDSDDE(iDOF,:) = 0.d0
!               DDSDDE(iDOF,iDOF) = coeff

!               iDOF = 6
!               coeff = DDSDDE(iDOF,iDOF)
!               DDSDDE(:,iDOF) = 0.d0
!               DDSDDE(iDOF,:) = 0.d0
!               DDSDDE(iDOF,iDOF) = coeff               
               
               maxR = maxval(dabs(R_stress(1:6)))
               write(901,*) dabs(R_stress(1:6)),' REF ',R_rxn
               if (maxR < tolR*abs(R_rxn)) forceBalance=.true.

               CALL solveAXb(DDSDDE,NTENS,R_stress,ierr)
               if (ierr /= 0) then
                  write(*,*) 'lapack error'
                  stop
               endif
               
               E_new = E_new - R_stress
               CALL getFfromE(F_new,E_new)
               
               nIter = nIter + 1
               
            enddo
                        
         endif
         
         if (cutBack) then

            nCutBack = nCutBack + 1
            timeIncr = timeIncr / 2.d0
            write(*,*) 'CUTBACK',nCutBack,timeIncr
            if (nCutBack < maxCutBack) then
               cycle
            else 
               write(*,*) 'too many cutbacks :('
               stop
            endif
         elseif(.not.forceBalance) then
            nCutBack = nCutBack + 1
            timeIncr = timeIncr / 2.d0
            write(*,*) 'forces not converge in ',nIter_max,' iterations: maxR',maxR,'F_rxn',R_rxn
            if (nCutBack < maxCutBack) then
               cycle
            else 
               write(*,*) 'too many cutbacks :('
               stop
            endif
         else
            ! reset the counter for consecutive cutbacks
            nCutBack = 0
            if(nIter < 5) timeIncr = timeIncr * 1.3d0
            if(nIter < 10) timeIncr = timeIncr * 1.1d0
            if (timeIncr > maxTimeIncr) timeIncr = maxTimeIncr
         endif
         
         if (DEBUG_UMAT) &         
            write(*,*) '# stress iterations', intData(4999),'# back-stress iterations',intData(5000)
         
         write(900,*) 'E - converged', E_new(:)
         CALL printMatrix(F_new,3,3,'F - converged',900)
         write(900,*) 'S - converged', stressPK(:)
         flush(900)
         
         if (DEBUG_UMAT) then
            write(904,*) ' CRSS', realData((intData(5000))*6+1:(intData(5000))*6+6), &      ! CRSS
                        ' g_sat', realData((intData(5000)-1)*6+2001:(intData(5000)-1)*6+2006), &   ! saturation stress
                        ' gamma_dot', (svarsNew(31:60)-svarsOld(31:60))/timeIncr ! delta gamma
            flush(904)
         endif
         
         ! EXPORT STRESS, STRAIN
         stressCauchy_prev = stressCauchy
         strainLog_prev = strainLog
         stressCauchy(:) = svarsNew(211:216) ! cauchy stress
         strainLog(:) = svarsNew(187:192) ! logarithmic strain
         
         totCrZZ = svarsNew(323)

         if (abs(stressCauchy(3)) > 50 .and. (.not.elasticModulusDetected)) then
         
            elasticModulusDetected = .true.
            elasticModulus = stressCauchy(3) / strainLog(3)
         
         endif
         

         if (elasticModulusDetected .and. .not.yieldStressDetected) then
         
            strainZ = strainLog(3)
            stressZ = stressCauchy(3)
            strainZ_prev = strainLog_prev(3)
            stressZ_prev = stressCauchy_prev(3)
            if (strainZ < 0.d0) then ! compression
               strainZ = - strainZ
               stressZ = - stressZ
               strainZ_prev = - strainZ_prev
               stressZ_prev = - stressZ_prev
            endif

            plasticStrainApprox = strainZ - stressZ/elasticModulus
         
            if (plasticStrainApprox > 0.002) then
            
               plasticModulus = (stressZ - stressZ_prev) / (strainZ - strainZ_prev)
               
               strainYield = 1.d0 / (elasticModulus - plasticModulus) * &
                  (-plasticModulus*strainZ_prev + stressZ_prev + 0.002*elasticModulus)
               stressYield = (strainYield - 0.002)*elasticModulus
            
               yieldStressDetected = .true.
            
            endif
            
         endif
         
         write(981,*) nStep, timeNew, stressCauchy(:), strainLog(:), totCrZZ
         flush(981)
         
         ! EXPORT TIME
         write(132,*) timeNew   ! time.out
         flush(132)
                  
         !GAMMA.out
         write(504,'(I0,F10.3,30ES10.2)') nStep,timeNew, svarsNew(31:60) !delta_gamma: 31--60
         flush(504)

         !FP
         write(506,'(I0,F10.3,9F12.7)')  nStep,timeNew, svarsNew(67:75) ! FP: 67--75   ! WARNING ! Fp - is stored ROW-MAJOR in STATE VARIABLES
         flush(506)

         !FP
         write(503,'(I0,F10.3,9F12.7)')  nStep,timeNew, F_New(1,1:3),F_New(2,1:3),F_New(3,1:3)
         flush(503)
         
         !Creep strain
         write(507,'(I0,F10.3,6F9.5)')  nStep,timeNew,svarsNew(321:326) ! G-L Creep Strain: 321-328
         flush(507)

                  
         ! ****** update displacement field solution ****** !
         ! update nodal positions,
         !     v_t <-- v_tau
         ! update variables for improving VTAU0
         F_Old          = F_New
         E_Old          = E_New
         timePrev       = timeOld
         timeOld        = timeNew

         nStep = nStep + 1

         ! *********** update state variables *********** !
         ! svarsOld: state variables from time t (converged) 
         ! svars2: trial state variables at time tau (not converged)
         ! svarsOld <- svars2
         svarsOld(:)=svarsNew(:)
         
!***********************************************************************************!
!        THE END OF CALCULATIONS FOR THIS TIME STEP. PROCEED TO THE NEXT.           !
!***********************************************************************************!
         if(timeNew > timeLast) IFLAG_OVER = 1

      ENDDO
      
701   CONTINUE ! program termination point

      if (yieldStressDetected) then
         write(*,*) 'yield stress:', stressYield
         write(*,*) 'elastic modulus:', elasticModulus
         write(*,*) 'plastic modulus at yield:', plasticModulus
      endif
      
      ! export yield stress to file
      OPEN(508,FILE=trim(path)//'yieldStress.out')
      write(508,*) stressYield
      flush(508)
      close(508)

      write(*,*) 'total stress iter ', totalIter
      write(*,*) 'total # cutBacks', nCutBackTotal
      write(*,*) 'total # increments', nStep

      
      ! End of calculation
     
      deallocate(svarsOld,svarsNew)   ! state variables array
      ! **************************************

      close(507)
      close(506)
      close(503)
      close(504)
      close(900)
      close(901)
      close(904)
      close(132)
      close(981)

      
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


      SUBROUTINE TMSTEP(timeLast,XIDT,XMXDT,XMINDT,NOSTP)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      DIMENSION FLAG(20)
      READ(LR,505) (FLAG(I),I=1,20)
      READ(LR,*) timeLast,XIDT,XMXDT,XMINDT,NOSTP
505   FORMAT(20A4)
      RETURN
      END SUBROUTINE

!*****************************************************************


! deniz -- read thermal loading information
      SUBROUTINE readTemp(info)
                          
      implicit none
      
      integer, intent(out) :: info 

      integer, parameter :: MAX_HISTORY_N = 1000
      
      integer :: temp_loading_type ! 1: temperature time history given
                                                ! 2: cyclic thermal loading around initial temp
                                                ! 3: cyclic thermal loading ABOVE initial temp
      real(8) :: temp_stress_free  ! initial/stress-free temperature
      real(8) :: temp_cyclic_amplitude   ! in Kelvins
      real(8) :: temp_cyclic_period      ! in seconds
      real(8) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
      integer :: temp_history_N                          ! for integration time steps

      ! info=1: temp loading imported, info=0: error, or no temp. loading
      
      integer(4) :: padding
      common/temperatureHistory/padding,temp_loading_type,temp_stress_free, &
            temp_cyclic_amplitude,temp_cyclic_period,temp_history_values, &
            temp_history_times,temp_history_N 
            
      character(100) :: path
      COMMON/pathVar/path
      
      ! locals
      integer :: iTemp
      logical :: fileExists
      CHARACTER(len=32) :: fileName
      
      fileName=trim(path)//'thermal_loading.inp'
      
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

!************************************************************
      SUBROUTINE tempInc(TIME, TEMP, TEMP_0)
      implicit none

      real(8),parameter :: PI=3.14159265359D0
      
      ! INPUTS
      real(8), intent(in) :: TIME(2)
      ! OUTPUTS
      real(8), intent(out) :: TEMP(2), TEMP_0

      integer, parameter :: MAX_HISTORY_N = 1000
      
      integer :: temp_loading_type ! 1: temperature time history given
                                               ! 2: cyclic thermal loading around initial temp
                                               ! 3: cyclic thermal loading ABOVE initial temp
      real(8) :: temp_stress_free  ! initial/stress-free temperature
      real(8) :: temp_cyclic_amplitude   ! in Kelvins
      real(8) :: temp_cyclic_period      ! in seconds
      real(8) :: temp_history_values(MAX_HISTORY_N) ! provide temp-time points. 
      real(8) :: temp_history_times(MAX_HISTORY_N)  ! program will interpolate
      integer :: temp_history_N                     ! for integration time steps
      
      integer*4 :: padding
      common/temperatureHistory/padding,temp_loading_type,temp_stress_free, &
            temp_cyclic_amplitude,temp_cyclic_period,temp_history_values, &
            temp_history_times,temp_history_N   
            
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
      
      function factorial(n) result(n_fact)
         integer, intent(in) :: n
         integer:: n_fact
         n_fact = 1
         do i=1,n
            n_fact=n_fact*i
         enddo
      end function factorial

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
      
      SUBROUTINE deformInc(timeNew, F_new, E_new)
      implicit none
      
      integer, parameter :: MAX_HISTORY_N = 1000
      
      real(8), intent(in)  :: timeNew
      real(8), intent(out) :: F_new(3,3)
      real(8), intent(out) :: E_new(6)
            
      integer :: load_type ! 1: F-controlled
                           ! 2: eps_z-controlled history
                           ! 3: eps_z-controlled cyclic
                           ! 4: S-controlled history
                           ! 5: S-controlled cyclic
      
      integer :: deform_history_repeat    ! 1: repeat provided F-time history, 0: keep u at the last provided value
      real(8) :: deform_history_values(9,MAX_HISTORY_N) ! provide F-time points. 
      real(8) :: deform_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: deform_history_N    ! number of F-time points
      
      real(8) :: epsZ_cyclic_max  ! in MPa
      real(8) :: epsZ_cyclic_min  ! in MPa
      real(8) :: epsZ_cyclic_period ! in seconds
      integer :: epsZ_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: epsZ_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: epsZ_history_repeat
      real(8) :: epsZ_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: epsZ_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: epsZ_history_N    ! number of load-time points
      
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
      
      common/deformationHistory/load_type,deform_history_repeat, &
            deform_history_values,deform_history_times,deform_history_N, &
            epsZ_cyclic_max,epsZ_cyclic_min,epsZ_cyclic_period,epsZ_dwell_load, &
            epsZ_dwell_ramptime,epsZ_history_repeat,epsZ_history_values, &
            epsZ_history_times,epsZ_history_N,P_cyclic_max,P_cyclic_min, &
            P_cyclic_period,P_dwell_load,P_dwell_ramptime,P_history_repeat, &
            P_history_values,P_history_times,P_history_N

      !locals
      integer :: iTime, t_or_tau
      real(8) :: w1,w2
      logical :: found
      real(8) :: epsZ_cyclic_amp, epsZ_cyclic_mean
      real(8) :: t_cyc, tramp_cyc
      real(8),parameter :: PI=3.14159265359D0

      IF(load_type.EQ.1)THEN       ! 1: F-controlled
      
         t_cyc = timeNew
         if(deform_history_repeat.EQ.1) then
            t_cyc = DMOD(timeNew,deform_history_times(deform_history_N))
         endif
         
         found=.FALSE.
         do iTime=1,deform_history_N
            if(t_cyc.LT.deform_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               F_new(1:3,1)=deform_history_values(1:3,1)
               F_new(1:3,2)=deform_history_values(4:6,2)
               F_new(1:3,3)=deform_history_values(7:9,3)
            else
               w1 = (deform_history_times(iTime)-t_cyc) / & 
            (deform_history_times(iTime)-deform_history_times(iTime-1))
               w2 = (t_cyc-deform_history_times(iTime-1)) / & 
            (deform_history_times(iTime)-deform_history_times(iTime-1))
               F_new(1:3,1) = w1*deform_history_values(1:3,iTime-1) & 
                         + w2*deform_history_values(1:3,iTime)
               F_new(1:3,2) = w1*deform_history_values(4:6,iTime-1) & 
                         + w2*deform_history_values(4:6,iTime)
               F_new(1:3,3) = w1*deform_history_values(7:9,iTime-1) & 
                         + w2*deform_history_values(7:9,iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            F_new(1:3,1)=deform_history_values(1:3,deform_history_N)
            F_new(1:3,2)=deform_history_values(4:6,deform_history_N)
            F_new(1:3,3)=deform_history_values(7:9,deform_history_N)
         endif
                     
      ELSEIF(load_type==2)THEN   ! 2: eps_z-controlled history
      
         t_cyc = timeNew
         if(epsZ_history_repeat.EQ.1) then
            t_cyc = DMOD(timeNew,epsZ_history_times(epsZ_history_N))
         endif
         
         found=.FALSE.
         do iTime=1,epsZ_history_N
            if(t_cyc.LT.epsZ_history_times(iTime))then
               found=.TRUE.
               exit
            endif
         enddo
         if(found)then
            if(iTime.EQ.1)then
               E_new(3) = epsZ_history_values(1)
            else
               w1 = (epsZ_history_times(iTime)-t_cyc) / & 
            (epsZ_history_times(iTime)-epsZ_history_times(iTime-1))
               w2 = (t_cyc-epsZ_history_times(iTime-1)) / & 
            (epsZ_history_times(iTime)-epsZ_history_times(iTime-1))
               E_new(3) = w1*epsZ_history_values(iTime-1) & 
                    + w2*epsZ_history_values(iTime)
            endif 
         else
            ! out of provided data range. assign the last provided value
            E_new(3)=epsZ_history_values(epsZ_history_N)
         endif      
         
      ELSEIF(load_type==3)THEN   ! 3: eps_z-controlled cyclic

         epsZ_cyclic_amp = 0.5D0*(epsZ_cyclic_max-epsZ_cyclic_min)
         epsZ_cyclic_mean = 0.5D0*(epsZ_cyclic_max+epsZ_cyclic_min)
         tramp_cyc = epsZ_cyclic_mean /  & 
                     (epsZ_cyclic_amp*2*PI/epsZ_cyclic_period)
  
         if (timeNew.LT.tramp_cyc) then
            E_new(3) = timeNew/tramp_cyc*epsZ_cyclic_mean
         else
            t_cyc = timeNew - tramp_cyc
            E_new(3) = epsZ_cyclic_mean +  & 
                     epsZ_cyclic_amp*SIN(2*PI*t_cyc/epsZ_cyclic_period)
         endif
         
      else
      
      
         write(*,*) 'called deformInc for unsupported load type:',load_type
         
      endif


         
      END SUBROUTINE
      
      SUBROUTINE readDeformationHistory(timeIncr, maxTimeIncr, timeLast)
      
      implicit none
                     
      integer, parameter :: MAX_HISTORY_N = 1000
      
      real(8), intent(out) :: timeIncr, maxTimeIncr, timeLast
      
      integer :: load_type ! 1: F-controlled
                           ! 2: eps_z-controlled history
                           ! 3: eps_z-controlled cyclic
                           ! 4: S-controlled history
                           ! 5: S-controlled cyclic
      
      integer :: deform_history_repeat    ! 1: repeat provided F-time history, 0: keep u at the last provided value
      real(8) :: deform_history_values(9,MAX_HISTORY_N) ! provide F-time points. 
      real(8) :: deform_history_times(MAX_HISTORY_N)    ! program will interpolate at integration time steps
      integer :: deform_history_N    ! number of F-time points
      
      real(8) :: epsZ_cyclic_max  ! in MPa
      real(8) :: epsZ_cyclic_min  ! in MPa
      real(8) :: epsZ_cyclic_period ! in seconds
      integer :: epsZ_dwell_load    ! 1: dwell type cycles, 0: sinusoidal cycles
      real(8) :: epsZ_dwell_ramptime! if dwell type cycles, indicates ramp time of dwell loading
      integer :: epsZ_history_repeat
      real(8) :: epsZ_history_values(MAX_HISTORY_N) ! provide load-time points. 
      real(8) :: epsZ_history_times(MAX_HISTORY_N)  ! program will interpolate
                                                 ! at integration time steps
      integer :: epsZ_history_N    ! number of load-time points
      
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
            
      common/deformationHistory/load_type,deform_history_repeat, &
            deform_history_values,deform_history_times,deform_history_N, &
            epsZ_cyclic_max,epsZ_cyclic_min,epsZ_cyclic_period,epsZ_dwell_load, &
            epsZ_dwell_ramptime,epsZ_history_repeat,epsZ_history_values, &
            epsZ_history_times,epsZ_history_N,P_cyclic_max,P_cyclic_min, &
            P_cyclic_period,P_dwell_load,P_dwell_ramptime,P_history_repeat, &
            P_history_values,P_history_times,P_history_N
            
      character(100) :: path
      COMMON/pathVar/path
            
      integer :: I
            
      OPEN(201,FILE=trim(path)//'deformationHistory.inp')
      
      CALL READSTR(201)
      
      read(201,*) timeIncr, maxTimeIncr, timeLast
      
      CALL READSTR(201)

      READ(201,*) load_type ! 1: F-controlled
                            ! 2: eps_z-controlled history
                            ! 3: eps_z-controlled cyclic
                            ! 4: S-controlled history
                            ! 5: S-controlled cyclic
                            
      IF(load_type.EQ.1)THEN       ! 1: F-controlled

         CALL READSTR(201)
         READ(201,*) deform_history_N, deform_history_repeat
         
          if(deform_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for strain_history to:', & 
                                                      deform_history_N
            STOP
         endif
         
         deform_history_times(:) = 0.0
         deform_history_values(:,:) = 0.0
         do I=1,deform_history_N
            READ(201,*) deform_history_times(I), deform_history_values(:,I)
         enddo
         
      ELSEIF(load_type==2)THEN   ! 2: eps_z-controlled history

         CALL READSTR(201)
         READ(201,*) epsZ_history_N, epsZ_history_repeat
      
          if(epsZ_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for strain_history to:', & 
                                                      epsZ_history_N
            STOP
         endif
         
         epsZ_history_times(:) = 0.0
         epsZ_history_values(:) = 0.0
         do I=1,epsZ_history_N
            READ(201,*) epsZ_history_times(I), epsZ_history_values(I)
         enddo
         
      ELSEIF(load_type==3)THEN   ! 3: eps_Z-controlled cyclic
      
         CALL READSTR(201)
         READ(201,*) epsZ_cyclic_max, epsZ_cyclic_min, & 
                     epsZ_cyclic_period, epsZ_dwell_load, &
                     epsZ_dwell_ramptime
      
      ELSEIF(load_type==4)THEN   ! 4: S-controlled history
      
         CALL READSTR(201)
         READ(201,*) P_history_N, P_history_repeat

          if(P_history_N.GT.MAX_HISTORY_N)then
            write(*,*) 'Increase MAX_HISTORY_N for strain_history to:', & 
                                                      P_history_N
            STOP
         endif
         
         P_history_times(:) = 0.0
         P_history_values(:) = 0.0
         do I=1,P_history_N
            READ(201,*) P_history_times(I), P_history_values(I)
         enddo
         
      ELSEIF(load_type==5)THEN   ! 5: S-controlled cyclic

         CALL READSTR(201)
         READ(201,*) P_cyclic_max, P_cyclic_min, & 
                     P_cyclic_period, P_dwell_load, &
                     P_dwell_ramptime
         
      ENDIF      

      CLOSE(201)
      
      END SUBROUTINE
      
                  
      SUBROUTINE getFfromE(F,E)
      implicit none
      real(8), intent(out) :: F(3,3)
      real(8), intent(in) :: E(6)
      
      F(1,1) = dsqrt(2.d0*E(1) + 1.d0)
      F(1,2) = 2.d0*E(4) / F(1,1)
      F(1,3) = 2.d0*E(5) / F(1,1)
      F(2,2) = dsqrt(2.d0*E(2)+1.d0-F(1,2)**2.d0)
      F(2,3) = (2.d0*E(6) - F(1,2)*F(1,3))/F(2,2)
      F(3,3) = dsqrt(2.d0*E(3) + 1.d0 - F(1,3)**2.d0 - F(2,3)**2.d0)
      
      END SUBROUTINE
      