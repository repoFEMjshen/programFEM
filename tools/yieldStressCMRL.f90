      program yieldStress
      
      implicit none
      
      logical :: elasticModulusDetected,yieldStressDetected,loadDirectionDetected
      real(8) :: elasticModulus,plasticModulus,plasticStrainApprox
      integer :: loadDir
      real(8) :: stressYield, strainYield, strainLoadDir, stressLoadDir, strainLoadDir_prev, stressLoadDir_prev
      integer :: err_read
      real(8) :: FINTME, sevolData(18), stress(3), strain(3), effPlasticStrain
      real(8) :: stress_prev(3), strain_prev(3)
      
      strain = 0.d0
      stress = 0.d0
      strainLoadDir = 0.d0
      stressLoadDir = 0.d0
      elasticModulusDetected = .false.
      yieldStressDetected = .false.
      loadDirectionDetected = .false.
      
      open(100,file='sevol.out')
      
      err_read = 0

      do while (err_read == 0) 
      
         read(100,*,iostat=err_read) FINTME, sevolData(:) ,effPlasticStrain
         
         strain(1) = sevolData(1)
         strain(2) = sevolData(4)
         strain(3) = sevolData(7)

         stress(1) = sevolData(3)
         stress(2) = sevolData(6)
         stress(3) = sevolData(9)
         
         if (.not.loadDirectionDetected) then
            loadDirectionDetected = .true.
            loadDir = maxloc(dabs(stress(1:3)),1)
            write(*,*) '** detected loading direction: ', char(ichar('X')+loadDir-1)
         endif
         
         strainLoadDir_prev = strainLoadDir
         stressLoadDir_prev = stressLoadDir
         
         strainLoadDir = strain(loadDir)
         stressLoadDir = stress(loadDir)

         
         write(*,*) strainLoadDir, stressLoadDir
      
         if (err_read /= 0) exit
         
         if (abs(stress(loadDir)) > 50.d0 .and. (.not.elasticModulusDetected)) then

            elasticModulusDetected = .true.
            elasticModulus = stress(loadDir) / strain(loadDir)
            
            write(*,*) '** detected elasticModulus: ',elasticModulus

         endif


         if (elasticModulusDetected .and. .not.yieldStressDetected) then
         
            if (strainLoadDir < 0.d0) then ! compression
               strainLoadDir = - strainLoadDir
               stressLoadDir = - stressLoadDir
               strainLoadDir_prev = - strainLoadDir_prev
               stressLoadDir_prev = - stressLoadDir_prev
            endif

            plasticStrainApprox = strainLoadDir - stressLoadDir/elasticModulus

            if (plasticStrainApprox > 0.002) then
            
            
               plasticModulus = (stressLoadDir - stressLoadDir_prev) / (strainLoadDir - strainLoadDir_prev)
               
               strainYield = 1.d0 / (elasticModulus - plasticModulus) * &
                  (-plasticModulus*strainLoadDir_prev + stressLoadDir_prev + 0.002*elasticModulus)
               stressYield = (strainYield - 0.002)*elasticModulus

!               stressYield = (elasticModulus*plasticModulus) / (elasticModulus + plasticModulus) * &
!                             (strainLoadDir_prev - stressLoadDir_prev/plasticModulus-0.002)
               
            
               yieldStressDetected = .true.
            
            endif
            
         endif
         
         if(yieldStressDetected) exit
         
      enddo
      
      close(100)
      

      if (elasticModulusDetected) then
         write(*,*) 'elastic modulus:'
         write(*,*) elasticModulus
      else
         write(*,*) 'data not sufficient to determine the elastic modulus'
      endif

      if (yieldStressDetected) then
         
         write(*,*) '0.2% yield Strain:'
         write(*,*) strainYield
         
         write(*,*) '0.2% yield Stress:'
         write(*,*) stressYield
               
      else
         write(*,*) 'data not sufficient to determine 0.2% yield stress'
      endif
      
      end program
      
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