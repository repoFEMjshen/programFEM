      program yieldStress
      
      implicit none
      
      logical :: elasticModulusDetected,yieldStressDetected,loadDirectionDetected
      real(8) :: elasticModulus,plasticModulus,plasticStrainApprox
      integer :: loadDir
      real(8) :: totStress_prev(6),totStrain_prev(6)
      real(8) :: stressLoadDir, strainLoadDir, stressLoadDir_prev, strainLoadDir_prev
      character(len=128) :: strDataFileName,strCommandLine
      logical :: readUserProvidedFile, fileExists, printOutputToFile
      real(8), parameter :: requestedStress_atStrain = 0.035d0
      logical :: requestedStressDetected
      real(8) :: stressYield, strainYield, requestedStress
      integer :: err_read
      real(8) :: FINTME, totVol, totStress(6), totStrain(6),totCrZZ,totWP
      
      totStress = 0.d0
      totStrain = 0.d0
      requestedStress = 0.d0
      plasticModulus = 0.d0
      stressYield = 0.d0
      elasticModulusDetected = .false.
      yieldStressDetected = .false.
      loadDirectionDetected = .false.
      readUserProvidedFile = .false.
      printOutputToFile = .false.
      requestedStressDetected = .false.
      
      if (iargc() >= 1) then                                   ! User specified stress-strain data file
         CALL getarg(1, strCommandLine)
         strDataFileName = trim(strCommandLine)
         inquire(file=trim(strDataFileName),exist=fileExists)
         if (fileExists) readUserProvidedFile = .true.
      endif      
      if (.not.readUserProvidedFile) strDataFileName = 'stress_strain_RVEavg.out'

      if (iargc() >= 2) then                                   ! User specified output file
         CALL getarg(2, strCommandLine)
         if (trim(strCommandLine(1:7))=='-tofile') printOutputToFile = .true.
      endif      
      
      write(*,*) 'reading from file:',trim(strDataFileName)
      
      open(100,file=trim(strDataFileName))
      CALL readstr(100)
      
      err_read = 0
      
      totStress = 0.d0
      totStrain = 0.d0

      do while (err_read == 0) 
      
         totStress_prev = totStress
         totStrain_prev = totStrain
         read(100,*,iostat=err_read) FINTME, totVol, totStress(:), totStrain(:),totCrZZ,totWP
         
         if (.not.loadDirectionDetected) then
            loadDirectionDetected = .true.
            loadDir = maxloc(dabs(totStress(1:3)),1)
            write(*,*) '** detected loading direction: ', char(ichar('X')+loadDir-1)
         endif
         
         write(*,*) totStress(loadDir), totStrain(loadDir)
      
         if (err_read /= 0) exit
         

         strainLoadDir = totStrain(loadDir)
         stressLoadDir = totStress(loadDir)
         strainLoadDir_prev = totStrain_prev(loadDir)
         stressLoadDir_prev = totStress_prev(loadDir)
         if (strainLoadDir < 0.d0) then ! compression
            strainLoadDir = - strainLoadDir
            stressLoadDir = - stressLoadDir
            strainLoadDir_prev = - strainLoadDir_prev
            stressLoadDir_prev = - stressLoadDir_prev
         endif

         
         if (abs(stressLoadDir) > 50.d0 .and. (.not.elasticModulusDetected)) then

            elasticModulusDetected = .true.
            elasticModulus = stressLoadDir / strainLoadDir
            
            write(*,*) 'elasticModulus',elasticModulus

         endif


         if (elasticModulusDetected .and. .not.yieldStressDetected) then

            plasticStrainApprox = strainLoadDir - stressLoadDir/elasticModulus

            if (plasticStrainApprox > 0.002) then
            
            
               plasticModulus = (stressLoadDir - stressLoadDir_prev) / (strainLoadDir - strainLoadDir_prev)
               
               strainYield = 1.d0 / (elasticModulus - plasticModulus) * &
                  (-plasticModulus*strainLoadDir_prev + stressLoadDir_prev + 0.002*elasticModulus)
               stressYield = (strainYield - 0.002)*elasticModulus

!               stressYield = (elasticModulus*plasticModulus) / (elasticModulus + plasticModulus) * &
!                             (strainLoadDir_prev - stressLoadDir_prev/plasticModulus-0.002)
               
               write(*,*) '** detected yield stress'
               yieldStressDetected = .true.
            
            endif
            
         endif
         
         if ( .not. requestedStressDetected .and. &
                    strainLoadDir_prev < requestedStress_atStrain .and. &
                    strainLoadDir > requestedStress_atStrain ) then
                    
            requestedStress = stressLoadDir_prev + (stressLoadDir-stressLoadDir_prev) / &
                  (strainLoadDir-strainLoadDir_prev)*(requestedStress_atStrain-strainLoadDir_prev)
            write(*,*) '** detected stress at strain:',requestedStress_atStrain
            requestedStressDetected = .true.
         
         endif
         
         !if(yieldStressDetected) exit
         
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

      if (printOutputToFile) then
         open(101,file='YieldStress_'//trim(strDataFileName))
         write(101,*) stressYield, plasticModulus, requestedStress
         close(101)
         write(*,*) 'yield stress printed to file:','YieldStress_'//trim(strDataFileName)
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