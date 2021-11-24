      PROGRAM calcEqvPlasticStrainMacroscopic

      IMPLICIT NONE

      real(8) :: EpVM, Ep(6), rateEpVM, Dp(6), prevTime, curTime, dt
      
   ! some dummy variables

      logical :: fileExists
      character(len=255) :: fileDp,fileEpVM

      integer :: nStepsExported
      integer :: error
      
      fileDp = 'dp_RVEavg.out'
      fileEpVM = 'eqvPlStrain_RVEavg.out'
      
      inquire(file=trim(fileDp),exist=fileExists)
      if (.not. fileExists) then
         write(*,*) 'cannot find the input file:',trim(fileDp)
         stop -1
      endif
      
      open(104,file=trim(fileDp))     ! input
      open(105,file=trim(fileEpVM))   ! output
      
      call READSTR(104) ! skip header
      write(105,*) 'time EpVM Ep(6)'
      
      EpVM = 0.d0
      Ep = 0.d0
      prevTime = 0.d0
      error = 0
      nStepsExported = 0
      do while (error==0)
 
         ! import Dp for this time step
         error = 0
         read(104,*,iostat=error) curTime, Dp
         if(error /= 0) exit ! EOF
         
         ! time step
         dt = curTime - prevTime
         
         ! calculate VM plastic strain rate
         CALL getVonMisesStrainVoigt_MathStrain(rateEpVM,Dp)
         Ep = Ep + Dp*dt
         EpVM = EpVM + rateEpVM*dt
         
         ! export total effective plastic strain
         write(105,*) curTime, EpVM, Ep(1:6)
         nStepsExported = nStepsExported + 1
         
         prevTime = curTime
      enddo
      
      ! close files
      close(104)
      close(105)
      
      write(*,*) 'exported total effective plastic strain to file: ',trim(fileEpVM)
      write(*,*) '# time steps:',nStepsExported
      write(*,*) 'last time step:',curTime
      
      
      END


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
      END
