      ! module task_Timing 
      ! tracking, analysis and reporting of CPU-time allocation of several tasks in a program
      ! usage: 
      ! N_TASKS : specify the predefined number of tasks
      ! taskIDs : chose constant variable names for tasks, and give indices starting from 1.
      ! taskNames: enter the names of tasks, padded with space to 30 characters
      !
      ! at the beginning of the program:
      ! CALL task_Timing_Initialize()
      !
      ! to mark the start of a task:
      ! CALL task_Start(taskID,MPI_WTIME())
      ! ... task ...
      !
      ! to mark the end of a task, and add to records:
      ! CALL task_Record(taskID,MPI_WTIME(),task_time,printTo)
      ! ... OUT: time passed during task is saved to task_time
      ! ... OUT: one line report is printed on screen if printTo = 0,
      !                                     on file if printTo = file handle number,
      !                                     no text output if < 0

      MODULE task_Timing

      integer, parameter :: N_TASKS = 12                     ! number of tasks
      
      integer, parameter :: task_PROGRAM = 1,        &      ! taskIDs
                            task_calcResidual = 2,   &
                            task_calcResidual_K = 3, &
                            task_LU_backsubs = 4,    &
                            task_LU_refactorize = 5, &
                            task_loadVC_comm = 6,    & 
                            task_loadVC_elst = 7,    &
                            task_lineSearch = 8,     &
                            task_BFGS_operation = 9, &
                            task_Broy_update = 10,   &
                            task_PostProcess = 11,   &
                            task_Export = 12

      integer, parameter :: LEN_TASKNAME = 30

      character(len=*), parameter :: &    ! task names
         taskNames = "PROGRAM                       &
                     &Residual Calculation          &
                     &Residual and K Calculation    &
                     &LU back-substitution          &
                     &LU refactorization            &
                     &loadVC communication for K    &
                     &loadVC - all ELST calls       &
                     &lineSearch                    &
                     &BFGS operations               &
                     &Broyden update                &
                     &Post processing               &
                     &Data Export                   "
                     
      integer :: taskCounter(N_TASKS)
      real(8) :: taskTimes(N_TASKS)
      real(8) :: task_sta(N_TASKS)
      real(8) :: task_end(N_TASKS)

      
      type testStruct
         integer :: var1
         integer :: var2
      end type testStruct
      
      contains 

         SUBROUTINE task_Timing_Initialize()
         implicit none
         
         taskCounter(:) = 0
         taskTimes(:) = 0.D0
         task_sta(:) = 0.D0
         task_end(:) = 0.D0
         
         END SUBROUTINE task_Timing_Initialize
      
         SUBROUTINE task_Start(taskID,wtime)
         implicit none
         integer, intent(in) :: taskID
         real(8), intent(in) :: wtime
         
         if (taskID.GT.N_TASKS) return

         task_sta(taskID) = wtime
         
         END SUBROUTINE task_Start
         
         SUBROUTINE task_Record(taskID,wtime,task_time,printTo)
         implicit none
         integer, intent(in) :: taskID
         real(8), intent(in) :: wtime
         real(8), intent(out) :: task_time
         integer, intent(in):: printTo
         
         integer :: name_sta, name_end
         
         if (taskID.GT.N_TASKS) return
         
         task_end(taskID) = wtime
         task_time = task_end(taskID) - task_sta(taskID)
         taskCounter(taskID) = taskCounter(taskID) + 1
         taskTimes(taskID) = taskTimes(taskID) + task_time

         ! print information
         name_sta = (taskID-1)*LEN_TASKNAME + 1
         name_end = (taskID)*LEN_TASKNAME
         if(printTo.GT.0) then
            write(printTo,*) 'Timer - ',taskNames(name_sta:name_end),task_time
         elseif(printTo.EQ.0) then
            write(*,*) 'Timer - ',taskNames(name_sta:name_end),task_time
         endif
         
         END SUBROUTINE task_Record   

         SUBROUTINE time_Read(taskID,taskName,nRecords,timeTotal,timeAvg)
         implicit none
         integer, intent(in) :: taskID
         character(len=*), intent(out) :: taskName
         integer, intent(out) :: nRecords
         real(8), intent(out) :: timeTotal, timeAvg
         
         integer :: name_sta, name_end
         
         if (taskID.GT.N_TASKS) return
         
         name_sta = (taskID-1)*LEN_TASKNAME + 1
         name_end = (taskID)*LEN_TASKNAME
         
         taskName = taskNames(name_sta:name_end)
         nRecords = taskCounter(taskID)
         timeTotal = taskTimes(taskID)
         timeAvg = taskTimes(taskID) / taskCounter(taskID)
         
         END SUBROUTINE time_Read   
         
         SUBROUTINE task_Report(printTo)
         implicit none
         integer, intent(in):: printTo
         
         integer :: taskID
         integer :: name_sta, name_end
         
         do taskID=1,N_TASKS
            if(taskCounter(taskID).EQ.0) cycle  !skip, if no records for this task
            ! print information
            name_sta = (taskID-1)*LEN_TASKNAME + 1
            name_end = (taskID)*LEN_TASKNAME
            if(printTo.GT.0) then
               write(printTo,*) taskNames(name_sta:name_end)
               write(printTo,*) '    # of calls:',taskCounter(taskID)
               write(printTo,*) '    total time:',taskTimes(taskID)
               write(printTo,*) '     avg. time:',taskTimes(taskID)/taskCounter(taskID)
            elseif(printTo.EQ.0) then
               write(*,*) taskNames(name_sta:name_end)
               write(*,*) '    # of calls:',taskCounter(taskID)
               write(*,*) '    total time:',taskTimes(taskID)
               write(*,*) '     avg. time:',taskTimes(taskID)/taskCounter(taskID)
            endif
         enddo
         
         
         END SUBROUTINE
      
      END MODULE
      